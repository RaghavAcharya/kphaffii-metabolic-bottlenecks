function model = addMisfoldingReactions(model, proteinInfo, rxnName, startMet)
% ADDMISFOLDINGREACTIONS Adds misfolding, refolding, ubiquitination,
% and degradation reactions for a translated protein with REALISTIC stoichiometry.
%
% INPUTS:
%   model       - GEM structure
%   proteinInfo - struct with .throughER (1 = ER-resident), .aaSequence, and optional .addRefoldingRxn
%   rxnName     - base gene/protein ID (e.g., GQ67_XXXX)
%   startMet    - starting metabolite ID for mature protein
%
% OUTPUT:
%   model       - updated GEM

%% Setup compartment and ID formats
if proteinInfo.throughER == 1
    comp = 'r';
else
    comp = 'c';
end

% Internal metabolite IDs (model.mets)
matureMetID      = sprintf('%s_mature_%s', rxnName, comp);
misfoldedMetID   = sprintf('%s_misfolded_%s', rxnName, comp);
ubiquitinatedMet = sprintf('%s_ub_%s', rxnName, comp);
ubMetID          = sprintf('GQ67_05145_mature_c');

% COBRA-style reaction IDs (for addRxns)
matureRxnID      = startMet;
misfoldedRxnID   = sprintf('%s_misfolded[%s]', rxnName, comp);
ubiquitinatedRxn = sprintf('%s_ub[%s]', rxnName, comp);
ubRxnID          = sprintf('GQ67_05145_mature[c]');

% Choose chaperone
if comp == 'r'
    chaperone = 'GQ67_00168_mature[r]';  % KAR2
else
    chaperone = 'GQ67_03278_mature[c]';  % SSA4
end

% Get protein length for realistic energy costs
proteinLength = length(proteinInfo.aaSequence);

%% 1. Misfolding Reaction (unchanged - this is reasonable)
rxnStruct1.rxns         = {sprintf('%s_misfold', rxnName)};
rxnStruct1.rxnNames     = rxnStruct1.rxns;
rxnStruct1.mets         = {{matureRxnID, misfoldedRxnID}};
rxnStruct1.stoichCoeffs = {[-1, 1]};
rxnStruct1.catBy        = {''};
rxnStruct1.eccodes      = {''};
rxnStruct1.lb           = [0];

metStruct1.mets         = {misfoldedMetID};
metStruct1.metNames     = {sprintf('%s_misfolded', rxnName)};
metStruct1.compartments = {comp};

model = addCatalyzedRxns(model, metStruct1, rxnStruct1, rxnName);

%% 2. Refolding Reaction (optional - controlled by proteinInfo.addRefoldingRxn)
% Check if refolding reactions should be added
addRefolding = false;
if isfield(proteinInfo, 'addRefoldingRxn') && proteinInfo.addRefoldingRxn
    addRefolding = true;
end

if addRefolding
    % Realistic cost: ~1 ATP per 50 amino acids for chaperone-assisted folding
    refold_atp_cost = max(1, ceil(proteinLength / 50));

    rxnStruct2.rxns     = {sprintf('%s_refold', rxnName)};
    rxnStruct2.rxnNames = rxnStruct2.rxns;
    rxnStruct2.mets     = {
        {misfoldedRxnID, chaperone, 'ATP[c]', 'H2O[c]', ...
         startMet, 'ADP[c]', 'Phosphate[c]', 'H+[c]'}
    };
    rxnStruct2.lb       = [0];
    rxnStruct2.stoichCoeffs = {
        [-1, -1, -refold_atp_cost, -refold_atp_cost, 1, refold_atp_cost, refold_atp_cost, refold_atp_cost]
    };
    rxnStruct2.catBy    = {''};
    rxnStruct2.eccodes  = {'3.6.4.10'};

    model = addCatalyzedRxns(model, struct('mets', {}, 'metNames', {}, 'compartments', {}), rxnStruct2, rxnName);
else
    % Skip refolding reaction - misfolded proteins can only be degraded
    % This prevents the model from "cheating" by refolding proteins to avoid degradation costs
end

%% 3. Ubiquitination (more realistic cost)
% Each ubiquitin attachment requires ATP, and larger proteins need more ubiquitin
% Use 1 ubiquitin per ~75 amino acids (minimum 4 for proteasomal recognition)
ubiquitin_required = max(4, ceil(proteinLength / 75));

rxnStruct3.rxns     = {sprintf('%s_ubiquitin', rxnName)};
rxnStruct3.rxnNames = rxnStruct3.rxns;
rxnStruct3.mets     = {
    {misfoldedRxnID, ubRxnID, 'ATP[c]', 'H2O[c]', ...
     ubiquitinatedRxn, 'ADP[c]', 'Phosphate[c]', 'H+[c]'}
};
rxnStruct3.lb       = [0];
rxnStruct3.stoichCoeffs = {
    [-1, -ubiquitin_required, -ubiquitin_required, -ubiquitin_required, ...
     1, ubiquitin_required, ubiquitin_required, ubiquitin_required]
};
rxnStruct3.catBy   = {''};
rxnStruct3.eccodes = {'2.3.2.27'};

metStruct3.mets         = {ubiquitinatedMet};
metStruct3.metNames     = {sprintf('%s_ub', rxnName)};
metStruct3.compartments = {comp};

model = addCatalyzedRxns(model, metStruct3, rxnStruct3, rxnName);

%% 4. CORRECTED Proteasomal Degradation
% Much more realistic stoichiometry:
% - ATP cost scales with protein size (literature: ~1 ATP per 4-6 amino acids)
% - Only partial amino acid recovery (some are lost/damaged during degradation)
% - Include realistic water requirements

% Parse sequence and count amino acids
aaSeq = proteinInfo.aaSequence;
aaList = 'ARNDCEQGHILKMFPSTWYV';
aaIDs = {
    'L-Alanine', 'L-Arginine', 'L-Asparagine', 'L-Aspartate', 'L-Cysteine', ...
    'L-Glutamate', 'L-Glutamine', 'Glycine', 'L-Histidine', 'L-Isoleucine', ...
    'L-Leucine', 'L-Lysine', 'L-Methionine', 'L-Phenylalanine', ...
    'L-Proline', 'L-Serine', 'L-Threonine', 'L-Tryptophan', ...
    'L-Tyrosine', 'L-Valine'
};

% Count residues
aaCounts = zeros(1, 20);
for i = 1:20
    aaCounts(i) = count(aaSeq, aaList(i));
end

% REALISTIC degradation costs:
% - ATP cost: 1 ATP per 5 amino acids (literature range: 4-6 aa per ATP)
degradation_atp_cost = ceil(proteinLength / 5);

% - Water cost: 1 H2O per peptide bond broken (n-1 peptide bonds)
water_cost = proteinLength - 1;

% - Amino acid recovery efficiency: 75% (some amino acids are damaged/lost)
recovery_efficiency = 0.75;
recoveredAaCounts = floor(aaCounts * recovery_efficiency);

% Generate COBRA-style metabolite IDs for recovered amino acids
aaMets = cellfun(@(id) sprintf('%s[c]', id), aaIDs, 'UniformOutput', false);

% Build reaction stoichiometry
inputs = {ubiquitinatedRxn, 'ATP[c]', 'H2O[c]'};
input_coeffs = [-1, -degradation_atp_cost, -water_cost];

outputs = {ubRxnID, 'ADP[c]', 'Phosphate[c]', 'H+[c]'};
output_coeffs = [ubiquitin_required, degradation_atp_cost, degradation_atp_cost, degradation_atp_cost];

% Add recovered amino acids (only if count > 0)
for i = 1:20
    if recoveredAaCounts(i) > 0
        outputs{end+1} = aaMets{i};
        output_coeffs(end+1) = recoveredAaCounts(i);
    end
end

% Combine all metabolites and coefficients
all_mets = [inputs, outputs];
all_coeffs = [input_coeffs, output_coeffs];

rxnStruct4.rxns         = {sprintf('%s_degrade', rxnName)};
rxnStruct4.lb           = [0];
rxnStruct4.rxnNames     = rxnStruct4.rxns;
rxnStruct4.mets         = {all_mets};
rxnStruct4.stoichCoeffs = {all_coeffs};
rxnStruct4.catBy        = {''};
rxnStruct4.eccodes      = {'3.4.25.1'};

model = addCatalyzedRxns(model, struct('mets', {}, 'metNames', {}, 'compartments', {}), rxnStruct4, rxnName);

% Print summary of what was added (optional debug output)
% if addRefolding
%     fprintf('Added misfolding reactions for %s (%d aa) WITH refolding\n', rxnName, proteinLength);
% else
%     fprintf('Added misfolding reactions for %s (%d aa) WITHOUT refolding\n', rxnName, proteinLength);
% end

end