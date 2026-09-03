function [model] = addTranslationReaction(model, aaSequence, dnaSequence, rxnName)
% ADDTRANSLATIONREACTION Adds a realistic translation reaction with proper stoichiometry
% 
% INPUTS:
%   model       - GEM structure
%   aaSequence  - amino acid sequence (single letter codes)
%   dnaSequence - DNA sequence (for validation, optional)
%   rxnName     - gene name (e.g., GQ67_XXXX)
%
% OUTPUT:
%   model       - updated GEM with translation reaction

% Define amino acid mapping
aaList = 'ARNDCEQGHILKMFPSTWYV';
aaIDs = {
    'L-Alanine', 'L-Arginine', 'L-Asparagine', 'L-Aspartate', 'L-Cysteine', ...
    'L-Glutamate', 'L-Glutamine', 'Glycine', 'L-Histidine', 'L-Isoleucine', ...
    'L-Leucine', 'L-Lysine', 'L-Methionine', 'L-Phenylalanine', ...
    'L-Proline', 'L-Serine', 'L-Threonine', 'L-Tryptophan', ...
    'L-Tyrosine', 'L-Valine'
};

% tRNA IDs (corresponding to amino acids)
trnaIDs = {
    'tRNA(Ala)', 'tRNA(Arg)', 'tRNA(Asn)', 'tRNA(Asp)', 'tRNA(Cys)', ...
    'tRNA (Glu)', 'tRNA(Gln)', 'tRNA(Gly)', 'tRNA(His)', 'tRNA(Ile)', ...
    'tRNA(Leu)', 'tRNA(Lys)', 'tRNA(Met)', 'tRNA(Phe)', ...
    'tRNA(Pro)', 'tRNA(Ser)', 'tRNA(Thr)', 'tRNA(Trp)', ...
    'tRNA(Tyr)', 'tRNA(Val)'
};

% Charged tRNA IDs
chargedTrnaIDs = {
    'L-Alanyl-tRNA(Ala)', 'L-Arginyl-tRNA(Arg)', 'L-Asparaginyl-tRNA(Asn)', ...
    'L-Aspartyl-tRNA(Asp)', 'L-Cysteinyl-tRNA(Cys)', 'L-Glutamyl-tRNA(Glu)', ...
    'L-Glutaminyl-tRNA(Gln)', 'Glycyl-tRNA(Gly)', 'L-Histidyl-tRNA(His)', ...
    'L-Isoleucyl-tRNA(Ile)', 'L-Leucyl-tRNA(Leu)', 'L-Lysine-tRNA (Lys)', ...
    'L-Methionyl-tRNA (Met)', 'L-Phenylalanyl-tRNA(Phe)', 'L-Prolyl-tRNA(Pro)', ...
    'L-Seryl-tRNA(Ser)', 'L-Threonyl-tRNA(Thr)', 'L-Tryptophanyl-tRNA(Trp)', ...
    'L-Tyrosyl-tRNA(Tyr)', 'L-Valyl-tRNA(Val)'
};

% Count amino acids in sequence
aaCounts = zeros(1, 20);
for i = 1:20
    aaCounts(i) = count(aaSequence, aaList(i));
end

proteinLength = length(aaSequence);

%% CORRECTED TRANSLATION STOICHIOMETRY
% Based on biochemical reality:
% - Each amino acid requires 1 charged tRNA
% - Each peptide bond requires 1 GTP for elongation + 1 GTP for translocation = 2 GTP
% - Translation initiation requires 1 GTP
% - Translation termination requires 1 GTP  
% - Each GTP hydrolysis should cost equivalent energy (use ATP for consistency)
% - Water is consumed for peptide bond formation (n-1 peptide bonds)

% Build metabolite list and stoichiometry
metabolites = {};
stoichCoeffs = [];

% 1. INPUTS (consumed, negative coefficients):

% ATP cost: 2 ATP per amino acid (for elongation) + 2 ATP for initiation/termination
atp_cost = 2 * proteinLength + 2;
metabolites{end+1} = 'ATP[c]';
stoichCoeffs(end+1) = -atp_cost;

% Water for peptide bonds (n-1 peptide bonds)
water_cost = proteinLength - 1;
metabolites{end+1} = 'H2O[c]';
stoichCoeffs(end+1) = -water_cost;

% Charged tRNAs (one for each amino acid in the sequence)
for i = 1:20
    if aaCounts(i) > 0
        metabolites{end+1} = sprintf('%s[c]', chargedTrnaIDs{i});
        stoichCoeffs(end+1) = -aaCounts(i);
    end
end

% 2. OUTPUTS (produced, positive coefficients):

% The protein product
metabolites{end+1} = sprintf('%s_aa[c]', rxnName);
stoichCoeffs(end+1) = 1;

% Free tRNAs (same number as charged tRNAs consumed)
for i = 1:20
    if aaCounts(i) > 0
        metabolites{end+1} = sprintf('%s[c]', trnaIDs{i});
        stoichCoeffs(end+1) = aaCounts(i);
    end
end

% ADP and phosphate from ATP hydrolysis
metabolites{end+1} = 'ADP[c]';
stoichCoeffs(end+1) = atp_cost;

metabolites{end+1} = 'Phosphate[c]';
stoichCoeffs(end+1) = atp_cost;

% Protons from ATP hydrolysis
metabolites{end+1} = 'H+[c]';
stoichCoeffs(end+1) = atp_cost;

%% Add the protein metabolite (if not already present)
proteinMet.mets = {sprintf('%s_aa_c', rxnName)};
proteinMet.metNames = {sprintf('%s_aa', rxnName)};
proteinMet.compartments = {'c'};

try
    model = addMets(model, proteinMet);
catch
    % Metabolite already exists
end

%% Add the translation reaction
rxnStruct.rxns = {sprintf('%s_aa', rxnName)};
rxnStruct.rxnNames = rxnStruct.rxns;
rxnStruct.mets = {metabolites};
rxnStruct.stoichCoeffs = {stoichCoeffs};
rxnStruct.lb = [0];  % Irreversible
rxnStruct.ub = [1000];  % Set reasonable upper bound
rxnStruct.catBy = {'Ribosome'};  % Catalyzed by ribosome
rxnStruct.eccodes = {''};  % EC for peptidyl transferase

model = addCatalyzedRxns(model, struct('mets', {}, 'metNames', {}, 'compartments', {}), rxnStruct, rxnName);
%fprintf('Added translation reaction for %s: %d amino acids, %d ATP cost\n', ...
%    rxnName, proteinLength, atp_cost);

end