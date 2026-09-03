function [model, startMet] = addNglyERReactions(model, proteinInfo, rxnName, startMet)
% ADDNGLYERREACTIONS Adds ER-localized N-linked glycosylation reactions.
%
% This function simulates the sequential addition and trimming of a glycan 
% in the ER lumen. Each reaction is catalyzed by a known complex.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - structure with numNGlyc (number of N-glyc sites)
%   rxnName     - base name of the protein/reaction
%   startMet    - starting peptide metabolite with [r] suffix
%
% OUTPUTS:
%   model       - updated GEM
%   startMet    - updated glycosylated metabolite name

% ------------------------------
% Reaction IDs (5-step canonical N-glycan maturation in ER)
% References: 
%   KEGG map00510: N-Glycan biosynthesis
%   
% ------------------------------
rxnStruct.rxns = {
    sprintf('%s_ngly_er1_transfer', rxnName), ...
    sprintf('%s_ngly_er2_trim1', rxnName), ...
    sprintf('%s_ngly_er3_trim2', rxnName), ...
    sprintf('%s_ngly_er4_trim3', rxnName), ...
    sprintf('%s_ngly_er5_demanosylation', rxnName)
};
rxnStruct.rxnNames = rxnStruct.rxns;

% ------------------------------
% Metabolite sequences (each intermediate is a trimmed form)
% Naming convention: G3M9 → G2M9 → G1M9 → M9 → M8
% ------------------------------
rxnStruct.mets = {
    {startMet, '(Glc)3(GlcNAc)2(Man)9(PP-Dol)1[r]', 'Dolichyl phosphate[r]', sprintf('%s_G3M9[r]', startMet(1:end-3))}, ...
    {sprintf('%s_G3M9[r]', startMet(1:end-3)), 'H2O[r]', 'D-Glucose[r]', sprintf('%s_G2M9[r]', startMet(1:end-3))}, ...
    {sprintf('%s_G2M9[r]', startMet(1:end-3)), 'H2O[r]', 'D-Glucose[r]', sprintf('%s_G1M9[r]', startMet(1:end-3))}, ...
    {sprintf('%s_G1M9[r]', startMet(1:end-3)), 'H2O[r]', 'D-Glucose[r]', sprintf('%s_M9[r]', startMet(1:end-3))}, ...
    {sprintf('%s_M9[r]', startMet(1:end-3)), 'H2O[r]', 'D-Mannose[r]', sprintf('%s_M8[r]', startMet(1:end-3))}
};
rxnStruct.lb  =  [0, 0, 0, 0, 0];

n = proteinInfo.numNGlyc;

% ------------------------------
% Stoichiometries: 1 mol glycoprotein per step × n glycosylation sites
% Dol-PP donor cost only occurs in first step
% Based on KEGG map00510 and UniProt enzyme activity
% ------------------------------
rxnStruct.stoichCoeffs = {
    [-1, -n, n, 1], ...
    [-1, n, n, 1], ...
    [-1, n, n, 1], ...
    [-1, n, n, 1], ...
    [-1, n, n, 1]
};

% ------------------------------
% Catalysts (ER enzymes)
% ------------------------------
rxnStruct.catBy = {
    'sec_OSTC_complex',        % oligosaccharyltransferase complex
    'sec_Cwh41p_complex',      % alpha-glucosidase I
    'sec_Rot2p_complex',       % alpha-glucosidase II
    'sec_Rot2p_complex',       % same as above
    'sec_Mns1p_complex'        % alpha-1,2-mannosidase
};

% ------------------------------
% EC numbers (from KEGG/BRENDA/UniProt)
% ------------------------------
rxnStruct.eccodes = {
    '2.4.99.18',   % OST complex
    '',   % glucosidase I
    '',   % glucosidase II
    '',   % glucosidase II
    '3.2.1.113'    % mannosidase
};

% ------------------------------
% Add intermediate glycoproteins
% ------------------------------
metsToAdd.mets = {
    sprintf('%s_G3M9_r', startMet(1:end-3)), ...
    sprintf('%s_G2M9_r', startMet(1:end-3)), ...
    sprintf('%s_G1M9_r', startMet(1:end-3)), ...
    sprintf('%s_M9_r', startMet(1:end-3)), ...
    sprintf('%s_M8_r', startMet(1:end-3))
};
metsToAdd.metNames = {
    sprintf('%s_G3M9', startMet(1:end-3)), ...
    sprintf('%s_G2M9', startMet(1:end-3)), ...
    sprintf('%s_G1M9', startMet(1:end-3)), ...
    sprintf('%s_M9', startMet(1:end-3)), ...
    sprintf('%s_M8', startMet(1:end-3))
};
metsToAdd.compartments = {'r', 'r', 'r', 'r', 'r'};

% ------------------------------
% Add all reactions and metabolites to the model
% ------------------------------
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Return final product (M8 form)
startMet = sprintf('%s_M8[r]', startMet(1:end-3));
end
