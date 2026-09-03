function [model, startMet] = addNglyGolgiReactions(model, proteinInfo, rxnName, startMet)
% ADDNGLYGOLGIREACTIONS Adds mannosylation reactions in the Golgi.
%
% This function simulates N-glycan extension in the Golgi apparatus after 
% ER processing, following yeast-type hypermannosylation.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - structure with numOGlyc (used as proxy for glycan extension level)
%   rxnName     - base reaction name (e.g., GQ67_XXXX)
%   startMet    - current glycoprotein in Golgi [g]
%
% OUTPUTS:
%   model       - updated GEM
%   startMet    - final extended glycoprotein metabolite

% ------------------------------
% Reaction cascade:
% GNG_G1 → GNG_G4
% Each step adds mannose residues via GDP-mannose donor
% Stoichiometry is estimated based on:
% - Och1p = +1 Man
% - MPOL I = +9 Man
% - MPOL II = +30 Man
% - Mnn1/2/5 = +1 Man
%
% Source: Orlean & Menon (2007), Biochim Biophys Acta
% ------------------------------

rxnStruct.rxns = {
    sprintf('%s_ngly_golgi1', rxnName), ...
    sprintf('%s_ngly_golgi2', rxnName), ...
    sprintf('%s_ngly_golgi3', rxnName), ...
    sprintf('%s_ngly_golgi4', rxnName)
};
rxnStruct.rxnNames = rxnStruct.rxns;
rxnStruct.lb       = [0, 0, 0, 0];

% Glycoform transition intermediates
rxnStruct.mets = {
    {startMet, 'GDP-D-mannose[g]', 'GDP[g]', sprintf('%s_GNG_G1[g]', startMet(1:end-3))}, ...
    {sprintf('%s_GNG_G1[g]', startMet(1:end-3)), 'GDP-D-mannose[g]', 'GDP[g]', sprintf('%s_GNG_G2[g]', startMet(1:end-3))}, ...
    {sprintf('%s_GNG_G2[g]', startMet(1:end-3)), 'GDP-D-mannose[g]', 'GDP[g]', sprintf('%s_GNG_G3[g]', startMet(1:end-3))}, ...
    {sprintf('%s_GNG_G3[g]', startMet(1:end-3)), 'GDP-D-mannose[g]', 'GDP[g]', sprintf('%s_GNG_G4[g]', startMet(1:end-3))}
};

n = proteinInfo.numOGlyc;  % proxy for glycan extensibility

rxnStruct.stoichCoeffs = {
    [-1, -n, n, 1], ...
    [-1, -9*n, 9*n, 1], ...
    [-1, -30*n, 30*n, 1], ...
    [-1, -n, n, 1]
};

% ------------------------------
% Catalysts: Golgi mannosyltransferases
% ------------------------------
rxnStruct.catBy = {
    'sec_Och1p_complex', ...
    'sec_MPOLI_complex', ...
    'sec_MPoLII_complex', ...
    'sec_Mnn1p_Mnn2p_Mnn5p_complex'
};

% ------------------------------
% EC codes for known transferases:
% Och1p        — EC 2.4.1.232  % corrected from previous
% MPOL I/II    — EC 2.4.1.47
% Mnn1p/2p/5p  — EC 2.4.1.132
% ------------------------------
rxnStruct.eccodes = {
    '2.4.1.232', ...   % Och1p (initiation of α‑1,6-mannose)
    '2.4.1.47',   ...  % MPOL I
    '2.4.1.47',   ...  % MPOL II
    '2.4.1.132'       % Mnn1p/2p/5p
};


% ------------------------------
% Add intermediate glycoprotein metabolites
% ------------------------------
metsToAdd.mets = {
    sprintf('%s_GNG_G1_g', startMet(1:end-3)), ...
    sprintf('%s_GNG_G2_g', startMet(1:end-3)), ...
    sprintf('%s_GNG_G3_g', startMet(1:end-3)), ...
    sprintf('%s_GNG_G4_g', startMet(1:end-3))
};
metsToAdd.metNames = {
    sprintf('%s_GNG_G1', startMet(1:end-3)), ...
    sprintf('%s_GNG_G2', startMet(1:end-3)), ...
    sprintf('%s_GNG_G3', startMet(1:end-3)), ...
    sprintf('%s_GNG_G4', startMet(1:end-3))
};
metsToAdd.compartments = {'g', 'g', 'g', 'g'};

% ------------------------------
% Add reactions and metabolites to model
% ------------------------------
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Final metabolite for downstream use
startMet = sprintf('%s_GNG_G4[g]', startMet(1:end-3));
end
