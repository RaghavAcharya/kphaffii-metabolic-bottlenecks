function [model, startMet] = addOglyGolgiReactions(model, proteinInfo, rxnName, startMet)
% ADDOGLYGOLGIREACTIONS Adds Golgi-localized O-linked mannosylation reactions.
%
% Simulates mannose addition in Golgi via α-1,2-mannosyltransferases
% (Kre2p/Mnt1p/Ktr proteins).
%
% INPUTS:
%   model         - GEM model structure
%   proteinInfo   - includes numOGlyc (number of O-glyc sites)
%   rxnName       - base reaction name (GQ67_XXXX)
%   startMet      - current glycoform in Golgi [g]
%
% OUTPUTS:
%   model         - updated model
%   startMet      - extended glycan form in Golgi [g]

% Define two-step extension
rxnStruct.rxns     = {
    sprintf('%s_ogly_golgi_ext1', rxnName), ...
    sprintf('%s_ogly_golgi_ext2', rxnName)
};
rxnStruct.rxnNames = rxnStruct.rxns;
rxnStruct.lb = [0, 0];

n = proteinInfo.numOGlyc;
rxnStruct.mets = {
    {startMet, 'GDP-D-mannose[g]', 'GDP[g]', sprintf('%s_GOG_G1[g]', startMet(1:end-3))}, ...
    {sprintf('%s_GOG_G1[g]', startMet(1:end-3)), 'GDP-D-mannose[g]', 'GDP[g]', sprintf('%s_GOG_G2[g]', startMet(1:end-3))}
};
rxnStruct.stoichCoeffs = {
    [-1, -3*n, 3*n, 1], ...
    [-1, -2*n, 2*n, 1]
};

% Catalyst: Kre2p/Mnt1p/Ktr1p/Ktr3p complex
rxnStruct.catBy = {
    'sec_Kre2p_ktr1p_ktr3p_complex', ...
    'sec_Kre2p_ktr1p_ktr3p_complex'
};

% Correct EC code: α-1,2-mannosyltransferase (unclassified)
rxnStruct.eccodes = {
    '2.4.1.-', ...
    '2.4.1.-'
};

% Add intermediate glycoform metabolites
metsToAdd.mets         = {
    sprintf('%s_GOG_G1_g', startMet(1:end-3)), ...
    sprintf('%s_GOG_G2_g', startMet(1:end-3))
};
metsToAdd.metNames     = {
    sprintf('%s_GOG_G1', startMet(1:end-3)), ...
    sprintf('%s_GOG_G2', startMet(1:end-3))
};
metsToAdd.compartments = {'g','g'};

% Add to model
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Update output metabolite
startMet = sprintf('%s_GOG_G2[g]', startMet(1:end-3));
end
