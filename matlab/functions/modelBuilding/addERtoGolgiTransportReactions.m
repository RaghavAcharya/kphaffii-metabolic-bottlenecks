function [model, startMet] = addERtoGolgiTransportReactions(model, proteinInfo, rxnName, startMet)
% ADDERTOGOLGITRANSPORTREACTIONS Adds COPII-mediated transport to Golgi.
%
% This function models ER-to-Golgi transport via COPII-coated vesicles.
% Transport is energy-dependent and can be conditionally routed based
% on whether the protein is transmembrane-anchored or GPI-anchored.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - structure with fields GPIAddition and transmembrane
%   rxnName     - base protein name
%   startMet    - peptide metabolite in ER lumen [r]
%
% OUTPUTS:
%   model       - updated GEM
%   startMet    - Golgi-localized form of the peptide

% ------------------------------
% Define reaction: GTP-driven vesicle transport from ER [r] to Golgi [g]
% Energy cost based on:
% - ~2 GTP hydrolyzed per vesicle coat cycle (Sec23/Sec24 + Sar1 GTPase)
% Ref: Barlowe et al. 1994 (Cell), Springer et al. 1999 (JCB)
% ------------------------------
rxnStruct.rxns     = sprintf('%s_ER2Golgi', rxnName);
rxnStruct.rxnNames = rxnStruct.rxns;

rxnStruct.mets = {
    startMet, ...
    'H2O[c]', ...
    'GTP[c]', ...
    'H+[c]', ...
    'Phosphate[c]', ...
    'GDP[c]', ...
    sprintf('%s[g]', startMet(1:end-3))
};
rxnStruct.lb = [0];
rxnStruct.stoichCoeffs = [-1, -2, -2, 2, 2, 2, 1];

% ------------------------------
% Select appropriate catalyst complex based on membrane anchoring
% Ref: Robinson et al., Nature Rev Mol Cell Biol 2015; COPII vesicle sorting
% ------------------------------
if proteinInfo.transmembrane == 1
    rxnStruct.catBy = 'sec_COPII_transmembrane_complex';
elseif proteinInfo.GPIAddition == 1
    rxnStruct.catBy = 'sec_COPII_GPI_complex';
else
    rxnStruct.catBy = 'sec_COPII_general_complex';
end

% ------------------------------
% EC Classification: Vesicle-mediated transport, no assigned EC
% Use placeholder EC class
% ------------------------------
rxnStruct.eccodes = {'-.-.-.-'};  % no EC for vesicle formation

% ------------------------------
% Add Golgi-localized metabolite
% ------------------------------
metsToAdd.mets         = {sprintf('%s_g', startMet(1:end-3))};
metsToAdd.metNames     = {sprintf('%s', startMet(1:end-3))};
metsToAdd.compartments = {'g'};

% ------------------------------
% Add transport reaction to model
% ------------------------------
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Return updated metabolite
startMet = sprintf('%s[g]', startMet(1:end-3));
end
