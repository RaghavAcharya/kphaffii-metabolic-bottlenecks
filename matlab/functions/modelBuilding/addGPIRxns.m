function [model, startMet] = addGPIRxns(model, proteinInfo, rxnName, startMet)
% ADDGPIRXNS Adds GPI anchor via GPI-transamidase in ER lumen.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - includes GPIsite flag (1 if GPI anchor present)
%   rxnName     - base reaction name (e.g., GQ67_XXXX)
%   startMet    - peptide metabolite ID in [r]
%
% OUTPUTS:
%   model       - updated GEM with GPI reaction added
%   startMet    - updated metabolite ID (GPI-anchored peptide)

% ------------------------------
% Reaction: covalent attachment of GPI anchor (single event)
% Catalyzed by: GPI transamidase complex (PIGK/Gaa1p/PIG-S/PIG-T)
% EC class: EC 2.6.1.- (protein-GPI transamidase) :contentReference[oaicite:6]{index=6}
% ------------------------------

rxnStruct.rxns     = {sprintf('%s_GPI1', rxnName)};
rxnStruct.rxnNames = rxnStruct.rxns;

rxnStruct.mets     = {
    startMet, ...
    'D-glucosaminylphosphatidylinositol[r]', ...
    'H2O[r]', ...
    sprintf('%s_GPI[r]', startMet(1:end-3))
};
rxnStruct.lb = [0];
rxnStruct.stoichCoeffs = {[-1, -1, 1, 1]};  % one anchor per protein

% ------------------------------
% Catalyst complex for GPI attachment
% PIG-K (cysteine protease) and associated subunits :contentReference[oaicite:7]{index=7}
% ------------------------------
rxnStruct.catBy    = {'sec_GPIR_complex'};

% ------------------------------
% Add EC code
% ------------------------------
rxnStruct.eccodes  = {'2.6.1.-'};

% ------------------------------
% Add metabolite for GPI-anchored peptide
% ------------------------------
metsToAdd.mets         = {sprintf('%s_GPI_r', startMet(1:end-3))};
metsToAdd.metNames     = {sprintf('%s_GPI', startMet(1:end-3))};
metsToAdd.compartments = {'r'};

% ------------------------------
% Add reaction and metabolite to model
% ------------------------------
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Update for downstream use
startMet = sprintf('%s_GPI[r]', startMet(1:end-3));
end
