function [model, foldedMet] = addFoldingReaction(model, inputMet, geneName)
% Add folding reaction using default chaperone
foldedBase = sprintf('%s_fold_1', inputMet(1:end-3));   % Strip [x] and append _fold_1
foldedMet  = sprintf('%s[r]', foldedBase);              % Fully qualified metabolite ID

% Define new metabolite to add
metsToAdd.mets         = {strrep(foldedMet, '[', '_')};  % ID without compartment brackets
metsToAdd.metNames     = {foldedBase};                   % Descriptive name
metsToAdd.compartments = {'r'};

% Add folding reaction
rxnStruct.rxns         = {sprintf('%s_folding_1', geneName)};
rxnStruct.rxnNames     = rxnStruct.rxns;
rxnStruct.mets         = {inputMet, 'H2O[r]', 'ATP[r]', 'H+[r]', 'Phosphate[r]', 'ADP[r]', foldedMet};
rxnStruct.stoichCoeffs = [-1, -2, -2, 2, 2, 2, 1];
rxnStruct.catBy        = {'GQ67_05237_mature[r]'};

model = addMets(model, metsToAdd);  % Explicitly add metabolite to prevent "m_0001"
model = addCatalyzedRxns(model, struct('mets', {}, 'metNames', {}, 'compartments', {}), rxnStruct, geneName);
end
