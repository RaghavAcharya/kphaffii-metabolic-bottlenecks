function [model] = addCatalyzedRxns(model, metStruct, rxnStruct, rxnName)
% ADDCATALYZEDRXNS Adds catalyzed reactions and tracks their catalysts.
% This function adds new metabolites and reactions to the model,
% while also tracking which metabolites catalyze which reactions using
% a custom field 'catBy'.

% INPUTS:
%   model     - the GEM model (RAVEN format assumed)
%   metStruct - structure of metabolites to add (compatible with addMets)
%   rxnStruct - structure of reactions to add (compatible with addRxns)
%   rxnName   - name of the reaction (for warnings)

% Add metabolites (if not already present)
try
    model = addMets(model, metStruct);
catch ME
    warning('%s: some or all metabolites may already exist. Skipping addMets.', rxnName);
end

% Make sure catBy is a cell array
if ischar(rxnStruct.catBy)
    rxnStruct.catBy = {rxnStruct.catBy};
elseif isstring(rxnStruct.catBy)
    rxnStruct.catBy = cellstr(rxnStruct.catBy);
end

% Initialize or extend 'catBy'
if ~isfield(model, 'catBy')
    model.catBy = cell(length(model.rxns), 1);
elseif length(model.catBy) < length(model.rxns)
    model.catBy(end+1:length(model.rxns)) = {[]};
end

% Extend 'catBy' to match new reactions
model.catBy(end+1:end+length(rxnStruct.catBy)) = rxnStruct.catBy;

% Add reactions to the model
% RAVEN-compatible call; will need adjustment for COBRA
% (e.g., using addReaction or batch conversion)
model = addRxns(model, rxnStruct, 3, [], true, true);

end
