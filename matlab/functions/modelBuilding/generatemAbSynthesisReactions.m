function [model] = generatemAbSynthesisReactions(model, rxnName, proteinInfoLC, proteinInfoHC, numInterchainDisulfide)
%GENERATEMABSYNTHESISREACTIONS
%   Adds synthesis, folding, assembly, and secretion reactions for a monoclonal antibody.
%
% INPUTS:
%   model         - GEM model
%   rxnName       - prefix (e.g., 'mAb1')
%   proteinInfoLC - structure for light chain (same format as generateSynthesisReactions4)
%   proteinInfoHC - structure for heavy chain
%   numInterchainDisulfide - (optional, default 0) number of interchain
%       disulfide bonds formed during assembly (e.g. IgG1 hinge + HC-LC
%       bonds - see addInterchainDSBReaction.m for the full rationale and
%       a literature-sourced default for standard IgG1 architecture).
%       Added July 2026 (shadow_exploration audit) as a SEPARATE,
%       glutathione-coupled reaction from the ATP-only intrachain folding
%       in addFoldingReaction.m/addDSBReactions.m below, because these
%       two bond classes are mechanistically and empirically distinct
%       (see addInterchainDSBReaction.m header). Defaults to 0 so
%       existing protein specs (trastuzumab, Adalimumab) are unaffected
%       until this is deliberately populated per protein - proteinInfoLC/
%       HC.numDisulfide currently represent each chain's total disulfide
%       count from ppPSIMv3_RP.csv, which has NOT yet been split into
%       intrachain vs. interchain; that split needs to happen before
%       nonzero values here would be meaningful rather than double- or
%       under-counting bonds already priced by addDSBReactions.m.
    if nargin < 5 || isempty(numInterchainDisulfide)
        numInterchainDisulfide = 0;
    end

%% Step 1: Add translation and PTMs for LC and HC

% Light Chain synthesis + processing
model = addTranslationReaction(model, proteinInfoLC.aaSequence, proteinInfoLC.dnaSequence, [rxnName '_LC']);
startMetLC = sprintf('%s_LC_aa[c]', rxnName);
if proteinInfoLC.throughER > 0
    % FIX (shadow_exploration audit, July 2026): translocate cytosolic
    % [c] peptide to the ER lumen [r] before addPTMs - see
    % addERTranslocationReactions.m header for the gap this closes.
    [model, startMetLC] = addERTranslocationReactions(model, proteinInfoLC, [rxnName '_LC'], startMetLC);
    [model, startMetLC] = addPTMs(model, proteinInfoLC, [rxnName '_LC'], startMetLC);
end

% Heavy Chain synthesis + processing
model = addTranslationReaction(model, proteinInfoHC.aaSequence, proteinInfoHC.dnaSequence, [rxnName '_HC']);
startMetHC = sprintf('%s_HC_aa[c]', rxnName);
if proteinInfoHC.throughER > 0
    [model, startMetHC] = addERTranslocationReactions(model, proteinInfoHC, [rxnName '_HC'], startMetHC);
    [model, startMetHC] = addPTMs(model, proteinInfoHC, [rxnName '_HC'], startMetHC);
end

%% Step 2: Add folding for each chain (simplified to 1-step)
[model, startMetLC] = addFoldingReaction(model, startMetLC, [rxnName '_LC']);
[model, startMetHC] = addFoldingReaction(model, startMetHC, [rxnName '_HC']);

%% Step 3: Assembly (LC + HC => intermediate complex)
intermediateName = sprintf('%s_LC_HC', rxnName);
intermediateComp = 'r';
intermediateFull = sprintf('%s[%s]', intermediateName, intermediateComp);

rxnStructAssembly.rxns         = {sprintf('%s_assembly', rxnName)};
rxnStructAssembly.rxnNames     = rxnStructAssembly.rxns;
rxnStructAssembly.mets         = {startMetLC, startMetHC, intermediateFull};
rxnStructAssembly.stoichCoeffs = [-1, -1, 1];
rxnStructAssembly.catBy        = {''};

metsToAdd.mets         = {strrep(intermediateFull, '[', '_')};
metsToAdd.metNames     = {intermediateName};
metsToAdd.compartments = {intermediateComp};
model = addCatalyzedRxns(model, metsToAdd, rxnStructAssembly, rxnName);

%% Step 3b: Interchain disulfide bond formation (added July 2026)
% Previously the LC+HC intermediate went straight into dimerization with
% no representation of interchain disulfide bond formation at all - see
% addInterchainDSBReaction.m for why that's the specific gap that kept
% the model from ever being able to predict a mAb-specific glutathione
% effect. No-ops (returns intermediateFull unchanged) when
% numInterchainDisulfide is 0, i.e. for any spec that hasn't opted in yet.
[model, intermediateFull] = addInterchainDSBReaction(model, numInterchainDisulfide, rxnName, intermediateFull);

%% Step 4: Dimerization (2x intermediate => mAb)
mAbName = sprintf('%s_mAb', rxnName);
mAbComp = 'r';
mAbFull = sprintf('%s[%s]', mAbName, mAbComp);

rxnStructDimer.rxns         = {sprintf('%s_dimerization', rxnName)};
rxnStructDimer.rxnNames     = rxnStructDimer.rxns;
rxnStructDimer.mets         = {intermediateFull, mAbFull};
rxnStructDimer.stoichCoeffs = [-2, 1];

rxnStructDimer.catBy        = {''};

metsToAdd.mets         = {strrep(mAbFull, '[', '_')};
metsToAdd.metNames     = {mAbName};
metsToAdd.compartments = {mAbComp};
model = addCatalyzedRxns(model, metsToAdd, rxnStructDimer, rxnName);

%% Step 5: ER to Golgi transport (if needed)
[model, startMet] = addERtoGolgiTransportReactions(model, proteinInfoHC, rxnName, mAbFull);

%% Step 6: Golgi N-glycosylation
[model, startMet] = addNglyGolgiReactions(model, proteinInfoHC, rxnName, startMet);

%% Step 7: Maturation
matureName = sprintf('%s_mature', rxnName);
matureFull = sprintf('%s[g]', matureName);
maturationMet = struct( ...
    'mets',         {{sprintf('%s_mature_g', rxnName)}}, ...
    'metNames',     {{matureName}}, ...
    'compartments', {{'g'}});

maturationRxn = struct( ...
    'rxns',         {{sprintf('%s_maturation', rxnName)}}, ...
    'rxnNames',     {{sprintf('%s_maturation', rxnName)}}, ...
    'mets',         {{startMet, matureFull}}, ...
    'stoichCoeffs', {[-1, 1]}, ...
    'catBy',        {{''}});

model = addMets(model, maturationMet);
model = addCatalyzedRxns(model, struct('mets', {}, 'metNames', {}, 'compartments', {}), maturationRxn, rxnName);

%% Step 8: Final transport to final compartment (if needed)
model = addTransportReactions(model, proteinInfoHC, rxnName, maturationMet);

end
