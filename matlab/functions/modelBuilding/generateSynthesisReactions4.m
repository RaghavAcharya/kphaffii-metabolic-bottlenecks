function [model] = generateSynthesisReactions4(model, rxnName, proteinInfo)
% GENERATESYNTHESISREACTIONS4
% Adds all relevant synthesis, processing, and degradation reactions for a gene.
%
% INPUTS:
%   model       - GEM structure
%   rxnName     - canonical gene ID (e.g., GQ67_XXXX)
%   proteinInfo - struct with required fields:
%       aaSequence     - 1-letter amino acid sequence
%       dnaSequence    - DNA coding sequence (ATGC)
%       throughER      - 0 = cytosolic; 1 = ER co-translational; 2 = ER post-translational
%       numNGlyc       - number of N-linked glycosylation sites
%       numDisulfide   - number of disulfide bonds
%       numOGlyc       - number of O-linked glycosylation sites
%       finalLocation  - compartment letter (e.g., 'c', 'g', 'e', etc.)

%% Step 1: Add core transcription, translation, and synthesis reactions
model = addTranslationReaction(model, ...
    proteinInfo.aaSequence, ...
    proteinInfo.dnaSequence, ...
    rxnName);

% Infer the peptide product generated at end of translation
startMet = sprintf('%s_aa[c]', rxnName);

%% Step 2: Add post-translational modifications (ER stage), ER->Golgi
%% transport, and Golgi-stage glycan extension
%
% FIX (shadow_exploration audit, July 2026): this step previously only
% called addPTMs (ER-lumen PTMs: disulfide, N-gly-ER, O-gly-ER) despite
% the old comment here claiming "Add ER and Golgi PTMs." The actual
% ER->Golgi COPII transport reaction (addERtoGolgiTransportReactions.m)
% and Golgi-stage glycan extension (addNglyGolgiReactions.m,
% addOglyGolgiReactions.m) were never called anywhere in this
% single-chain pipeline - only generatemAbSynthesisReactions.m called the
% first two. Real single-chain RPs used throughout this project (HSA,
% Etanercept, EPO, HRP, Factor Xa - see data/ppPSIMv3_RP.csv) have
% nonzero N-gly/O-gly counts, so this was a real, consequential gap: the
% model silently jumped straight from an ER-lumen peptide to a
% Golgi-tagged "mature[g]" product via one uncatalyzed, cost-free
% reaction in Step 3 below, with no COPII GTP cost, no COPII
% enzyme-capacity requirement, and no Golgi-stage glycan extension at
% all.
%
% This now mirrors generatemAbSynthesisReactions.m's Steps 5-6: ER->Golgi
% transport, then Golgi N-gly extension. addNglyGolgiReactions.m is
% called unconditionally after transport (not gated on numNGlyc), exactly
% matching how the mAb pipeline calls it - including that function's
% pre-existing "use numOGlyc as a glycan-extensibility proxy" quirk
% (see its own header), which is preserved as-is here rather than
% redesigned, since that specific quirk wasn't part of this fix request.
% At numOGlyc==0 this reduces to a zero-substrate passthrough/renaming
% step, same as it already does in the mAb pipeline.
%
% addOglyGolgiReactions.m has no precedent call site anywhere in the
% codebase (confirmed via grep) - it is genuinely new here, gated on
% numOGlyc > 0, since several real single-chain RPs in this project do
% have nonzero O-glyc counts and the function was already fully built,
% just unused.
%
% NOTE: addERtoGolgiTransportReactions.m reads proteinInfo.transmembrane
% and proteinInfo.GPIAddition unconditionally (no isfield guard inside
% that function). Every proteinInfo-building call site in this codebase
% (buildProteinInfo.m, extractProteinInfo.m, runComprehensiveValidation.m,
% and every local re-implementation, verified via grep) populates both
% fields, so this is safe in practice; defensive defaults are set below
% anyway in case a stripped-down proteinInfo struct is ever passed in.
if proteinInfo.throughER > 0
    % FIX (shadow_exploration audit, July 2026): translocate the peptide
    % from the cytosol [c] (where addTranslationReaction.m always emits
    % it) into the ER lumen [r] BEFORE calling addPTMs - see
    % addERTranslocationReactions.m's header for the full gap this closes
    % (every ER-lumen PTM function assumed [r] input; nothing ever
    % produced it).
    [model, startMet] = addERTranslocationReactions(model, proteinInfo, rxnName, startMet);

    % Add ER-lumen PTMs (disulfide, N-gly-ER, O-gly-ER, GPI anchor)
    [model, startMet] = addPTMs(model, proteinInfo, rxnName, startMet);

    % Defensive defaults for fields addERtoGolgiTransportReactions.m
    % reads unconditionally (see NOTE above).
    if ~isfield(proteinInfo, 'transmembrane') || isempty(proteinInfo.transmembrane)
        proteinInfo.transmembrane = 0;
    end
    if ~isfield(proteinInfo, 'GPIAddition') || isempty(proteinInfo.GPIAddition)
        proteinInfo.GPIAddition = 0;
    end

    % ER -> Golgi COPII-vesicle transport (GTP-costed, enzyme-capacitated)
    [model, startMet] = addERtoGolgiTransportReactions(model, proteinInfo, rxnName, startMet);

    % Golgi-stage N-glycan extension (mirrors mAb pipeline; unconditional
    % after transport, per that pipeline's own convention)
    [model, startMet] = addNglyGolgiReactions(model, proteinInfo, rxnName, startMet);

    % Golgi-stage O-glycan extension (new: no prior call site existed)
    if isfield(proteinInfo, 'numOGlyc') && proteinInfo.numOGlyc > 0
        [model, startMet] = addOglyGolgiReactions(model, proteinInfo, rxnName, startMet);
    end
end

%% Step 3: Maturation Step
% Define the mature protein (final product after processing)
if proteinInfo.throughER > 0
    finalCompartment = 'g';  % protein exits Golgi
else
    finalCompartment = 'c';  % cytoplasmic mature form
end

maturationMet.mets         = {sprintf('%s_mature_%s', rxnName, finalCompartment)};
maturationMet.metNames     = {sprintf('%s_mature', rxnName)};
maturationMet.compartments = {finalCompartment};

maturationRxn.rxns         = {sprintf('%s_maturation', rxnName)};
maturationRxn.rxnNames     = maturationRxn.rxns;
maturationRxn.lb           = [0];
maturationRxn.mets         = {startMet, sprintf('%s_mature[%s]', rxnName, finalCompartment)};
maturationRxn.stoichCoeffs = [-1, 1];
maturationRxn.catBy        = {''};           % Required for addCatalyzedRxns
maturationRxn.eccodes      = {''};           % Optional metadata field

% Add the mature protein metabolite (if not already present)
try
    model = addMets(model, maturationMet);
catch
    warning('%s is already present as a metabolite.', rxnName)
end

% Add the maturation reaction
model = addCatalyzedRxns(model, struct('mets', {}, 'metNames', {}, 'compartments', {}), maturationRxn, rxnName);


% Update the starting point for downstream functions
startMet = sprintf('%s_mature[%s]', rxnName, finalCompartment);

%% Step 4: Add protein transport reactions (handled separately)
model = addTransportReactions(model, proteinInfo, rxnName, maturationMet);

%% Step 5: Add misfolding/refolding/degradation reactions 
model = addMisfoldingReactions(model, proteinInfo, rxnName, startMet);

end
