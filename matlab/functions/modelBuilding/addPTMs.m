function [model, startMet] = addPTMs(model, proteinInfo, rxnName, startMet)
% ADDPTMS
% Adds post-translational modifications to a protein (ER-lumen stage).
%
% INPUTS:
%   model       - the GEM structure
%   proteinInfo - struct with fields like numNGlyc, numDisulfide, GPIAddition, etc.
%   rxnName     - gene or protein ID string (e.g., GQ67_XXXX_HC)
%   startMet    - ID of the translated peptide metabolite
%
% FIX (shadow_exploration audit, July 2026): added GPI-anchor attachment
% (addGPIRxns.m). That function already existed, fully built, but was
% never called from anywhere in the codebase - neither this shared
% ER-stage orchestrator (used by both the single-chain and mAb synthesis
% pipelines) nor either pipeline directly. GPI-transamidase attachment is
% biologically an ER-lumen event (occurs before ER exit), so it belongs
% here, gated on proteinInfo.GPIAddition (the field name this codebase
% actually populates from ppPSIMv3_RP.csv's GPIsite column - see
% buildProteinInfo.m/extractProteinInfo.m/runComprehensiveValidation.m).
% Called last among the ER PTMs (order among these uncoupled, sequential
% mass-balance conversions doesn't matter for a coarse-grained
% stoichiometric model). Per data/ppPSIMv3_RP.csv, GPIsite=1 mostly hits
% native/host yeast genes (GPI-anchored cell wall proteins), not the
% recombinant single-chain proteins investigated elsewhere in this
% project - but generateSynthesisReactions4.m is called for every
% catalytically relevant gene during model expansion, not just RPs, so
% this gap was real and consequential for those native genes too.
%
% NOTE: like addDSBReactions.m/addNglyERReactions.m/addOglyERReactions.m,
% addGPIRxns.m assumes its input metabolite already carries an [r]
% (ER lumen) compartment suffix. That assumption is now actually correct
% at every call site: both generateSynthesisReactions4.m and
% generatemAbSynthesisReactions.m call addERTranslocationReactions.m
% immediately before addPTMs whenever throughER > 0, which performs the
% cytosol->ER translocation (previously missing entirely - translation
% always emitted "..._aa[c]" and nothing ever converted it to [r]) and
% always returns an [r]-tagged metabolite regardless of route.

% Initialize current metabolite
currentMet = startMet;

% Add disulfide bonds
if isfield(proteinInfo, 'numDisulfide') && proteinInfo.numDisulfide > 0
    [model, currentMet] = addDSBReactions(model, proteinInfo, rxnName, currentMet);
end

% Add ER N-glycosylation
if isfield(proteinInfo, 'numNGlyc') && proteinInfo.numNGlyc > 0
    [model, currentMet] = addNglyERReactions(model, proteinInfo, rxnName, currentMet);
end

% Add ER O-glycosylation
if isfield(proteinInfo, 'numOGlyc') && proteinInfo.numOGlyc > 0
    [model, currentMet] = addOglyERReactions(model, proteinInfo, rxnName, currentMet);
end

% Add GPI anchor (FIX: was never called anywhere in the codebase)
if isfield(proteinInfo, 'GPIAddition') && proteinInfo.GPIAddition == 1
    [model, currentMet] = addGPIRxns(model, proteinInfo, rxnName, currentMet);
end

startMet = currentMet;
end
