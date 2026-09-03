function [model, startMet] = addERTranslocationReactions(model, proteinInfo, rxnName, startMet)
% ADDERTRANSLOCATIONREACTIONS Moves a cytosolic peptide into the ER lumen.
%
% FIX (shadow_exploration audit, July 2026): this function did not exist
% before. Every ER-lumen PTM function in this codebase
% (addDSBReactions.m, addNglyERReactions.m, addOglyERReactions.m, and the
% newly-wired-in addGPIRxns.m - see addPTMs.m's FIX note) explicitly
% documents its input as "peptide metabolite with [r] suffix" / "in
% [r]" - i.e. already inside the ER lumen. But addTranslationReaction.m
% ALWAYS emits "<rxnName>_aa[c]" (cytosolic), regardless of
% proteinInfo.throughER, and nothing anywhere in the codebase ever
% converted a [c]-tagged peptide to an [r]-tagged one. So whenever a
% throughER>0 protein had ALL PTM counts equal to zero (no disulfides, no
% N-gly, no O-gly, no GPI), addPTMs.m was a complete no-op and the
% "ER-processed" peptide silently stayed tagged [c] the entire time -
% and even when at least one PTM did fire, the model never actually
% priced or enzyme-capacitated the act of getting the peptide across the
% ER membrane in the first place.
%
% This was especially wasteful because data/complexData.xlsx already had
% full, literature-cited complex definitions for BOTH translocation
% routes sitting completely unused (verified via grep - no reaction
% anywhere in the codebase referenced any of these catalyst names before
% this fix):
%   Co-translational (SRP-dependent) route:
%     sec_SRPC                    - Signal Recognition Particle complex
%     coTranslationTransloc_complex - wraps sec_SRPC (the catalyst used below)
%   Post-translational (BiP/Hsp70-dependent) route:
%     sec_RAC_complex              - Ribosome-Associated Complex (holds
%                                     nascent chain unfolded after release)
%     sec_Ssa1_Ydj1_Snl1_complex   - cytosolic Hsp70 co-chaperone escort
%     postTranslationTransloc_Complex - wraps sec_Ssa1_Ydj1_Snl1_complex
%                                     (the catalyst used below)
%     sec_SEC61SEC63C_complex      - the Sec61/Sec63 translocon channel
%     sec_BIP_NEFS_complex         - lumenal BiP + nucleotide exchange
%                                     factors that ratchet the chain in
% All of these already have real gene sets and TN values in
% complexData.xlsx (with citations), so using them as catBy targets here
% lets expandGEMv6.m's existing resolveComplexGeneSet/catalyst-building
% machinery pick them up with no changes needed there or to
% complexData.xlsx.
%
% proteinInfo.throughER selects the route (per generateSynthesisReactions4.m's
% own docstring: 1 = co-translational, 2 = post-translational). In the
% project's actual data (ppPSIMv3.csv / ppPSIMv3_RP.csv) ThroughER only
% ever takes values 0 or 1 today - the post-translational branch below is
% implemented for completeness/correctness (and because complexData.xlsx
% clearly anticipates it) but is not currently exercised by any real gene
% in this dataset. Any throughER value other than 1 or 2 (but > 0) falls
% back to the co-translational route, since that's what 100% of the
% real nonzero-throughER data uses.
%
% ENERGETIC COSTS: fixed, per-event GTP/ATP costs are used for the
% chaperone/GTPase steps (matching the existing convention in this
% codebase - e.g. addDSBReactions.m's fixed 2.775 ATP/bond - rather than
% inventing a new per-residue scaling scheme with no precedent). The
% single exception is the BiP/NEF lumenal ratchet step, which the
% literature ties fairly directly to chain length (each Hsp70-ATPase
% binding/release cycle threads roughly one ~30-40 residue segment; see
% Matlack, Mothes & Rapoport, Cell 1998/1999 review of BiP-driven
% post-translational import), so that one step scales with
% ceil(length/40) ATP. The SRP-GTPase co-translational targeting cost (2
% GTP: one each for SRP54 and the SRP receptor) is a commonly cited
% fixed, per-targeting-event figure (Rapoport, Nature 2007, "Protein
% translocation across the eukaryotic ER and bacterial plasma
% membrane"), not a per-residue one, since SRP engages once per
% targeting event rather than once per residue. All of these are
% order-of-magnitude, literature-informed estimates, not exact figures -
% flagged here for future refinement, same spirit as this project's other
% "aggregate literature estimate" costs.
%
% The Sec61/Sec63 channel-crossing step itself is modeled as
% energy-free (a passive protein-conducting channel - the actual driving
% force is the lumenal BiP ratchet, priced in the following step), so
% that reaction only carries the sec_SEC61SEC63C_complex enzyme-capacity
% requirement and a compartment change, no cofactor cost.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - struct with throughER (1 or 2) and aaSequence
%   rxnName     - base protein/gene name (e.g., GQ67_XXXX or GQ67_XXXX_HC)
%   startMet    - translated peptide metabolite, "<rxnName>_aa[c]"
%
% OUTPUTS:
%   model       - updated GEM
%   startMet    - ER-lumen peptide metabolite, "<rxnName>_er[r]" -
%                 downstream addPTMs.m/addERtoGolgiTransportReactions.m
%                 etc. all expect exactly this: an [r]-tagged metabolite.

if ~isfield(proteinInfo, 'throughER') || proteinInfo.throughER <= 0
    % Nothing to do - caller should not have called this for a cytosolic
    % protein, but no-op defensively rather than error.
    return;
end

route = proteinInfo.throughER;
if route ~= 1 && route ~= 2
    route = 1; % fallback: co-translational, matches 100% of real data today
end

if route == 1
    %% Co-translational (SRP-dependent) route: one lumped step
    rxnStruct.rxns     = {sprintf('%s_cotranslocation', rxnName)};
    rxnStruct.rxnNames = rxnStruct.rxns;
    rxnStruct.lb       = [0];

    rxnStruct.mets = {
        startMet, ...
        'GTP[c]', ...
        'H2O[c]', ...
        'GDP[c]', ...
        'Phosphate[c]', ...
        'H+[c]', ...
        sprintf('%s_er[r]', rxnName)
    };
    % 2 GTP: SRP54 + SRP-receptor GTPase cycle per targeting event
    % (Rapoport, Nature 2007)
    rxnStruct.stoichCoeffs = [-1, -2, -2, 2, 2, 2, 1];

    rxnStruct.catBy   = {'coTranslationTransloc_complex'};
    rxnStruct.eccodes = {'-.-.-.-'}; % translocation machinery, no assigned EC

    metsToAdd.mets         = {sprintf('%s_er_r', rxnName)};
    metsToAdd.metNames     = {sprintf('%s_er', rxnName)};
    metsToAdd.compartments = {'r'};

    model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);
    startMet = sprintf('%s_er[r]', rxnName);

else
    %% Post-translational (BiP/Hsp70-dependent) route: four sequential steps
    % Not exercised by any gene in the current dataset (ThroughER never
    % equals 2 today) but implemented in full since complexData.xlsx
    % already anticipates it.

    proteinLength = 0;
    if isfield(proteinInfo, 'aaSequence') && ~isempty(proteinInfo.aaSequence)
        proteinLength = length(proteinInfo.aaSequence);
    end

    % --- Step A: RAC holds the released nascent chain unfolded (no
    % cofactor cost - RAC is a passive holdase) ---
    rxnStructA.rxns         = {sprintf('%s_rac_hold', rxnName)};
    rxnStructA.rxnNames     = rxnStructA.rxns;
    rxnStructA.lb           = [0];
    rxnStructA.mets         = {startMet, sprintf('%s_rac[c]', rxnName)};
    rxnStructA.stoichCoeffs = [-1, 1];
    rxnStructA.catBy        = {'sec_RAC_complex'};
    rxnStructA.eccodes      = {''};

    metsToAddA.mets         = {sprintf('%s_rac_c', rxnName)};
    metsToAddA.metNames     = {sprintf('%s_rac', rxnName)};
    metsToAddA.compartments = {'c'};

    model = addCatalyzedRxns(model, metsToAddA, rxnStructA, rxnName);
    metA = sprintf('%s_rac[c]', rxnName);

    % --- Step B: Ssa1p/Ydj1p/Snl1p Hsp70 co-chaperone escort to the ER
    % membrane (fixed ATP cost per event, same convention as
    % addDSBReactions.m's fixed per-bond ATP cost) ---
    n_escort = 3; % fixed, approximate Hsp70 ATPase-cycle cost per escort event
    rxnStructB.rxns     = {sprintf('%s_chaperone_escort', rxnName)};
    rxnStructB.rxnNames = rxnStructB.rxns;
    rxnStructB.lb       = [0];
    rxnStructB.mets = {
        metA, ...
        'ATP[c]', ...
        'H2O[c]', ...
        'ADP[c]', ...
        'Phosphate[c]', ...
        'H+[c]', ...
        sprintf('%s_escorted[c]', rxnName)
    };
    rxnStructB.stoichCoeffs = [-1, -n_escort, -n_escort, n_escort, n_escort, n_escort, 1];
    rxnStructB.catBy        = {'postTranslationTransloc_Complex'};
    rxnStructB.eccodes      = {''};

    metsToAddB.mets         = {sprintf('%s_escorted_c', rxnName)};
    metsToAddB.metNames     = {sprintf('%s_escorted', rxnName)};
    metsToAddB.compartments = {'c'};

    model = addCatalyzedRxns(model, metsToAddB, rxnStructB, rxnName);
    metB = sprintf('%s_escorted[c]', rxnName);

    % --- Step C: Sec61/Sec63 channel crossing (energy-free; passive
    % channel, driving force is the lumenal BiP ratchet in Step D) ---
    rxnStructC.rxns         = {sprintf('%s_translocon_crossing', rxnName)};
    rxnStructC.rxnNames     = rxnStructC.rxns;
    rxnStructC.lb           = [0];
    rxnStructC.mets         = {metB, sprintf('%s_translocon[r]', rxnName)};
    rxnStructC.stoichCoeffs = [-1, 1];
    rxnStructC.catBy        = {'sec_SEC61SEC63C_complex'};
    rxnStructC.eccodes      = {''};

    metsToAddC.mets         = {sprintf('%s_translocon_r', rxnName)};
    metsToAddC.metNames     = {sprintf('%s_translocon', rxnName)};
    metsToAddC.compartments = {'r'};

    model = addCatalyzedRxns(model, metsToAddC, rxnStructC, rxnName);
    metC = sprintf('%s_translocon[r]', rxnName);

    % --- Step D: BiP/NEF lumenal ratchet completes import (ATP cost
    % scales with chain length - see header note; ~1 ATP per ~40
    % residues threaded, Matlack/Mothes/Rapoport) ---
    n_ratchet = max(1, ceil(proteinLength / 40));
    rxnStructD.rxns     = {sprintf('%s_bip_ratchet', rxnName)};
    rxnStructD.rxnNames = rxnStructD.rxns;
    rxnStructD.lb       = [0];
    rxnStructD.mets = {
        metC, ...
        'ATP[r]', ...
        'H2O[r]', ...
        'ADP[r]', ...
        'Phosphate[r]', ...
        'H+[r]', ...
        sprintf('%s_er[r]', rxnName)
    };
    rxnStructD.stoichCoeffs = [-1, -n_ratchet, -n_ratchet, n_ratchet, n_ratchet, n_ratchet, 1];
    rxnStructD.catBy        = {'sec_BIP_NEFS_complex'};
    rxnStructD.eccodes      = {''};

    metsToAddD.mets         = {sprintf('%s_er_r', rxnName)};
    metsToAddD.metNames     = {sprintf('%s_er', rxnName)};
    metsToAddD.compartments = {'r'};

    model = addCatalyzedRxns(model, metsToAddD, rxnStructD, rxnName);
    startMet = sprintf('%s_er[r]', rxnName);
end

end
