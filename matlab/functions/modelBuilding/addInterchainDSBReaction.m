function [model, startMet] = addInterchainDSBReaction(model, numInterchainDisulfide, rxnName, startMet)
% ADDINTERCHAINDSBREACTION Adds INTERCHAIN disulfide bond formation
% reaction for multi-chain proteins (e.g. mAb heavy+light chain
% assembly), coupled to the ER glutathione redox couple.
%
% Added July 2026 (shadow_exploration project) alongside the
% expandGEMv6.m GPR fixes. Context: the phase-5 shadow-price pilot found
% that ER glutathione regeneration (gthrd_r) ranked as the #1 internal
% dual for all 10 pilot proteins under all 3 limiting-nutrient axes, but
% traced this to a structural artifact rather than genuine biology - no
% reaction in the reconstruction actually coupled a protein's disulfide
% burden to glutathione at all (addDSBReactions.m costs disulfides in
% ATP only), and the mAb assembly/dimerization steps in
% generatemAbSynthesisReactions.m were completely free (uncatalyzed, no
% cofactor, no gene). This meant the model had no mechanism by which it
% could ever predict a differential, mAb-specific glutathione effect -
% exactly the gap that showed up when bench data (glutathione dosing
% helping Trastuzumab but not most other pilot proteins) had no
% counterpart in the model's own predictions.
%
% This function is the fix: it adds a SEPARATE reaction for interchain
% bond formation (called once per assembled multi-chain intermediate, not
% per individual chain - see generatemAbSynthesisReactions.m, inserted
% between the HC+LC assembly step and the dimerization step), consuming
% oxidized glutathione (gthox_r) and producing reduced glutathione
% (gthrd_r) in the textbook 1:2 ratio for thiol-disulfide exchange
% (1 GSSG oxidizes one pair of free cysteine thiols to 1 disulfide bond,
% releasing 2 GSH - standard redox stoichiometry, see e.g. Gilbert 1990,
% Adv Enzymol Relat Areas Mol Biol, for the general thiol-exchange
% chemistry). This now matches the model's native PDI_reoxidation_GSSG
% reaction (baselineModifications.xlsx), which used the same GSSG-consumed
% : GSH-produced roles but an inverted 2:1 ratio - and, as found in the
% same July 2026 audit, had its PDI_mature/PDI_mature_ox reactant/product
% roles swapped as well (the equation as originally written had two
% simultaneous reductions - PDI_ox->PDI_reduced and GSSG->GSH - with no
% compensating oxidation, which isn't balanceable at any ratio). Both
% reactions now consistently model reoxidation (reduced->oxidized,
% consuming an electron acceptor) with the correct 1 GSSG : 2 GSH ratio -
% see Cuozzo & Kaiser 1999, Nat Cell Biol, PMID 10559898, for the GSSG-
% mediated PDI reoxidation pathway itself.
%
% Catalyzed by sec_PDI1_ERV2_Ero1p_complex, whose gene set was corrected
% in the same audit (complexData.xlsx already had the right genes -
% GQ67_05237/PDI1, GQ67_02568/ERO1, GQ67_03397/ERV2 - expandGEMv6.m just
% wasn't copying them into model.catalysts; see resolveComplexGeneSet in
% expandGEMv6.m). Without that fix this reaction would be enzyme-capacity
% -unconstrained the same way PDI_reoxidation_GSSG currently is.
%
% Standard human IgG1 architecture (for reference when populating
% numInterchainDisulfide per protein): 4 interchain disulfide bonds total
% (2 in the hinge, H-H; 1 per H-L pair, x2) plus 12 intrachain disulfides
% across the 12 Ig domains - see Liu & May 2012, mAbs 4(1):17-23,
% "Disulfide bond structures of IgG molecules"
% (https://doi.org/10.4161/mabs.4.1.18347). This count is NOT assumed by
% this function - it must be supplied by the caller (see
% generatemAbSynthesisReactions.m and buildShadowExplorationModels.m) so
% existing protein specs are unaffected until deliberately updated.
%
% INPUTS:
%   model                  - GEM model structure
%   numInterchainDisulfide - number of interchain disulfide bonds formed
%                            at this assembly step (0 => no reaction is
%                            added, and startMet is returned unchanged)
%   rxnName                - base name for the reaction/metabolite IDs
%                            (e.g. 'trastuzumab')
%   startMet                - input assembled-intermediate metabolite ID
%                            with [r] suffix (e.g.
%                            'trastuzumab_LC_HC[r]', free-thiol form)
%
% OUTPUTS:
%   model    - updated model with the interchain disulfide reaction
%   startMet - updated metabolite ID (disulfide-bonded form), or the
%              unchanged input if numInterchainDisulfide == 0

if nargin < 2 || isempty(numInterchainDisulfide) || numInterchainDisulfide == 0
    % Backward-compatible no-op: existing protein specs (built before
    % this field existed) are unaffected until numInterchainDisulfide is
    % explicitly supplied.
    return
end

bondedMet = sprintf('%s_dsb[r]', startMet(1:end-3));

rxnStruct.rxns       = sprintf('%s_interchain_disulfide', rxnName);
rxnStruct.rxnNames   = rxnStruct.rxns;
rxnStruct.mets       = {
    startMet, ...
    'Oxidized glutathione[r]', ...
    bondedMet, ...
    'Reduced glutathione[r]'
};
rxnStruct.lb        = 0;
n = numInterchainDisulfide;
% 1 GSSG oxidizes one free-thiol pair to 1 disulfide bond, releasing 2
% GSH (Sevier & Kaiser 2002) - see header note above re: the model's
% existing PDI_reoxidation_GSSG reaction using a different ratio.
rxnStruct.stoichCoeffs = [-1, -n, 1, 2*n];

% Catalyst complex for interchain disulfide bond formation - same
% PDI1/ERV2/Ero1p machinery as intrachain folding (addDSBReactions.m),
% now gene-linked (see expandGEMv6.m audit notes above).
rxnStruct.catBy = 'sec_PDI1_ERV2_Ero1p_complex';
rxnStruct.eccodes = {'5.3.4.1'};

metsToAdd.mets         = {sprintf('%s_dsb_r', startMet(1:end-3))};
metsToAdd.metNames     = {sprintf('%s_dsb', startMet(1:end-3))};
metsToAdd.compartments = {'r'};

model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

startMet = bondedMet;
end
