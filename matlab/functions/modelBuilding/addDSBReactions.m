function [model, startMet] = addDSBReactions(model, proteinInfo, rxnName, startMet)
% ADDDSBREACTIONS Adds INTRACHAIN disulfide bond formation reaction in the ER lumen.
%
% This function models oxidative folding of a single polypeptide chain via
% disulfide bond formation, catalyzed by the ER-resident complex
% PDI1/ERV2/Ero1p (sec_PDI1_ERV2_Ero1p_complex - gene set corrected July
% 2026, see expandGEMv6.m audit notes: GQ67_05237/PDI1, GQ67_02568/ERO1,
% GQ67_03397/ERV2). The ATP cost is an aggregate literature estimate per
% disulfide bond and is NOT coupled to the glutathione redox couple.
%
% SCOPE NOTE (added July 2026, shadow_exploration project): this function
% is called once per polypeptide chain - for a single-chain protein, or
% for EACH chain (HC, LC) of a mAb individually via
% generatemAbSynthesisReactions.m - so proteinInfo.numDisulfide here
% should represent only that chain's own INTRAchain bonds. mAb
% INTERCHAIN bonds (joining separate HC/LC polypeptides, e.g. the IgG1
% hinge and HC-LC bonds) are added separately by
% addInterchainDSBReaction.m, called from generatemAbSynthesisReactions.m
% between chain assembly and dimerization, because interchain bonds are
% mechanistically and empirically distinct: they are ~100-1000x more
% susceptible to reduction than intrachain cysteines due to solvent
% exposure (Zhang et al. 2015, mAbs, "Structural and functional impact of
% chromatography-based interchain disulfide bond reduction," PMC10741012
% for the heterogeneity/solvent-accessibility finding), and their
% oxidation state is documented to be sensitive to extracellular/process
% glutathione and thioredoxin redox potential in a way plain intrachain
% folding is not (Kao et al. 2010, Handlogten et al. 2020, mAbs
% bioprocessing literature - see shadow_exploration Trastuzumab/gthrd_r
% review, July 2026). Splitting the two lets the model represent the
% mAb-specific interchain redox dependency that a flat per-protein
% disulfide count could not.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - structure with field numDisulfide (number of INTRACHAIN
%                 disulfide bonds in this chain)
%   rxnName     - base name of the protein/reaction (e.g., GQ67_XXXX)
%   startMet    - input peptide metabolite with [r] suffix
%
% OUTPUTS:
%   model       - updated model with disulfide bond reaction
%   startMet    - updated metabolite ID (now oxidized form)

% ------------------------------
% Define stoichiometry
% Based on:
% - 1 mol protein input
% - ~2.775 mol ATP per disulfide bond (aggregate estimate from literature)
% Ref: Liebermeister et al. 2010 (PNAS), BioNumbers BNID 111770
% ------------------------------

rxnStruct.rxns       = sprintf('%s_disulfide', rxnName);
rxnStruct.rxnNames   = rxnStruct.rxns;
rxnStruct.mets       = {
    startMet, ...
    'H2O[r]', ...
    'ATP[r]', ...
    'H+[r]', ...
    'Phosphate[r]', ...
    'ADP[r]', ...
    sprintf('%s_dsb[r]', startMet(1:end-3))
};
rxnStruct.lb        = [0];
n = 2.775 * proteinInfo.numDisulfide;
rxnStruct.stoichCoeffs = [-1, -n, -n, n, n, n, 1];

% ------------------------------
% Catalyst complex for disulfide bond formation
% Ref: Tu & Weissman 2002 (Cell), ER folding pathway
% ------------------------------
rxnStruct.catBy = 'sec_PDI1_ERV2_Ero1p_complex';

% ------------------------------
% EC number for disulfide isomerase (PDI1): EC 5.3.4.1
% Note: Ero1p acts as the oxidizing partner, not EC-classified
% ------------------------------
rxnStruct.eccodes = {'5.3.4.1'};

% ------------------------------
% Create the oxidized peptide metabolite
% ------------------------------
metsToAdd.mets         = {sprintf('%s_dsb_r', startMet(1:end-3))};
metsToAdd.metNames     = {sprintf('%s_dsb', startMet(1:end-3))};
metsToAdd.compartments = {'r'};

% ------------------------------
% Add to model
% ------------------------------
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Return updated metabolite for downstream use
startMet = sprintf('%s_dsb[r]', startMet(1:end-3));
end
