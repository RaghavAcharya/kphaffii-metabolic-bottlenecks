function [model, startMet] = addOglyERReactions(model, proteinInfo, rxnName, startMet)
% ADDOGLYERREACTIONS Adds ER-localized O-linked glycosylation to a protein.
%
% This function simulates the initial mannosylation of Ser/Thr residues
% on nascent polypeptides in the ER lumen using Dolichyl-phosphate-mannose.
%
% INPUTS:
%   model       - GEM model structure
%   proteinInfo - structure with field numOGlyc (number of O-glyc sites)
%   rxnName     - base reaction name (e.g., GQ67_XXXX)
%   startMet    - starting peptide metabolite in [r]
%
% OUTPUTS:
%   model       - updated model
%   startMet    - updated O-glycosylated peptide

% ------------------------------
% Reaction represents initial mannosylation by PMT family (Pmt1–7)
% Ref: Gentzsch & Tanner, Biochim Biophys Acta, 1996; EC 2.4.1.221
% ------------------------------

rxnStruct.rxns     = sprintf('%s_ogly_er', rxnName);
rxnStruct.rxnNames = rxnStruct.rxns;

% Stoichiometry:
% One mannose from Dol-P-Man added per glycosylation site
% Dol-P regenerated
% ------------------------------
n = proteinInfo.numOGlyc;
rxnStruct.mets = {
    startMet, ...
    'Dolichyl phosphate D-mannose[r]', ...
    'Dolichyl phosphate[r]', ...
    sprintf('%s_ogM1[r]', startMet(1:end-3))
};
rxnStruct.stoichCoeffs = [-1, -n, n, 1];
rxnStruct.lb = [0];

% Catalyst: PMT complex (includes Pmt1p–Pmt7p)
% Ref: Strahl-Bolsinger et al., 1999 (Yeast)
rxnStruct.catBy = {'sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_Pmt3p_complex'};

% EC for O-mannosyltransferase
rxnStruct.eccodes = {'2.4.1.221'};

% ------------------------------
% Add new metabolite: mono-O-mannosylated form
% ------------------------------
metsToAdd.mets         = {sprintf('%s_ogM1_r', startMet(1:end-3))};
metsToAdd.metNames     = {sprintf('%s_ogM1', startMet(1:end-3))};
metsToAdd.compartments = {'r'};

% ------------------------------
% Add to model
% ------------------------------
model = addCatalyzedRxns(model, metsToAdd, rxnStruct, rxnName);

% Update final metabolite
startMet = sprintf('%s_ogM1[r]', startMet(1:end-3));
end
