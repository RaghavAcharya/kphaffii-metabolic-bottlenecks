function [model] = addTransportReactions(model, proteinInfo, rxnName, maturationMet)
% ADDTRANSPORTREACTIONS Adds transport and dilution reactions for mature proteins
% based on their final cellular destination(s).
%
% INPUTS:
%   model         - GEM structure
%   proteinInfo   - struct with .finalLocation (cell array of letters)
%   rxnName       - gene name (e.g. GQ67_XXXX)
%   maturationMet - struct with .metNames (cell), .compartments (cell)
%
% OUTPUT:
%   model         - updated GEM

% Extract base metabolite name and source compartment
srcMetName = maturationMet.metNames{1};      % unwrap from cell
srcComp = maturationMet.compartments{1};
srcMet = sprintf('%s[%s]', srcMetName, srcComp);

% Ensure finalLocation is a cell array
if ischar(proteinInfo.finalLocation)
    proteinInfo.finalLocation = {proteinInfo.finalLocation};
end


if length(proteinInfo.finalLocation)==1 
    if isempty(proteinInfo.finalLocation{1})
        proteinInfo.finalLocation = {'c'};
    end
elseif length(proteinInfo.finalLocation)>1
    if any(contains(proteinInfo.finalLocation,'t')) && any(contains(proteinInfo.finalLocation,'e'))
        proteinInfo.finalLocation(contains(proteinInfo.finalLocation,'e')) = [];
    end
end

for i = 1:length(proteinInfo.finalLocation)
    originalDest = proteinInfo.finalLocation{i};

    % Redirect 't' to 'e' for modeling, but preserve catalyst logic
    if strcmp(originalDest, 't')
        dest = 'e';  % Use extracellular compartment
        isTransmembrane = true;
    else
        dest = originalDest;
        isTransmembrane = false;
    end

    destMet = sprintf('%s[%s]', srcMetName, dest);
    baseID = srcMetName;

    % Default GTP-based transport logic
    transRxn.rxnID       = sprintf('%s_to%s', rxnName, dest);
    transRxn.metabolites = {srcMet, 'H2O[c]', 'GTP[c]', 'H+[c]', 'Phosphate[c]', 'GDP[c]', destMet};
    transRxn.stoich      = [-1, -1, -1, 1, 1, 1, 1];
    transRxn.lb          = 0;
    transRxn.catBy       = '';
    transRxn.EC          = '';

    % Dilution (sink) reaction
    diluteRxn.rxnID       = sprintf('%s_dilute_%s', rxnName, dest);
    diluteRxn.metabolites = {destMet};
    diluteRxn.stoich      = [-1];
    diluteRxn.lb          = 0;

    switch dest
        case {'c'}
            % Cytosolic proteins: only dilute if already in cytoplasm
            if srcComp == 'c'
                transRxn = [];  % No transport
            else
                % Example: re-imported protein (unusual case)
                transRxn.catBy = '';
            end

        case {'g'}
            % Golgi via COPII vesicles: Sar1/Sec23–24/Sec13–31
            % Ref: Barlowe et al. (1994), Cell
            if srcComp ~= 'g'
                transRxn.catBy = 'sec_COPII_general_complex';
                transRxn.stoich = [-1, -1, -1, 1, 1, 1, 1]; % 1 GTP
            else
                transRxn = [];
            end

        case 'e'
            if isTransmembrane
                % Transmembrane protein via LDSV → extracellular face
                % Ref: TerBush et al. (1996), J Cell Biol
                transRxn.catBy = 'sec_LDSV_complex';
                transRxn.stoich = [-1, -1, -1, 1, 1, 1, 1];  % 1 GTP
            else
                % Secreted via HDSV vesicles
                % Ref: Schekman et al. (1997), PNAS
                transRxn.catBy = 'sec_HDSV_complex';
                transRxn.stoich = [-1, -2, -2, 2, 2, 2, 1];  % 2 GTP
            end
        case {'r'}
            % ER targeting (retrograde): COPI
            % Ref: Lee & Goldberg (2010), Nat Rev Mol Cell Biol
            transRxn.catBy = 'sec_COPI_complex';
            transRxn.stoich = [-1, -2, -2, 2, 2, 2, 1];  % 2 GTP

        case {'n'}
            % Nucleus: importin-α/β and RanGTP
            % Ref: Görlich & Kutay (1999), Annu Rev Cell Dev Biol
            transRxn.catBy = 'sec_Importin_Ran_complex';
            transRxn.stoich = [-1, -1, -1, 1, 1, 1, 1];  % 1 GTP

        case {'m'}
            % Mitochondria: TOM/TIM complex
            % Ref: Neupert & Herrmann (2007), Annu Rev Biochem
            transRxn.catBy = 'sec_TOM_TIM_complex';
            transRxn.stoich = [-1, -1, -1, 1, 1, 1, 1];  % ATP + GTP cost

        case {'i'}
            % Mitochondrial intermembrane space
            transRxn.catBy = 'sec_TOM_complex';

        case {'v'}
            % Vacuole via CPY pathway
            % Ref: Bowers & Stevens (2005), JCB
            transRxn.catBy = 'sec_CPY_complex_combined';
            transRxn.metabolites = {srcMet, 'H2O[c]', 'ATP[c]', 'GTP[c]', 'H+[c]', ...
                                    'Phosphate[c]', 'ADP[c]', 'GDP[c]', destMet};
            transRxn.stoich = [-1, -5, -1, -4, 5, 5, 1, 4, 1];

        case {'vm'}
            % Vacuole membrane via ALP pathway
            % Ref: Reggiori & Pelham (2001), JCB
            transRxn.catBy = 'sec_ALP_complex';
            transRxn.stoich = [-1, -4, -4, 4, 4, 4, 1];

        case {'x'}
            % Peroxisome import via Pex machinery
            % Ref: Platta et al. (2005), Biochim Biophys Acta
            transRxn.catBy = 'sec_PEX_import_complex';
            transRxn.stoich = [-1, -2, -2, 2, 2, 2, 1];

        otherwise
            warning('Unknown destination code: %s', dest);
            continue;
    end

    % Add destination metabolite to model
    metStruct.mets         = {sprintf('%s_%s', baseID, dest)};
    metStruct.metNames     = {baseID};
    metStruct.compartments = {dest};

    if ~isempty(transRxn)
        % Add both transport and dilution
        rxnStruct.rxns         = {transRxn.rxnID, diluteRxn.rxnID};
        rxnStruct.rxnNames     = rxnStruct.rxns;
        rxnStruct.mets         = {transRxn.metabolites, diluteRxn.metabolites};
        rxnStruct.stoichCoeffs = {transRxn.stoich, diluteRxn.stoich};
        rxnStruct.catBy        = {transRxn.catBy, ''};
        rxnStruct.lb           = [transRxn.lb, diluteRxn.lb];

        model = addCatalyzedRxns(model, metStruct, rxnStruct, rxnName);
    else
        % Add only dilution
        rxnStruct.rxns         = {diluteRxn.rxnID};
        rxnStruct.rxnNames     = rxnStruct.rxns;
        rxnStruct.mets         = {diluteRxn.metabolites};
        rxnStruct.stoichCoeffs = {diluteRxn.stoich};
        rxnStruct.catBy        = {''};
        rxnStruct.lb           = diluteRxn.lb;

        model = addRxns(model, rxnStruct, 3, [], true, true);
    end
end
end
