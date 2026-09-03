clear
close all
clc
% === Portable setup (added for GitHub release) ===
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(fullfile(repoRoot,'matlab','functions')));
cd(scriptDir);                              % colocated inputs/outputs resolve here
run(fullfile(repoRoot,'setup_solver.m'));   % edit setup_solver.m for your Gurobi/COBRA/RAVEN paths
% =================================================
%% Import base iMT1026v3 model
iMT1026v3 = load(fullfile(repoRoot,'models','iMT1026v3.mat'));
iMT1026v3 = iMT1026v3.iMT1026v3;
%% Add Necessary Metabolites and Reactions

% Excel file consisting of two sheets - one with new metabolites and one
% with the corresponding necessary reactions. 
% Define the file name
fileName = fullfile(repoRoot,'data','baselineModifications.xlsx');
% Read the "New Reactions" sheet into the variable newRxns
newRxns = readtable(fileName, 'Sheet', 'New Reactions');
% Read the "New Metabolites" sheet into the variable newMets
newMets = readtable(fileName, 'Sheet', 'New Metabolites');

metsToAdd.mets = newMets.MetName;
metsToAdd.metNames = newMets.MetName_1;
metsToAdd.compartments = newMets.Compartment;
iMT1026v4 = addMets(iMT1026v3,metsToAdd);

rxnsToAdd.rxns = newRxns.RxnID';
rxnsToAdd.rxnNames = rxnsToAdd.rxns;
rxnsToAdd.equations = newRxns.Equation';
rxnsToAdd.rRules = newRxns.GPR';
iMT1026v4 = addRxns(iMT1026v4,rxnsToAdd,3,[],true,true);

numMetRxns = length(iMT1026v4.rxns);
ppTable = readtable(fullfile(repoRoot,'data','ppPSIMv3_RP.csv'));
%% Determine costs associated with production/translation/etc of each protein
% Set up Carbon Source Data 
cSourceData.Name = ["Citrate", "NoCofeed", "Rhamnose", "Sorbitol", "Glycerol", "Xylose"];
cSourceData.ExrxnName = ["Ex_cit", "Ex_meoh", "Ex_rmn", "Ex_sbt_D", "Ex_glyc", "Ex_xyl_D"];
iMT1026v3_constr = iMT1026v4;
[exRxns, exRxnIds] = getExchangeRxns(rmfield(iMT1026v3_constr,"unconstrained"));
cSourceData.ExrxnIdx = cell(size(cSourceData.ExrxnName));
iMT1026v3_constr = iMT1026v4;
for k = 1:length(cSourceData.ExrxnName)
    cSourceData.ExrxnIdx{k} = find(strcmp(iMT1026v3_constr.rxns, cSourceData.ExrxnName{k}));
end

numGenes = height(ppTable);
% Set up Matrices to store synth/trafficking cost ratios

% Preallocate synthCost with the appropriate size
uptakeRxns = {'Ex_h2o','Ex_nh4','Ex_pi','Ex_h','Ex_o2','Ex_so4', 'Ex_fe2', 'Ex_btn'}';
necUptakeRxnIds = cell(size(uptakeRxns));
for k = 1:length(necUptakeRxnIds)
    necUptakeRxnIds{k} = find(strcmp(iMT1026v3_constr.rxns, uptakeRxns{k}));
end

synthCost = zeros(numGenes, length(cSourceData.Name), length(necUptakeRxnIds)+3); %gene x csource x flux from (nec uptake molecules, carbon source, protein, biomass)
synthCostBiomass = synthCost;

for c = 1:length(cSourceData.Name)
    tic
    sprintf("Carbon Source: %s", cSourceData.Name{c})
    % Return to 'minimal media' conditions
    iMT1026v3_constr.lb(exRxnIds) = 0;
    iMT1026v3_constr = setParam(iMT1026v3_constr,'lb',uptakeRxns,zeros(length(uptakeRxns),1)-1000);
    % iMT1026v3_constr = setParam(iMT1026v3_constr,'lb','Ex_o2',-10);

    % Open exchange reaction for specific csource of interest
    iMT1026v3_constr = setParam(iMT1026v3_constr,'lb',cSourceData.ExrxnName{c},-10); 
    
    % Establish baseline biomass load
    iMT1026v3_constr = setParam(iMT1026v3_constr,'lb', 2041, -1000);
    iMT1026v3_constr = setParam(iMT1026v3_constr,'ub', 2041, 1000);
    iMT1026v3_constr = setParam(iMT1026v3_constr,'obj', 2041, 1);
    sol_biomass=solveLP(iMT1026v3_constr);
    
    
    for i = 1:numGenes
        if mod(i,500)==0
            sprintf("Currently on Gene %i", i)
        end
        proteinInfo.aaSequence = ppTable.Sequence{i};
        proteinInfo.dnaSequence = ppTable.dnaSequence{i};
        proteinInfo.throughER = ppTable.ThroughER(i);
        proteinInfo.numNGlyc = ppTable.NumNGlyc(i);
        proteinInfo.numDisulfide = ppTable.NumDisulfides(i);
        proteinInfo.transmembrane = ppTable.Transmembrane(i);
        proteinInfo.numOGlyc = ppTable.NumOGlyc(i);
        proteinInfo.GPIAddition = ppTable.GPIsite(i);
        proteinInfo.finalLocation = strsplit(ppTable.Localization{i},',');
        
        warning ('off','all');
        [printOut, iMT1026v4] = evalc('generateSynthesisReactions4(iMT1026v3_constr,ppTable.Gene{i},proteinInfo)');
        % Find synthesis reaction, run FBA
        iMT1026v4 = setParam(iMT1026v4,'lb', 2041, -1000);
        synthRxnIdx = find(strcmp(iMT1026v4.rxnNames, sprintf("%s_maturation",ppTable.Gene{i})));

        %iMT1026v4 = setParam(iMT1026v4,'ub', synthRxnIdx, 1);
        iMT1026v4 = setParam(iMT1026v4,'obj', synthRxnIdx, 1);
        sol_syn=solveLP(iMT1026v4);
        
        if ~isempty(sol_syn.x)
            for k = 1:length(necUptakeRxnIds)
                synthCost(i,c,k) = sol_syn.x(necUptakeRxnIds{k});
            end
            
            synthCost(i,c,length(necUptakeRxnIds)+1) = sol_syn.x(cSourceData.ExrxnIdx{c});
            synthCost(i,c,length(necUptakeRxnIds)+2) = sol_syn.x(synthRxnIdx);
            synthCost(i,c,length(necUptakeRxnIds)+3) = sol_syn.x(2041);
        end
        
        iMT1026v4 = setParam(iMT1026v4,'lb', 2041, -0.1*sol_biomass.f);
        sol_synbiomass=solveLP(iMT1026v4);
        if ~isempty(sol_synbiomass.x)
            for k = 1:length(necUptakeRxnIds)
                synthCostBiomass(i,c,k) = sol_synbiomass.x(necUptakeRxnIds{k});
            end
            synthCostBiomass(i,c,length(necUptakeRxnIds)+1) = sol_synbiomass.x(cSourceData.ExrxnIdx{c});
            synthCostBiomass(i,c,length(necUptakeRxnIds)+2) = sol_synbiomass.x(synthRxnIdx);
            synthCostBiomass(i,c,length(necUptakeRxnIds)+3) = sol_synbiomass.x(2041);
        end
    end
    toc
end
%%
synthCostTable = zeros(numGenes, length(cSourceData.Name)*(length(necUptakeRxnIds)+3));
synthCostBiomassTable = synthCostTable;
for k = 1:(length(necUptakeRxnIds)+3)
    synthCostTable(:,((k-1)*length(cSourceData.Name)+1):k*length(cSourceData.Name)) = synthCost(:,:,k);
    synthCostBiomassTable(:,((k-1)*length(cSourceData.Name)+1):k*length(cSourceData.Name)) = synthCostBiomass(:,:,k);
end
scT = table(synthCostTable, RowNames = ppTable.Gene);
scTB = table(synthCostBiomassTable, RowNames = ppTable.Gene);
writetable(scT, 'synthCost_Climited.csv')
writetable(scTB, 'synthCostBiomas_Climited.csv')