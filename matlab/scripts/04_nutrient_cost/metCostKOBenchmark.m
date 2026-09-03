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
%% Determine potential improvement in productivity when removing load of each protein

% Set uptake rxn limits 
uptakeRxns = {'Ex_h2o','Ex_nh4','Ex_pi','Ex_h','Ex_o2','Ex_so4', 'Ex_fe2', 'Ex_btn'}';
iMT1026v4 = setParam(iMT1026v4,'lb',uptakeRxns,zeros(length(uptakeRxns),1)-1000);
iMT1026v4 = setParam(iMT1026v4,'lb',"Ex_meoh",-1);

% Add recombinant protein production reactions
rpIdx = 5179;
rpproteinInfo.aaSequence = ppTable.Sequence{rpIdx};
rpproteinInfo.dnaSequence = ppTable.dnaSequence{rpIdx};
rpproteinInfo.throughER = ppTable.ThroughER(rpIdx);
rpproteinInfo.numNGlyc = ppTable.NumNGlyc(rpIdx);
rpproteinInfo.numDisulfide = ppTable.NumDisulfides(rpIdx);
rpproteinInfo.transmembrane = ppTable.Transmembrane(rpIdx);
rpproteinInfo.numOGlyc = ppTable.NumOGlyc(rpIdx);
rpproteinInfo.GPIAddition = ppTable.GPIsite(rpIdx);
rpproteinInfo.finalLocation = strsplit(ppTable.Localization{rpIdx},',');
iMT1026v4 = generateSynthesisReactions4(iMT1026v4,ppTable.Gene{rpIdx},rpproteinInfo);
rpProdIdx = 2267;
            
% Set up tables for information storage 
numSteps = 20;
fluxData = zeros(length(iMT1026v4.rxns),5120,5); % Full flux dataset will be stored in case it can be mined later.
rpProductionFluxData = zeros(5120,numSteps); % Subset of flux dataset containing only the production fluxes of RP under different loads.
geneProductionFluxData = zeros(5120,numSteps+1); % Subset of flux dataset containing only the production fluxes of each gene tested. The last column is the maximum production flux of the gene.

% Run FBA w biomass accumulation maximization
iMT1026v4 = setParam(iMT1026v4,'obj', 2041, 1);
sol_biomass=solveLP(iMT1026v4,3);

% Constrain model to have lb on biomass accumulation, check max RP
% productivity
% iMT1026v4 = setParam(iMT1026v4,'lb', 2041, sol_biomass.x(2041)*0.1);
iMT1026v4 = setParam(iMT1026v4,'obj', rpProdIdx, 1);
sol_RP=solveLP(iMT1026v4,3);

%% Loop through each protein 
reverseStr = '';
for i = 1:5120
    % Create model with protein production reactions
    proteinInfo.aaSequence = ppTable.Sequence{i};
    proteinInfo.dnaSequence = ppTable.dnaSequence{i};
    proteinInfo.throughER = ppTable.ThroughER(i);
    proteinInfo.numNGlyc = ppTable.NumNGlyc(i);
    proteinInfo.numDisulfide = ppTable.NumDisulfides(i);
    proteinInfo.transmembrane = ppTable.Transmembrane(i);
    proteinInfo.numOGlyc = ppTable.NumOGlyc(i);
    proteinInfo.GPIAddition = ppTable.GPIsite(i);
    proteinInfo.finalLocation = strsplit(ppTable.Localization{i},',');
    [printOut, iMT1026v4_i] = evalc('generateSynthesisReactions4(iMT1026v4,ppTable.Gene{i},proteinInfo)');

    % Run FBA for native protein maximization
    synthRxnIdx = find(strcmp(iMT1026v4_i.rxnNames, sprintf("%s_maturation",ppTable.Gene{i})));
    iMT1026v4_i = setParam(iMT1026v4_i,'obj', synthRxnIdx, 1);
    [printOut, sol_syn] = evalc('solveLP(iMT1026v4_i, 3)');
    geneProductionFluxData(i, end) = sol_syn.x(synthRxnIdx);
    
    % Set new objective function to be RP production
    iMT1026v4_i = setParam(iMT1026v4_i, 'obj', rpProdIdx, 1);
    % Loop through multiple lower bounds native protein maximization - 0%
    % to ~90% of the maximum of native protein production
    synthFluxes = linspace(0,0.99*sol_syn.x(synthRxnIdx), numSteps);
    for j = 1:numSteps
        % Run FBA for recombinant protein maximization
        iMT1026v4_i_j = setParam(iMT1026v4_i,'lb', synthRxnIdx, synthFluxes(j));
        [printOut, sol_cost] = evalc('solveLP(iMT1026v4_i_j)');

        % Store flux results
        fluxData(:,i,j) = sol_cost.x(1:length(iMT1026v4.rxns));
        rpProductionFluxData(i, j) = sol_cost.x(rpProdIdx);
        geneProductionFluxData(i, j) = sol_cost.x(synthRxnIdx);
    end

   percentDone = 100 * i / 5120;
   msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

%% Save the data 
save("fluxData.mat","fluxData")
writematrix(geneProductionFluxData, 'geneProductionFluxData_trastuzumabHC.csv')
writematrix(rpProductionFluxData, 'rpProductionFluxData_trastuzumabHC.csv')

%% Load data 
clearvars -except fluxData geneProductionFluxData rpProductionFluxData iMT1026v4 sol_RP sol_biomass ppTable numSteps
close all
clc 

load("fluxData.mat")
geneProductionFluxData = readmatrix('geneProductionFluxData_trastuzumabHC.csv');
rpProductionFluxData = readmatrix('rpProductionFluxData_trastuzumabHC.csv');
%% Initial analysis of data 
close all
figure 
hold on 
for i = 1:5120
    plot([rpProductionFluxData(i,:), 0], geneProductionFluxData(i,:),'Color',[0.75, 0.75, 0.75])
end
xlabel('Recombinant Protein Production')
ylabel('Gene Production')
hold off

%% Analysis of all Runs w




% %% Exploration of the weird ones
% strangeIdxs = find(rpProductionFluxData(:,1)./rpProductionFluxData(1,1)>(1.000001));
% strangeFluxes = fluxData(:,strangeIdxs,1:numSteps);
% DEstrangeFluxes = strangeFluxes;
% for i = 1:numSteps
%     for j = 1:length(strangeIdxs)
%         DEstrangeFluxes(:,j,i) = DEstrangeFluxes(:,j,i)-sol_RP.x;
%     end
% end
% DEfluxes = find(any(abs(DEstrangeFluxes(:,:,1))>10,2));
% baseDEfluxes = DEstrangeFluxes(:,:,1);
% iMT1026v4.rxnNames(DEfluxes)