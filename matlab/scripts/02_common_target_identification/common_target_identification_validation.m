%% HSA Production Model Analysis
% This script creates and analyzes a genome-scale metabolic model for
% Human Serum Albumin (HSA) production in Pichia pastoris and validates
% against experimental data.
% 
% Steps:
% 1. Import and configure base metabolic model
% 2. Add HSA production pathway
% 3. Validate model with experimental data
% 4. Analyze correlation between predicted and experimental productivity
% 5. Run FSEOF and identify reactions predicted to improve productivity

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
%% 1. Create Base Model
% Add Necessary Metabolites and Reactions
% Excel file consisting of two sheets - one with new metabolites and one
% with the corresponding necessary reactions. 
% Define the file name

% Import base iMT1026v3 model
fprintf('Importing base metabolic model...\n');
iMT1026v3 = importModel(fullfile(repoRoot,'models','mbt212871sup0003AppendixS3.xml'), false); % Obtained from https://www.ebi.ac.uk/biomodels/MODEL1612130000#Files
iMT1026v3.rev = double(iMT1026v3.lb < 0);  % Set reversibility flag based on lower bounds
iMT1026v3 = rmfield(iMT1026v3, "unconstrained");  % Remove unconstrained field


fprintf('Modifying base metabolic model...\n');
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

% Switch biomass reaction to more general biomass reaction
% Less accurate for glycerol specifically, but better for different carbon sources
iMT1026v4.ub(2233) = 1000;  % Disable specific biomass reaction
iMT1026v4.ub(1305) = 000;  % Enable general biomass reaction

iMT1026v4_original = iMT1026v4;

numMetRxns = length(iMT1026v4.rxns);
%% 2. 
% Define carbon sources and exchange reactions
carbonExRxns = {'Ex_sbt_D', 'Ex_arab_D', 'Ex_fru', 'Ex_gal', 'Ex_sucr', 'Ex_xyl_D', 'Ex_rmn', 'Ex_cit', 'Ex_glyc'};
carbons = {'Sorbitol', 'D-Arabinose', 'Fructose', 'Galactose', 'Sucrose', 'D-Xylose', 'Rhamnose', 'Citrate', 'Glycerol'};
[~, carbonExRxnsIdxs] = ismember(carbonExRxns, iMT1026v4_original.rxns);

metadata = readtable('carbonScreen_metadata.xlsx');
RPs = unique(metadata.Molecule);

ppTable = readtable(fullfile(repoRoot,'data','ppPSIMv3_RP.csv'));

% Calculate growth rate from experimental data (cellgrowth-5)/24/7 in 1/hr
metadata.GrowthRate = (metadata.Cellgrowth - 0) / 48 / 7;
metadata.ProdRateEst = zeros(height(metadata), 1);
metadata.timeToConsumption = zeros(height(metadata),1);
metadata.CofeedConsumptionRate = zeros(height(metadata), 1);
metadata.MeOHConsumptionRate = zeros(height(metadata), 1);
models = {};
for RP = 1:length(RPs)
    if strcmp(RPs{RP},'trastuzumab')

        RP_mass = 0;
        %With PTMs

        [~, RP_idx] = ismember('trastuzumab_LC', ppTable.Gene);
        proteinInfoLC.aaSequence = ppTable.Sequence{RP_idx};
        proteinInfoLC.dnaSequence = ppTable.dnaSequence{RP_idx};
        proteinInfoLC.throughER = ppTable.ThroughER(RP_idx);
        proteinInfoLC.numNGlyc = ppTable.NumNGlyc(RP_idx);
        proteinInfoLC.numDisulfide = ppTable.NumDisulfides(RP_idx);
        proteinInfoLC.transmembrane = ppTable.Transmembrane(RP_idx);
        proteinInfoLC.numOGlyc = ppTable.NumOGlyc(RP_idx);
        proteinInfoLC.GPIAddition = ppTable.GPIsite(RP_idx);
        proteinInfoLC.finalLocation = strsplit(ppTable.Localization{RP_idx},',');
        RP_mass = 2*ppTable.proteinSize(RP_idx)+RP_mass; % FI

        [~, RP_idx] = ismember('trastuzumab_HC', ppTable.Gene);
        proteinInfoHC.aaSequence = ppTable.Sequence{RP_idx};
        proteinInfoHC.dnaSequence = ppTable.dnaSequence{RP_idx};
        proteinInfoHC.throughER = ppTable.ThroughER(RP_idx);
        proteinInfoHC.numNGlyc = ppTable.NumNGlyc(RP_idx);
        proteinInfoHC.numDisulfide = ppTable.NumDisulfides(RP_idx);
        proteinInfoHC.transmembrane = ppTable.Transmembrane(RP_idx);
        proteinInfoHC.numOGlyc = ppTable.NumOGlyc(RP_idx);
        proteinInfoHC.GPIAddition = ppTable.GPIsite(RP_idx);
        proteinInfoHC.finalLocation = strsplit(ppTable.Localization{RP_idx},',');
        RP_mass = 2*ppTable.proteinSize(RP_idx)+RP_mass; % FI
        iMT1026v4_RP = generatemAbSynthesisReactions(iMT1026v4,'trastuzumab',proteinInfoLC, proteinInfoHC);

        % No PTMs
        % [~, RP_idx] = ismember('trastuzumab_LC', ppTable.Gene);
        % proteinInfoLC.aaSequence = ppTable.Sequence{RP_idx};
        % proteinInfoLC.dnaSequence = ppTable.dnaSequence{RP_idx};
        % proteinInfoLC.throughER = 0;
        % proteinInfoLC.numNGlyc = 0;
        % proteinInfoLC.numDisulfide = 0;
        % proteinInfoLC.transmembrane = 0;
        % proteinInfoLC.numOGlyc = 0;
        % proteinInfoLC.GPIAddition = 0;
        % proteinInfoLC.finalLocation = strsplit(ppTable.Localization{RP_idx},',');
        % RP_mass = 2*ppTable.proteinSize(RP_idx)+RP_mass; % FI
        % 
        % [~, RP_idx] = ismember('trastuzumab_HC', ppTable.Gene);
        % proteinInfoHC.aaSequence = ppTable.Sequence{RP_idx};
        % proteinInfoHC.dnaSequence = ppTable.dnaSequence{RP_idx};
        % proteinInfoHC.throughER = 0;
        % proteinInfoHC.numNGlyc = 0;
        % proteinInfoHC.numDisulfide = 0;
        % proteinInfoHC.transmembrane = 0;
        % proteinInfoHC.numOGlyc = 0;
        % proteinInfoHC.GPIAddition = 0;
        % proteinInfoHC.finalLocation = strsplit(ppTable.Localization{RP_idx},',');
        % RP_mass = 2*ppTable.proteinSize(RP_idx)+RP_mass; % FI
        % iMT1026v4_RP = generatemAbSynthesisReactions(iMT1026v4,'trastuzumab',proteinInfoLC, proteinInfoHC);

        iMT1026v4_RP_original = iMT1026v4_RP;  % Store original model for reuse
        
    else
        % With PTMs
        [~, RP_idx] = ismember(RPs{RP}, ppTable.Gene);
        proteinInfo.aaSequence = ppTable.Sequence{RP_idx};
        proteinInfo.dnaSequence = ppTable.dnaSequence{RP_idx};
        proteinInfo.throughER = ppTable.ThroughER(RP_idx);
        proteinInfo.numNGlyc = ppTable.NumNGlyc(RP_idx);
        proteinInfo.numDisulfide = ppTable.NumDisulfides(RP_idx);
        proteinInfo.transmembrane = ppTable.Transmembrane(RP_idx);
        proteinInfo.numOGlyc = ppTable.NumOGlyc(RP_idx);
        proteinInfo.GPIAddition = ppTable.GPIsite(RP_idx);
        proteinInfo.finalLocation = strsplit(ppTable.Localization{RP_idx},',');

        % No PTMs
        % proteinInfo.aaSequence = ppTable.Sequence{RP_idx};
        % proteinInfo.dnaSequence = ppTable.dnaSequence{RP_idx};
        % proteinInfo.throughER = 0;
        % proteinInfo.numNGlyc = 0;
        % proteinInfo.numDisulfide = 0;
        % proteinInfo.transmembrane = 0;
        % proteinInfo.numOGlyc = 0;
        % proteinInfo.GPIAddition = 0;
        % proteinInfo.finalLocation = strsplit(ppTable.Localization{RP_idx},',');
        
        iMT1026v4_RP = generateSynthesisReactions4(iMT1026v4,ppTable.Gene{RP_idx},proteinInfo);
        iMT1026v4_RP_original = iMT1026v4_RP;  % Store original model for reuse
        
        RP_mass = ppTable.proteinSize(RP_idx);
    end
    %% 2a. Validate Model with Different Carbon Source
    iMT1026v4_RP = iMT1026v4_RP_original;
    iMT1026v4_RP = setParam(iMT1026v4_RP, 'obj', sprintf('%s_toe',RPs{RP}), 1);  % Set HSA synthesis as objective
    iMT1026v4_RP = setParam(iMT1026v4_RP, 'lb', 'Ex_glyc', 0);  % Zero glycerol uptake by default
    iMT1026v4_RP = setParam(iMT1026v4_RP, 'lb', 'Ex_tre', -0.1);
    iMT1026v4_RP = setParam(iMT1026v4_RP, 'lb', 'Ex_meoh', -1.5);  % Set methanol uptake %metadata.MaxUptakeRate(10))
    iMT1026v4_RP = setParam(iMT1026v4_RP, 'lb', 'Ex_o2', -3);

    metadataIdxs = find(strcmp(metadata.Molecule,RPs{RP}));
    for i = 1:length(metadataIdxs)
        % Determine max protein produced under MeOH
        
        % Handle different carbon sources
        if (~strcmp(metadata.Carbon{i}, 'MeOH'))
            % Find index of current carbon source
            carbIdx = ismember(carbons, metadata.Carbon{metadataIdxs(i)});
            
            % Set biomass production constraint and carbon uptake rate
            iMT1026v4_RP_i = setParam(iMT1026v4_RP, 'lb', {'Ex_biomass', carbonExRxns{carbIdx}}, ...
                                      [metadata.GrowthRate(i) , -metadata.MaxUptakeRate(i)]); %metadata.GrowthRate(i) 
             
            % Solve linear programming problem
            sol_bm_load = solveLP(iMT1026v4_RP_i, 1);
    
            % Determine max protein produced 
            intakeFlux = sol_bm_load.x(carbonExRxnsIdxs(carbIdx)); % mmol/gDCW/hr
            metadata.CofeedConsumptionRate(metadataIdxs(i)) = intakeFlux;
            metadata.timeToConsumption (metadataIdxs(i)) = (0.5 * 3 / 1000 * 1000) / (intakeFlux) / (metadata.Cellgrowth(metadataIdxs(i)) / 0.25 * 0.8);
            % (mol/L * 3 mL / 1000 mL/L * 1000 mmol/mol) / (mmol/gDCW/hr) / (OD
            % / 0.25 OD/gWCW * 0.8 gDCW/gWCW) == hr
        else
            % For methanol, only set biomass constraint
            iMT1026v4_RP_i = setParam(iMT1026v4_RP, 'lb', 'Ex_biomass', metadata.GrowthRate(i));
            sol_bm_load = solveLP(iMT1026v4_RP_i, 1);
            intakeFlux = -Inf;
        end
        intakeMeOHFlux = sol_bm_load.x(1985); % mmol/gDCW/hr
        metadata.MeOHConsumptionRate(metadataIdxs(i)) = intakeMeOHFlux;
        if (intakeFlux == 0)
            limitFlux = -intakeMeOHFlux;
        else
            limitFlux = min(-intakeMeOHFlux, -intakeFlux);
        end
        metadata.timeToConsumption (metadataIdxs(i)) = (1.5* 3 * 0.791 / 32.04 * 1000) / (limitFlux) / (metadata.Cellgrowth(metadataIdxs(i)) / 0.25 * 0.8);
        % (1.5 mL/mL * 3 mL * 0.791 g/mL / 32.04 g/mol * 1000 mmol/mol) 
        % / (mmol/gDCW/hr) / (OD / 0.25 OD/gWCW * 0.8 gDCW/gWCW) == hr
        % Calculate HSA production rate estimate (g/gDCW/day)
    
        metadata.ProdRateEst(metadataIdxs(i)) = sol_bm_load.f * RP_mass *  24;
    end
    
    
    models = [models, iMT1026v4_RP_original];
end
%% 4. Plot Correlation Between Predicted and Experimental Data
fprintf('Creating correlation plot with separate legends and molecule-specific fits...\n');
% Select data points to plot (excluding problematic or incomplete data)
plotVec = 1:height(metadata);
plotVec = find(~(strcmp(metadata.Carbon,'Fructose') | strcmp(metadata.Carbon,'Glycerol') | strcmp(metadata.Carbon,'MeOH')));%
% Extract predicted production rates
x_values = metadata.ProdRateEst(plotVec);
% Convert experimental values to matching units (g/gDCW/day)
y_values = metadata.specific_productivity(plotVec) * 0.25 / 0.8 / 1000 * 3 / 1000;
% Get carbon source names for color mapping
carbonNames = metadata.Carbon(plotVec);
% Get Molecule names for shape mapping
moleculeNames = metadata.Molecule(plotVec);
% Close any existing figures to prevent interference
close all;
% Create figure with appropriate size
figure('Position', [100, 100, 800, 600]);
mainAxes = gca; % Get handle to the main axes
hold(mainAxes, 'on'); % Ensure hold is on for the main plot
% --- Handle Carbon Legend ---
uniqueCarbons = unique(carbonNames);
numCarbons = length(uniqueCarbons);
colorMap = jet(numCarbons); % Create colormap
carbonHandles = gobjects(numCarbons, 1); % Pre-allocate handles for carbon legend
carbonLegendAxes = axes('Position', [0.95, 0.85, 0.01, 0.01], 'Visible', 'off'); % Tiny axes
hold(carbonLegendAxes, 'on');
for i = 1:numCarbons
    carbonHandles(i) = scatter(carbonLegendAxes, [], [], 100, colorMap(i,:), 'o', 'filled', 'MarkerEdgeColor', 'k');
end
carbonLegend = legend(carbonLegendAxes, carbonHandles, uniqueCarbons, 'FontSize', 9, 'Interpreter', 'none');
title(carbonLegend, 'Carbon Source');
set(carbonLegendAxes, 'Visible', 'off');
% --- Handle Molecule Legend and Molecule-Specific Fits ---
uniqueMolecules = unique(moleculeNames);
numMolecules = length(uniqueMolecules);
markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x', '+', '*'};
if numMolecules > length(markerShapes)
    markerShapes = repmat(markerShapes, 1, ceil(numMolecules/length(markerShapes)));
end
moleculeHandles = gobjects(numMolecules, 1); % Pre-allocate handles for molecule legend
moleculeLegendAxes = axes('Position', [0.95, 0.70, 0.01, 0.01], 'Visible', 'off'); % Tiny axes
hold(moleculeLegendAxes, 'on');
for i = 1:numMolecules
    % Get data specific to the current molecule
    moleculeIdxs = strcmp(moleculeNames, uniqueMolecules{i});
    x_mol = x_values(moleculeIdxs);
    y_mol = y_values(moleculeIdxs);
    % Calculate and plot line of best fit for this molecule
    if length(x_mol) > 1 % Need at least two points for a fit
        p_mol = polyfit(x_mol, y_mol, 1);
        x_fit_mol = linspace(min(x_mol), max(x_mol), 100);
        y_fit_mol = polyval(p_mol, x_fit_mol);
        plot(mainAxes, x_fit_mol, y_fit_mol, 'Color', [0.1 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1); % Gray dashed line
    end
    % Calculate Spearman correlation
    if length(x_mol) > 1
        rho = corr(x_mol, y_mol, 'Type', 'Spearman');
        legendLabel = sprintf('%s ($R_s$ = %.2f, slope = %.2f)', uniqueMolecules{i}, rho, p_mol(1));
    else
        legendLabel = uniqueMolecules{i};
    end
    moleculeHandles(i) = scatter(moleculeLegendAxes, [], [], 100, 'k', markerShapes{i}, 'filled');
    set(moleculeHandles(i), 'DisplayName', legendLabel); % Use DisplayName for legend
end
moleculeLegend = legend(moleculeLegendAxes, moleculeHandles, 'FontSize', 9, 'Interpreter', 'latex');
title(moleculeLegend, 'Molecule Name');
set(moleculeLegendAxes, 'Visible', 'off');
% --- Plot the actual data points on the main axes ---
for i = 1:length(x_values)
    carbonIdx = find(strcmp(uniqueCarbons, carbonNames{i}));
    moleculeIdx = find(strcmp(uniqueMolecules, moleculeNames{i}));
    currentShape = markerShapes{moleculeIdx};
    scatter(mainAxes, x_values(i), y_values(i), 100, colorMap(carbonIdx,:), currentShape, 'filled', 'MarkerEdgeColor', 'k');
end
% Add overall line of best fit
p_overall = polyfit(x_values, y_values, 1);
x_fit_overall = linspace(min(x_values), max(x_values), 100);
y_fit_overall = polyval(p_overall, x_fit_overall);
plot(mainAxes, x_fit_overall, y_fit_overall, 'k-', 'LineWidth', 2);
% Calculate overall R² and p-value
[r_overall, p_value_overall] = corr(x_values, y_values);
r_squared_overall = r_overall^2;
% Display overall R² and p-value on the main axes
text_string = sprintf('Overall R = %.4f\nOverall R² = %.4f\nOverall p-value = %.4f', r_overall, r_squared_overall, p_value_overall);
text(mainAxes, min(x_values) + 0.1*(max(x_values)-min(x_values)), ...
    max(y_values) - 0.1*(max(y_values)-min(y_values)), ...
    text_string, 'FontSize', 12, 'FontWeight', 'bold');
% Set axis limits with some padding for the main axes
x_range = max(x_values) - min(x_values);
y_range = max(y_values) - min(y_values);
xlim(mainAxes, [min(x_values) - 0.1*x_range, max(x_values) + 0.1*x_range]);
ylim(mainAxes, [min(y_values) - 0.1*y_range, max(y_values) + 0.1*y_range]);
% Make plot more professional for the main axes
set(mainAxes, 'Box', 'on', 'LineWidth', 1.5); % Add box and make it more visible
set(mainAxes, 'TickDir', 'out'); % Make tick marks point outward
set(mainAxes, 'FontSize', 10); % Increase font size of tick labels
% set(mainAxes, 'Xscale', 'log', 'Yscale', 'log'); % Change to log for easier viewing 
grid(mainAxes, 'on'); % Add grid lines
% Make x and y axis labels for the main axes
xlabel(mainAxes, 'Predicted Protein Flux ($\frac{g}{g_{DCW} \cdot day}$)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel(mainAxes, 'Experimental Specific Productivity ($\frac{g}{g_{DCW} \cdot day}$)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
% Optional: Add title to the main axes
% title(mainAxes, 'Correlation between Predicted and Experimental HSA Production', 'FontSize', 14, 'FontWeight', 'bold');
hold(mainAxes, 'off');

%% 5. FSEOF 
% Initialize an empty table to store the collated results
collatedTable = table(); % Initialize an empty table

for RP = 1:length(RPs)
    model = models{RP};
    model = setParam(model, 'lb', 'Ex_glyc', 0);
    model = setParam(model, 'lb', 'Ex_meoh', -1.5);
    model = setParam(model, 'lb', 'Ex_tre', -1.5);
    model = setParam(model, 'lb', 'Ex_o2', -0.5);
    targets = FSEOF_custom(model, 'Ex_biomass', sprintf('%s_toe',RPs{RP}), 100, 0.99, sprintf('FSEOF_results_%s.csv',RPs{RP}));
    slope_direction = double(abs(targets.solns(:,1))<abs(targets.solns(:,end)));
    slope_direction(slope_direction==0) = -(abs(targets.solns(slope_direction==0,1))>abs(targets.solns(slope_direction==0,end)));
    rxn_direction = double(targets.solns(:,1) > 0);
    rxn_direction(rxn_direction==0) = -(targets.solns(rxn_direction==0,1) < 0);
    FSEOFtargets = table(targets.slope, targets.logical, model.rev,...
        (1:length(targets.logical))', model.rxnNames, rxn_direction,...
        slope_direction, model.grRules,...
        'VariableNames',{'slope', 'target', 'reversibility', 'index', 'reaction name',...
        'rxndirection', 'slopeDirection','grRules'}, 'RowNames',model.rxns);
    lbModel = model;
    lbModel = setParam(lbModel, 'obj', 'Ex_biomass', 1);
    lbSol = solveLP(lbModel,1);
    ubModel = model;
    ubModel = setParam(ubModel, 'obj', sprintf('%s_toe', RPs{RP}), 1);
    ubSol = solveLP(ubModel,1);
    dsFSEOFtargets = FSEOFtargets(FSEOFtargets.slope~=0,:);
    dsFSEOFtargets = sortrows(dsFSEOFtargets,"slope","descend");
    dsFSEOFtargets = dsFSEOFtargets(dsFSEOFtargets.grRules ~= "",:);
    geneTargets = createGeneTargetsTable(dsFSEOFtargets);
    disp(['Number of rows in geneTargets before adding columns (RP = ' num2str(RP) '): ' num2str(height(geneTargets))]); % Debugging
    geneTargets = addPredictionColumns(geneTargets, model);
    geneTargets = sortrows(geneTargets,"predicted_change","descend");
    geneTargets = addAnnotationColumn(geneTargets, ppTable);
    geneTargets = addEssentialityAndTargetColumns(geneTargets, ppTable);

    geneTargets.Molecule = repmat(RPs{RP}, height(geneTargets), 1);
    geneTargets.Molecule = cellstr(repmat(RPs{RP}, height(geneTargets), 1));
    if RP == 1
        collatedTable = geneTargets;
    else
        collatedTable.Molecule = cellstr(collatedTable.Molecule); % Ensure consistent type
        geneTargets.Molecule = cellstr(geneTargets.Molecule);     % Ensure consistent type
        collatedTable = [collatedTable; geneTargets];
    end
    
end
collatedTable = collatedTable(collatedTable.Str_Eng_Target==1,:);
collatedTable = sortrows(collatedTable,"predicted_change","descend");
% Display the final collated table
disp(collatedTable);
%% 6. Target Summary
% Original code to create the main table
collatedTable = collatedTable(collatedTable.GQ_names~="",:);
% Get all unique genes and molecules
uniqueGenes = unique(collatedTable.gene);
uniqueMolecules = unique(collatedTable.Molecule);
% Initialize the new table structure
numUniqueGenes = numel(uniqueGenes);
commonColumns = {'gene', 'indices', 'category', 'FinalAnnotation', 'GQ_names', 'Essentiality', 'proteinType'};
% Create column names for molecule-specific columns
moleculeColumns = cell(1, length(uniqueMolecules));
for i = 1:length(uniqueMolecules)
    moleculeColumns{i} = [char(uniqueMolecules{i}), '_predicted_change'];
end
% Combine all column names
allColumnNames = [commonColumns, moleculeColumns];
% Initialize a cell array to hold all our data
allData = cell(numUniqueGenes, length(allColumnNames));

% Initialize a cell array to store information about genes with multiple categories
multiCatGenesInfo = {};

% Populate the data matrix
for i = 1:numUniqueGenes
    currentGene = uniqueGenes{i};
    geneRows = strcmp(collatedTable.gene, currentGene);
    % Get the first occurrence of this gene
    firstOccurrence = find(geneRows, 1, 'first');
    
    % Fill in the common columns
    for j = 1:length(commonColumns)
        allData{i, j} = collatedTable.(commonColumns{j})(firstOccurrence);
    end
    
    % Fill in molecule-specific predicted_change values
    for j = 1:length(uniqueMolecules)
        currentMolecule = uniqueMolecules{j};
        moleculeRows = strcmp(collatedTable.Molecule, currentMolecule);
        % Find the intersection of rows with this gene and this molecule
        targetRows = find(geneRows & moleculeRows);
        if ~isempty(targetRows)
            % If there's a match, store the predicted_change value
            allData{i, length(commonColumns) + j} = collatedTable.predicted_change(targetRows(1));
        else
            % If this gene wasn't tested with this molecule, store empty
            allData{i, length(commonColumns) + j} = [];
        end
    end
    
    % Check if the category is the same for all occurrences of this gene
    allCategories = collatedTable.category(geneRows);
    uniqueCategories = unique(allCategories);
    if length(uniqueCategories) > 1
        fprintf('Gene %s has multiple categories: ', currentGene);
        for k = 1:length(uniqueCategories)
            if k > 1
                fprintf(', ');
            end
            fprintf('%s', uniqueCategories{k});
        end
        fprintf('\n');
        
        % Collect information for the multi-category genes table
        geneIndices = find(geneRows);
        for idx = 1:length(geneIndices)
            rowIdx = geneIndices(idx);
            % Handle different possible data types for each column
            geneValue = collatedTable.gene(rowIdx);
            if iscell(geneValue)
                geneValue = geneValue{1};
            end
            
            moleculeValue = collatedTable.Molecule(rowIdx);
            if iscell(moleculeValue)
                moleculeValue = moleculeValue{1};
            end
            
            categoryValue = collatedTable.category(rowIdx);
            if iscell(categoryValue)
                categoryValue = categoryValue{1};
            end
            
            predictedChangeValue = collatedTable.predicted_change(rowIdx);
            % Don't use cell indexing for predicted_change, just access directly
            
            multiCatGenesInfo = [multiCatGenesInfo; {
                geneValue, ...
                moleculeValue, ...
                categoryValue, ...
                predictedChangeValue
            }];
        end
    end
end

% Create the main table in one go
outputTable = cell2table(allData, 'VariableNames', allColumnNames);

% Create the multi-category genes table
if ~isempty(multiCatGenesInfo)
    multiCatColumns = {'gene', 'Molecule', 'category', 'predicted_change'};
    multiCatTable = cell2table(multiCatGenesInfo, 'VariableNames', multiCatColumns);
    
    % Optional: Sort the table by gene name for better readability
    multiCatTable = sortrows(multiCatTable, 'gene');
    
    % Display summary information
    fprintf('\nSummary of genes with multiple categories:\n');
    multiCatGenes = unique(multiCatTable.gene);
    fprintf('Found %d genes with multiple categories\n', length(multiCatGenes));
    
    % Save the multi-category table to a separate variable
    % This makes it accessible in the MATLAB workspace
    assignin('base', 'multiCategoryGenesTable', multiCatTable);
end


% Check which genes have data for all molecules
hasAllMolecules = true(height(outputTable), 1);
for i = 1:height(outputTable)
    for j = 1:length(moleculeColumns)
        currentValue = outputTable.(moleculeColumns{j})(i);
        % Check if it's empty (handle both cell and direct values)
        if iscell(currentValue) && isempty(currentValue{1}) || ~iscell(currentValue) && isempty(currentValue)
            hasAllMolecules(i) = false;
            break;
        end
    end
end

% Final filtered table with only genes present in all molecules
finalOutputTable = outputTable(hasAllMolecules, :);

% 2. Add columns for average and variance of predicted_change
numRows = height(finalOutputTable);
avgPredictedChange = zeros(numRows, 1);
varPredictedChange = zeros(numRows, 1);

for i = 1:numRows
    % Extract predicted_change values for this gene
    predictedValues = zeros(1, length(moleculeColumns));
    
    for j = 1:length(moleculeColumns)
        % Rank order
        [sorted,indxs] = sort(cell2mat(finalOutputTable.(moleculeColumns{j})),'descend');
        [~,ranks] = sort(indxs);
        currentValue = ranks(i);
        if iscell(currentValue)
            predictedValues(j) = currentValue{1};
        else
            predictedValues(j) = currentValue;
        end
    end
    
    % Calculate average and variance
    avgPredictedChange(i) = mean(predictedValues);
    varPredictedChange(i) = var(predictedValues);
end

% Add these columns to the table
finalOutputTable.avg_predicted_change = avgPredictedChange;
finalOutputTable.var_predicted_change = varPredictedChange;

% 3. Sort the table in descending order based on average predicted_change
[~, sortIdx] = sort(avgPredictedChange, 'ascend');
finalOutputTable = finalOutputTable(sortIdx, :);

writetable(finalOutputTable,'commonTargets.csv')
%% Helper functions
% Create geneTargets table with unique genes and their corresponding data
function geneTargets = createGeneTargetsTable(dsFSEOFtargets)
    % Initialize a map to collect gene data
    geneMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Process each row in the original table
    for i = 1:height(dsFSEOFtargets)
        grRule = dsFSEOFtargets.grRules{i};
        currentIndex = dsFSEOFtargets.index(i);
        currentSlope = dsFSEOFtargets.slope(i);
        currentSlopeDirection = dsFSEOFtargets.slopeDirection(i);
        
        % Skip if grRule is empty
        if isempty(grRule)
            continue;
        end
        
        % Extract genes from the rule
        % Replace logical operators with spaces
        cleanRule = grRule;
        cleanRule = regexprep(cleanRule, '\(|\)', ' '); % Remove parentheses
        cleanRule = regexprep(cleanRule, ' and | or ', ' '); % Replace logical operators with spaces
        
        % Split into individual genes and remove any empty strings
        geneList = strsplit(strtrim(cleanRule));
        geneList = geneList(~cellfun('isempty', geneList));
        
        % For each gene in this rule
        for j = 1:length(geneList)
            gene = geneList{j};
            
            % Add or update the gene in the map
            if isKey(geneMap, gene)
                currentData = geneMap(gene);
                currentData.indices = [currentData.indices; currentIndex];
                currentData.slopes = [currentData.slopes; currentSlope];
                currentData.slopeDirections = [currentData.slopeDirections; currentSlopeDirection];
                geneMap(gene) = currentData;
            else
                geneMap(gene) = struct('indices', currentIndex, ...
                                      'slopes', currentSlope, ...
                                      'slopeDirections', currentSlopeDirection);
            end
        end
    end
    
    % Convert the map to a table
    uniqueGenes = keys(geneMap)';
    numGenes = length(uniqueGenes);
    
    geneIndicesList = cell(numGenes, 1);
    geneSlopesList = cell(numGenes, 1);
    geneSlopeDirectionsList = cell(numGenes, 1);
    
    for i = 1:numGenes
        gene = uniqueGenes{i};
        data = geneMap(gene);
        geneIndicesList{i} = data.indices;
        geneSlopesList{i} = data.slopes;
        geneSlopeDirectionsList{i} = data.slopeDirections;
    end
    
    % Create the geneTargets table
    geneTargets = table(uniqueGenes, geneIndicesList, geneSlopesList, geneSlopeDirectionsList, ...
                        'VariableNames', {'gene', 'indices', 'slopes', 'slopeDirections'});
end

% Function to add category and predicted_change columns to geneTargets
function geneTargets = addPredictionColumns(geneTargets, iMT1026v4_HSA)
    % Initialize new columns
numGenes = height(geneTargets);
category = cell(numGenes, 1);
predicted_change = zeros(numGenes, 1);
proteinType = cell(numGenes, 1); % New column for protein types

% Process each row
for i = 1:numGenes
    indices = geneTargets.indices{i};
    slopes = geneTargets.slopes{i};
    slopeDirections = geneTargets.slopeDirections{i};
    currentGene = geneTargets.gene{i}; % Get current gene name
    
    % Case 1: Only one index
    if length(indices) == 1
        predicted_change(i) = slopes;
        if slopeDirections > 0
            category{i} = 'Upregulation';
        else
            category{i} = 'Knockout';
        end
    % Case 2: Multiple indices
    else
        % Check if all slopeDirections are the same
        allPositive = all(slopeDirections > 0);
        allNegative = all(slopeDirections < 0);
        
        % Case 2.1 & 2.2: All slopeDirections are the same
        if allPositive
            category{i} = 'Upregulation';
            predicted_change(i) = sum(slopes);
        elseif allNegative
            category{i} = 'Knockout';
            predicted_change(i) = sum(slopes);
        % Case 2.3: Mixed slopeDirections
        else
            % Multiply slopeDirections with slopes and find sum
            weightedSum = sum(slopeDirections .* slopes);
            if weightedSum > 0
                category{i} = 'Mixed, Upregulation';
            else
                category{i} = 'Mixed, Knockout';
            end
            predicted_change(i) = abs(weightedSum);
        end
    end
    
    % Determine protein type for current gene
    proteinType{i} = determineProteinType(currentGene, iMT1026v4_HSA.grRules);
end

% Add new columns to the table
geneTargets.category = category;
geneTargets.predicted_change = predicted_change;
geneTargets.proteinType = proteinType; % Add the new proteinType column
end

% Function to determine the protein type for a gene
function typeResult = determineProteinType(geneName, grRules)
    % Initialize array to store roles found for this gene
    roles = {};
    
    % Find all entries in grRules where this gene appears
    for i = 1:length(grRules)
        currentRule = grRules{i};
        
        % Check if the gene appears in this rule
        if contains(currentRule, geneName)
            % Determine the role for this rule
            if contains(currentRule, ' or ')
                % There is an 'or' in the rule
                
                % Check if there's an 'and' associated with this gene
                % Need to analyze the rule structure more carefully
                
                % Simple check: if gene is within parentheses and there's an 'and'
                % within those parentheses, it's likely a complex subunit
                
                % Extract all subexpressions in parentheses
                subexprPatterns = extractBetween(currentRule, '(', ')');
                
                isComplexForParallel = false;
                isParallel = false;
                
                % Check each subexpression
                for j = 1:length(subexprPatterns)
                    if contains(subexprPatterns{j}, geneName) && contains(subexprPatterns{j}, ' and ')
                        isComplexForParallel = true;
                        break;
                    end
                end
                
                % If not found in a complex (parentheses with 'and')
                if ~isComplexForParallel
                    % Check if it's a direct part of an 'or' expression
                    isParallel = true;
                end
                
                if isComplexForParallel
                    roles{end+1} = 'Complex Subunit for Parallel Reaction';
                else
                    roles{end+1} = 'Parallel Reaction';
                end
            else
                % No 'or' in the rule
                if contains(currentRule, ' and ')
                    roles{end+1} = 'Complex Subunit';
                else
                    roles{end+1} = 'Single-Gene Enzyme';
                end
            end
        end
    end
    
    % Determine the overall role based on the hierarchy
    if isempty(roles)
        typeResult = 'Unknown'; % Gene not found in any rule
    else
        if any(strcmp(roles, 'Complex Subunit for Parallel Reaction'))
            typeResult = 'Complex Subunit for Parallel Reaction';
        elseif any(strcmp(roles, 'Complex Subunit')) && any(strcmp(roles, 'Parallel Reaction'))
            typeResult = 'Complex Subunit for Parallel Reaction';
        elseif any(strcmp(roles, 'Complex Subunit'))
            typeResult = 'Complex Subunit';
        elseif any(strcmp(roles, 'Parallel Reaction'))
            typeResult = 'Parallel Reaction';
        else
            typeResult = 'Single-Gene Enzyme';
        end
    end
end

% Function to add FinalAnnotation from geneConversion to geneTargets
function geneTargets = addAnnotationColumn(geneTargets, geneConversion)
    % Initialize new column for FinalAnnotation
    numGenes = height(geneTargets);
    finalAnnotations = cell(numGenes, 1);
    GQ_names = cell(numGenes,1);
    % Process each row in geneTargets
    for i = 1:numGenes
        currentGene = geneTargets.gene{i};
        
        % Find matching entry in geneConversion where PAS_names matches the gene
        matchIdx = find(strcmp(geneConversion.PAS_names, currentGene));
        
        % If a match is found, get the corresponding FinalAnnotation
        if ~isempty(matchIdx)
            % If multiple matches exist, take the first one (you may want to handle this differently)
            if length(matchIdx) > 1
                matchIdx = matchIdx(1);
            end
            finalAnnotations{i} = geneConversion.FinalAnnotation{matchIdx};
            GQ_names{i} = geneConversion.Gene{matchIdx};
        else
            % No match found
            finalAnnotations{i} = '';  % or 'N/A' or other placeholder as preferred
            GQ_names{i} = '';
        end

        
    end
    
    % Add the new column to geneTargets
    geneTargets.FinalAnnotation = finalAnnotations;
    geneTargets.GQ_names = GQ_names;
end

% Function to add Essentiality from essentialityTable and Str_Eng_Target column to geneTargets
function geneTargets = addEssentialityAndTargetColumns(geneTargets, essentialityTable)
    % Initialize new column for Essentiality
    numGenes = height(geneTargets);
    essentiality = NaN(numGenes, 1);  % Using NaN as default for numeric values
    
    % Process each row in geneTargets to add Essentiality
    for i = 1:numGenes
        % Get the GQ_name
        gqName = geneTargets.GQ_names{i};
        
        % Skip if empty
        if isempty(gqName)
            continue;
        end
        
        % Find matching entry in essentialityTable where geneID matches the GQ_name
        matchIdx = find(strcmp(essentialityTable.Gene, gqName));
        
        % If a match is found, get the corresponding x_logpval
        if ~isempty(matchIdx)
            % If multiple matches exist, take the first one
            if length(matchIdx) > 1
                matchIdx = matchIdx(1);
            end
            essentiality(i) = -log(essentialityTable.pval(matchIdx));
        end
    end
    
    % Add the Essentiality column to geneTargets
    geneTargets.Essentiality = essentiality;
    
    % Now add the Str_Eng_Target column based on criteria
    str_eng_target = zeros(numGenes, 1);  % Default to 0
    
    for i = 1:numGenes
        category = geneTargets.category{i};
        
        % Check conditions
        if strcmp(category, 'Upregulation') || strcmp(category, 'Mixed, Upregulation')
            str_eng_target(i) = 1;
        elseif (strcmp(category, 'Knockout') || strcmp(category, 'Mixed, Knockout')) && ...
               ~isnan(geneTargets.Essentiality(i)) && geneTargets.Essentiality(i) < 0.5
            str_eng_target(i) = 1;
        end
    end
    
    % Add the Str_Eng_Target column to geneTargets
    geneTargets.Str_Eng_Target = str_eng_target;
end

% This function generates perceptually distinct colors
% Credit: Timothy E. Holy, 2011
function colors = distinguishable_colors(n)
    if n <= 0
        colors = [];
        return;
    end
    
    % Generate the seed colors - basic RGB combinations
    baseColors = [
        1 0 0;  % red
        0 1 0;  % green
        0 0 1;  % blue
        1 1 0;  % yellow
        1 0 1;  % magenta
        0 1 1;  % cyan
        0.5 0.5 0.5;  % gray
        1 0.5 0;  % orange
        0.5 0 0.5;  % purple
        0 0.5 0.5;  % teal
        0.5 0.5 0;  % olive
        1 0.5 0.5;  % pink
        0.5 0 0;    % maroon
        0 0.5 0;    % forest green
        0 0 0.5;    % navy
    ];
    
    % If we need more colors than base colors, interpolate
    if n <= size(baseColors, 1)
        colors = baseColors(1:n, :);
    else
        % For additional colors, use HSV space with evenly spaced hues
        additionalNeeded = n - size(baseColors, 1);
        extraColors = zeros(additionalNeeded, 3);
        
        for i = 1:additionalNeeded
            h = (i-1)/additionalNeeded;
            s = 0.7 + 0.3*mod(i,2); % Alternating saturation
            v = 0.7 + 0.3*mod(floor((i-1)/2),2); % Alternating value
            extraColors(i,:) = hsv2rgb([h s v]);
        end
        
        colors = [baseColors; extraColors];
    end
end

function targets=FSEOF_custom(model,biomassRxn,targetRxn,iterations,coefficient,outputFile)
% FSEOF: implements the algorithm of Flux Scanning based on Enforced Objective Flux.
%
%   model           a model structure
%   biomassRxn      string with reaction ID of the biomass formation or
%                   growth reaction
%   targetRxn       string with reaction ID of target reaction
%   iterations      double indicating number of iterations (opt, default 10)
%   coefficient     double indicating ratio of optimal target reaction
%                   flux, must be less than 1 (opt, default 0.9)
%   outputFile      string with output filename (opt, default prints to
%                   command window)
%
%   targets         structure with target identifying information for each reaction
%       logical     logical array indicating FSEOF identified target reaction
%       slope       double array with FSEOF calculated slopes for each reaction

%OUTPUTS
%   This function writes an tab-delimited file or prints to command window. If an output
%   has been specified (targets), it will also generate a structure indicating for
%   each reaction whether it is identified by FSEOF as a target and the slope of the
%   reaction when switching from biomass formation to product formation.
%
%   This function has been modified from the version published in the RAVEN
%   toolbox V2.0.3. Modified version also outputs the raw linear
%   programming solutions at each iteration of the FSEOF algorithm.
%
%   Usage: targets=FSEOF(model,biomassRxn,targetRxn,iterations,coefficient,outputFile)

biomassRxn=char(biomassRxn);
targetRxn=char(targetRxn);

if nargin<4
    iterations=10;
    coefficient=0.9;
end

if nargin <5
    coefficient=0.9;
end

if nargin == 6
    output=1;
else
    output=0;
end

%Find out the maximum theoretical yield of target reaction
model=setParam(model,'obj',targetRxn,1);
sol=solveLP(model,1);
targetMax=abs(sol.f*coefficient);   % 90 percent of the theoretical yield

model=setParam(model,'obj',biomassRxn,1);

fseof.results=zeros(length(model.rxns),iterations);
fseof.target=zeros(length(model.rxns),1);
rxnDirection=zeros(length(model.rxns),1);

%Enforce objective flux iteratively
for i=1:iterations
    n=i*targetMax/iterations;
    model=setParam(model,'lb',targetRxn,n);
    
    sol=solveLP(model,1);
    
    fseof.results(:,i)=sol.x;
    
    %Loop through all fluxes and identify the ones that increased upon the
    %enforced objective flux
    for j=1:length(fseof.results)
        if fseof.results(j,1) > 0   %Check the positive fluxes
            
            if i == 1   %The initial round
                rxnDirection(j,1)=1;
                fseof.target(j,1)=1;
            else
                
                if (fseof.results(j,i) > fseof.results(j,i-1)) & fseof.target(j,1)
                    fseof.target(j,1)=1;
                else
                    fseof.target(j,1)=0;
                end
            end
            
        elseif fseof.results(j,1) < 0 %Check the negative fluxes
            
            if i == 1   %The initial round
                rxnDirection(j,1)=-1;
                fseof.target(j,1)=1;
            else
                if (fseof.results(j,i) < fseof.results(j,i-1)) & fseof.target(j,1)
                    fseof.target(j,1)=1;
                else
                    fseof.target(j,1)=0;
                end
            end
            
        end
        
    end
end

%Generating output
formatSpec='%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
if output == 1    %Output to a file
    outputFile=char(outputFile);
    fid=fopen(outputFile,'w');
    fprintf(fid,formatSpec,'Slope','rowID','Enzyme ID','Enzyme Name','Subsystems','Direction','Gr Rule');
else              %Output to screen
    fprintf(formatSpec,'Slope','rowID','Enzyme ID','Enzyme Name','Subsystems','Direction','Gr Rule');
end

for num=1:length(fseof.target)
    if fseof.target(num,1) == 1
        A0=num2str(abs(fseof.results(num,iterations)-fseof.results(num,1))/abs(targetMax-targetMax/iterations)); %Slope calculation
        A1=num2str(num);                                  %row ID
        A2=char(model.rxns(num));                         %enzyme ID
        A3=char(model.rxnNames(num));                     %enzyme Name
        if isfield(model,'subSystems') && ~isempty(model.subSystems{num});
            A4=char(strjoin(model.subSystems{num,1},';'));                   %Subsystems
        else
            A4='';
        end
        A5=num2str(model.rev(num)*rxnDirection(num,1));   %reaction Dirction
        A6=char(model.grRules(num));                      %Gr Rule
        if output == 1    %Output to a file
            fprintf(fid,formatSpec,A0,A1,A2,A3,A4,A5,A6);
        else              %Output screen
            fprintf(formatSpec,A0,A1,A2,A3,A4,A5,A6);
        end
    end
end

if output == 1    %Output to a file
    fclose(fid);
end

if nargout == 1
    targets.logical=logical(fseof.target);
    targets.slope=abs(fseof.results(:,iterations)-fseof.results(:,1))/abs(targetMax-targetMax/iterations); %Slope calculation
    targets.solns=fseof.results;
end
end
