2%% HSA Production Model Analysis
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

figure('Name', 'FSEOF Summary Visualization', 'Position', [100, 100, 1400, 1000]);
t = tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');
%% 1. Create Base Model and Add HSA Production
% Add Necessary Metabolites and Reactions
% Excel file consisting of two sheets - one with new metabolites and one
% with the corresponding necessary reactions. 
% Define the file name

% Import base iMT1026v3 model
fprintf('Importing base metabolic model...\n');
iMT1026v3 = importModel(fullfile(repoRoot,'models','MODEL1612130000.xml'), false); % Obtained from https://www.ebi.ac.uk/biomodels/MODEL1612130000#Files
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

numMetRxns = length(iMT1026v4.rxns);


fprintf('Generating HSA synthesis reactions...\n');
ppTable = readtable(fullfile(repoRoot,'data','ppPSIMv3_RP.csv'));

HSAidx = 5194; 
% With PTMs
proteinInfo.aaSequence = ppTable.Sequence{HSAidx};
proteinInfo.dnaSequence = ppTable.dnaSequence{HSAidx};
proteinInfo.throughER = ppTable.ThroughER(HSAidx);
proteinInfo.numNGlyc = ppTable.NumNGlyc(HSAidx);
proteinInfo.numDisulfide = ppTable.NumDisulfides(HSAidx);
proteinInfo.transmembrane = ppTable.Transmembrane(HSAidx);
proteinInfo.numOGlyc = ppTable.NumOGlyc(HSAidx);
proteinInfo.GPIAddition = ppTable.GPIsite(HSAidx);
proteinInfo.finalLocation = strsplit(ppTable.Localization{HSAidx},',');

% Without PTMs
% proteinInfo.aaSequence = ppTable.Sequence{HSAidx};
% proteinInfo.dnaSequence = ppTable.dnaSequence{HSAidx};
% proteinInfo.throughER = 0;
% proteinInfo.numNGlyc = 0;
% proteinInfo.numDisulfide = 0;
% proteinInfo.transmembrane = 0;
% proteinInfo.numOGlyc = 0;
% proteinInfo.GPIAddition = 0;
% proteinInfo.finalLocation = strsplit(ppTable.Localization{HSAidx},',');

iMT1026v4_HSA = generateSynthesisReactions4(iMT1026v4,ppTable.Gene{HSAidx},proteinInfo);
iMT1026v4_HSA_original = iMT1026v4_HSA;  % Store original model for reuse
% Check model structure and production capability

produced = canProduce(iMT1026v4_HSA, 'HSA_mature_e');
fprintf('HSA production possible: %d\n', produced);

%% 3. Validate Model with Different Carbon Sources
% Define carbon sources and exchange reactions
fprintf('Validating model with experimental data...\n');
carbonExRxns = {'Ex_sbt_D', 'Ex_arab_D', 'Ex_fru', 'Ex_gal', 'Ex_sucr', 'Ex_xyl_D', 'Ex_rmn', 'Ex_cit', 'Ex_glyc'};
carbons = {'Sorbitol', 'D-Arabinose', 'Fructose', 'Galactose', 'Sucrose', 'D-Xylose', 'Rhamnose', 'Citrate', 'Glycerol'};
[~, carbonExRxnsIdxs] = ismember(carbonExRxns, iMT1026v4_HSA_original.rxns);

% Load experimental data
metadata = readtable('HSA_carbonScreen_metadata.xlsx');

% Calculate growth rate from experimental data (cellgrowth-5)/24/7 in 1/hr
metadata.GrowthRate = (metadata.Cellgrowth - 5) / 24 / 7;

% Reset model to original state and prepare for simulations
iMT1026v4_HSA = iMT1026v4_HSA_original;
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'obj', 'HSA_toe', 1);  % Set HSA synthesis as objective
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_glyc', 0);  % Zero glycerol uptake by default
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_tre', -0.1);
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_meoh', -1.5);  % Set methanol uptake %metadata.MaxUptakeRate(10))
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_o2', -3);

% Initialize production rate estimates
metadata.ProdRateEst = zeros(height(metadata), 1);
metadata.timeToConsumption = zeros(height(metadata),1);
metadata.CofeedConsumptionRate = zeros(height(metadata), 1);
metadata.MeOHConsumptionRate = zeros(height(metadata), 1);
% Calculate per gram molarity for uptake rate scaling
metadata.perGramMolarity = -.5 .* 1000 ./ (metadata.GrowthRate ./ 0.25 .* 0.8) ./ 24;

% Run simulations for each carbon source
for i = 1:height(metadata)
    % Determine max protein produced under MeOH
    
    % Handle different carbon sources
    if (~strcmp(metadata.Carbon{i}, 'MeOH'))
        % Find index of current carbon source
        carbIdx = ismember(carbons, metadata.Carbon{i});
        %iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_tre', -0.1);
        % Set biomass production constraint and carbon uptake rate
        iMT1026v4_HSA_i = setParam(iMT1026v4_HSA, 'lb', {'Ex_biomass', carbonExRxns{carbIdx}}, ...
                                  [metadata.GrowthRate(i), -metadata.MaxUptakeRate(i)]);
        
        % Solve linear programming problem
        sol_bm_load = solveLP(iMT1026v4_HSA_i, 1);

        % Determine max protein produced 
        intakeFlux = sol_bm_load.x(carbonExRxnsIdxs(carbIdx)); % mmol/gDCW/hr
        metadata.CofeedConsumptionRate(i) = intakeFlux;
        metadata.timeToConsumption (i) = (0.5 * 3 / 1000 * 1000) / (intakeFlux) / (metadata.Cellgrowth(i) / 0.25 * 0.8);
        % (mol/L * 3 mL / 1000 mL/L * 1000 mmol/mol) / (mmol/gDCW/hr) / (OD
        % / 0.25 OD/gWCW * 0.8 gDCW/gWCW) == hr
    else
        % For methanol, only set biomass constraint
        iMT1026v4_HSA_i = setParam(iMT1026v4_HSA, 'lb', 'Ex_biomass', metadata.GrowthRate(i));
        sol_bm_load = solveLP(iMT1026v4_HSA_i, 1);
        intakeFlux = -Inf;
    end
    intakeMeOHFlux = sol_bm_load.x(1985); % mmol/gDCW/hr
    metadata.MeOHConsumptionRate(i) = intakeMeOHFlux;
    if (intakeFlux == 0)
        limitFlux = -intakeMeOHFlux;
    else
        limitFlux = min(-intakeMeOHFlux, -intakeFlux);
    end
    metadata.timeToConsumption (i) = (1.5* 3 * 0.791 / 32.04 * 1000) / (limitFlux) / (metadata.Cellgrowth(i) / 0.25 * 0.8);
    % (1.5 mL/mL * 3 mL * 0.791 g/mL / 32.04 g/mol * 1000 mmol/mol) 
    % / (mmol/gDCW/hr) / (OD / 0.25 OD/gWCW * 0.8 gDCW/gWCW) == hr
    % Calculate HSA production rate estimate (g/gDCW/day)

    metadata.ProdRateEst(i) = sol_bm_load.f / ppTable.proteinSize(HSAidx)  * 24;
end

%% 4. Plot Correlation Between Predicted and Experimental Data
fprintf('Creating correlation plot...\n');

% Select data points to plot (excluding problematic or incomplete data)
plotVec = [4:9, 13:24, 28:30];
plotVec = [4:24, 28:30];
% plotVec = [1:30]

% Extract predicted production rates
x_values = metadata.ProdRateEst(plotVec);

% Convert experimental values to matching units (g/gDCW/day)
% mg/L/OD/day * 0.25 OD/gWCW / 0.8 gDCW/gWCW / 1000mg/g * 3mL * 1L/1000mL
y_values = metadata.specific_productivity(plotVec) * 0.25 / 0.8 / 1000 * 3 / 1000;

% Get carbon source names for color mapping
carbonNames = metadata.Carbon(plotVec);

% Create figure with appropriate size
nexttile(t, 1, [2, 2]);

% Get unique carbon sources for coloring
uniqueCarbons = unique(carbonNames);
numCarbons = length(uniqueCarbons);
% Define carbon source color mapping from R
csource_colors = containers.Map(...
    uniqueCarbons, ...
    {'#FFFFB3', '#80B1D3', '#FCCDE5', '#8DD3C7', ...
     '#BEBADA', '#FB8072', '#FDB462', '#B3DE69'});

% Convert hex colors to RGB
hex2rgb = @(hex) sscanf(hex(2:end),'%2x%2x%2x',[1 3])/255;

% Generate colorMap in the order of uniqueCarbons
colorMap = zeros(length(uniqueCarbons), 3);

for i = 1:length(uniqueCarbons)
    colorMap(i,:) = hex2rgb(csource_colors(uniqueCarbons{i}));
end


% Plot each point with color based on carbon source
hold on;
for i = 1:length(x_values)
    carbonIdx = find(strcmp(uniqueCarbons, carbonNames{i}));
    scatter(x_values(i), y_values(i), 100, colorMap(carbonIdx,:), 'filled', 'MarkerEdgeColor', 'k');
end

% Add line of best fit
p = polyfitZero(x_values, y_values, 1);
x_fit = linspace(min(x_values), max(x_values), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);

% Calculate R² and p-value
[r, p_value] = corr(x_values, y_values, 'Type', 'Spearman');

% Display R² and p-value on the plot
text_string = sprintf('slope = %.2f\nR = %.2f\np-value = %.4f', p(1), r, p_value);
text(min(x_values) + 0.1*(max(x_values)-min(x_values)), ...
    max(y_values) - 0.05*(max(y_values)-min(y_values)), ...
    text_string, 'FontSize', 12, 'FontWeight', 'bold');

% Create a legend with dummy points
h = zeros(length(uniqueCarbons), 1);
for i = 1:length(uniqueCarbons)
    h(i) = scatter(-1, -1, 100, colorMap(i,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', uniqueCarbons{i});
end
legend(h, 'Location', 'northwest', 'FontSize', 10);

% Set axis limits with some padding
x_range = max(x_values) - min(x_values);
y_range = max(y_values) - min(y_values);
xlim([min(x_values) - 0.1*x_range, max(x_values) + 0.1*x_range]);
ylim([min(y_values) - 0.1*y_range, max(y_values) + 0.1*y_range]);

% Make plot more professional
set(gca, 'Box', 'on', 'LineWidth', 1.5);  % Add box and make it more visible
set(gca, 'TickDir', 'out');  % Make tick marks point outward
set(gca, 'FontSize', 10);  % Increase font size of tick labels
grid on;  % Add grid lines

% Make x and y axis labels
xlabel('Predicted HSA Flux ($\frac{g}{g_{DCW} \cdot day}$)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Experimental HSA specific productivity ($\frac{g}{g_{DCW} \cdot day}$)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');

hold off;

%% 5. Run FSEOF on HSA model
fprintf('Running FSEOF... \n')
iMT1026v4_HSA = iMT1026v4_HSA_original;
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_glyc', 0);  % Zero glycerol uptake by default
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_meoh', -1.5);  % Set methanol uptake %metadata.MaxUptakeRate(10))
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_tre', -1.5);
iMT1026v4_HSA = setParam(iMT1026v4_HSA, 'lb', 'Ex_o2', -0.5);
targets = FSEOF_custom(iMT1026v4_HSA, 'Ex_biomass', 'HSA_toe', 500, 0.99,'FSEOF_results.csv');

slope_direction = double(abs(targets.solns(:,1))<abs(targets.solns(:,end)));
slope_direction(slope_direction==0) = -(abs(targets.solns(slope_direction==0,1))>abs(targets.solns(slope_direction==0,end)));
%Slope_direction is +1 if the absolute flux of the reaction is increasing
%with increasing target; - 1 if the absolute flux of the reaction is
%decreasing with increasing target.

rxn_direction = double(targets.solns(:,1) > 0);
rxn_direction(rxn_direction==0) = -(targets.solns(rxn_direction==0,1) < 0);
%rxn_direction is +1 if the initial flux is positive, -1 if the initial
%flux is negative 

FSEOFtargets = table(targets.slope, targets.logical, iMT1026v4_HSA_original.rev,...
    (1:length(targets.logical))', iMT1026v4_HSA_original.rxnNames, rxn_direction,...
    slope_direction, iMT1026v4_HSA_original.grRules,...
    'VariableNames',{'slope', 'target', 'reversibility', 'index', 'reaction name',...
    'rxndirection', 'slopeDirection','grRules'}, 'RowNames',iMT1026v4_HSA_original.rxns);

lbiMT1026v4_HSA = iMT1026v4_HSA;
lbiMT1026v4_HSA = setParam(lbiMT1026v4_HSA, 'obj', 'Ex_biomass', 1);
lbSol = solveLP(lbiMT1026v4_HSA,1);
ubiMT1026v4_HSA = iMT1026v4_HSA;
ubiMT1026v4_HSA = setParam(ubiMT1026v4_HSA, 'obj', 'HSA_toe', 1);
ubSol = solveLP(ubiMT1026v4_HSA,1);

%% 6. Analyze FSEOF, Identify Genes
dsFSEOFtargets = FSEOFtargets(FSEOFtargets.slope~=0,:);
dsFSEOFtargets = sortrows(dsFSEOFtargets,"slope","descend");
dsFSEOFtargets = dsFSEOFtargets(dsFSEOFtargets.grRules ~= "",:);

geneTargets = createGeneTargetsTable(dsFSEOFtargets);
geneTargets = addPredictionColumns(geneTargets, iMT1026v4_HSA);
geneTargets = sortrows(geneTargets,"predicted_change","descend");
geneTargets = addAnnotationColumn(geneTargets, ppTable);
geneTargets = addEssentialityAndTargetColumns(geneTargets, ppTable);
geneTargets = sortrows(geneTargets,"Str_Eng_Target","descend");

writetable(geneTargets(geneTargets.Str_Eng_Target==1,[1,5:9]),'geneTargets.csv')
xvals = [lbSol.x(end), targets.solns(end,:), ubSol.x(end)];
yvals = ([lbSol.x, targets.solns, ubSol.x]); % - lbSol.x
slopes = targets.slope;
% Find the top 5% of lines based on slope values
[~, sortedIndices] = sort(slopes, 'descend');
numLinesToHighlight = ceil(0.01 * length(slopes)); % Calculate 5% of total lines
topIndices = sortedIndices(1:numLinesToHighlight);

% % Create a new figure with appropriate size
% hFig = figure('Position', [100, 100, 1000, 600]);
% hAx = axes('Parent', hFig);
% hold(hAx, 'on');

% % Plot handles for the legend
% plotHandles = zeros(numLinesToHighlight, 1);
% 
% % First plot all lines in transparent grey
% for i = 1:size(yvals, 1)
%     if ~ismember(i, topIndices) % Only plot the lines that are not in top 5%
%         plot(hAx, xvals, yvals(i,:), 'Color', [0.7, 0.7, 0.7, 0.4]); % Light grey with alpha=0.1
%     end
% end

% Generate distinct colors for the top lines
colors = distinguishable_colors(numLinesToHighlight);

% Create a cell array to store legend entries
legendEntries = cell(numLinesToHighlight, 1);

% % Plot the top lines with distinct colors and add to legend
% for i = 1:numLinesToHighlight
%     lineIndex = topIndices(i);
%     plotHandles(i) = plot(hAx, xvals, yvals(lineIndex,:), 'Marker', '.', 'LineStyle', 'none', 'LineWidth', 2, 'Color', colors(i,:));
% 
%     % Get the reaction name for this line and store for legend
%     rxnName = iMT1026v4_HSA_original.rxnNames{lineIndex};
%     legendEntries{i} = rxnName;
% end

% % Add labels and title
% xlabel(hAx, 'Flux HSA (mmol/g_{DCW}/hr)');
% ylabel(hAx, '\Delta Flux Reaction');
% title(hAx, 'Metabolic Reaction Profiles with Top 1% Slopes Highlighted');

% % Add legend - using the plot handles for proper legend creation
% hLegend = legend(hAx, plotHandles, legendEntries, 'Location', 'eastoutside', 'FontSize', 8);
% 
% % Add grid
% grid(hAx, 'on');
% % Set box property safely
% set(hAx, 'Box', 'on', 'FontSize', 12);
% 
% % Optional: Add a note about the highlighted lines
% annotation(hFig, 'textbox', [0.15, 0.01, 0.7, 0.05], 'String', ...
%     ['The highlighted lines represent the top 1% reactions based on slope (' num2str(numLinesToHighlight) ' out of ' num2str(length(slopes)) ' total reactions)'], ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% 
% hold off;

%% === Barplot of Top 1% Reactions by Slope, Colored by Subsystem ===
fprintf('Plotting barplot of top reactions by FSEOF slope...\n')

% Extract top reaction names, slopes, and subsystems
topRxnNames = iMT1026v4_HSA_original.rxnNames(topIndices);
topSlopes = slopes(topIndices);
topSubsystems = iMT1026v4_HSA.subSystems(topIndices);

% Ensure subsystem names are consistent for coloring
% Ensure subsystems are character vectors
% If they're nested cells, extract the first string or use 'Unknown'
topSubsystemsFlat = cell(size(topSubsystems));
for i = 1:numel(topSubsystems)
    if iscell(topSubsystems{i}) && ~isempty(topSubsystems{i})
        topSubsystemsFlat{i} = topSubsystems{i}{1};
    elseif ischar(topSubsystems{i})
        topSubsystemsFlat{i} = topSubsystems{i};
    else
        topSubsystemsFlat{i} = 'Unknown';
    end
end

[uniqueSubsystems, ~, subsystemIdx] = unique(topSubsystemsFlat);


% Generate distinguishable colors for subsystems
subsystemColors = distinguishable_colors(numel(uniqueSubsystems));

% Create color map for bar plot
barColors = subsystemColors(subsystemIdx, :);

% Create figure
figure

% Create bar plot
b = bar(topSlopes, 'FaceColor', 'flat');
b.CData = barColors;

% Add labels
xticks(1:numel(topRxnNames));
xticklabels(topRxnNames);
xtickangle(45);
ylabel('Slope');
title('Top 1% Reactions by FSEOF Slope Colored by Subsystem');
grid on;

% Create legend for subsystems
% Create dummy plot handles for legend (ensure color matches subsystem)
hold on;
legendHandles = gobjects(numel(uniqueSubsystems), 1);
for i = 1:numel(uniqueSubsystems)
    legendHandles(i) = plot(nan, nan, 's', 'MarkerFaceColor', subsystemColors(i,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8);
end
legend(legendHandles, uniqueSubsystems, 'Location', 'northeast', 'FontSize', 8);
hold off;


%% === Histogram of All Slopes ===
% Compute histogram data first
[counts, edges] = histcounts(slopes, 50);
barCenters = edges(1:end-1) + diff(edges)/2;
figure
% Create nested axes for broken y-axis within current tile
parentPos = get(gca, 'Position');  % Save tile position
delete(gca);  % Clear tile

% Custom axes positions
breakHeight = 25;   % Maximum height for bottom axis
yTopMax = max(counts);
bottomFrac = 0.8;
topFrac = 1 - bottomFrac; % Fraction for top axis

gap = 0.02;
totalHeight = parentPos(4);
axBot = axes('Position', [parentPos(1), parentPos(2), parentPos(3), totalHeight * bottomFrac]);
axTop = axes('Position', [parentPos(1), parentPos(2) + totalHeight * bottomFrac + gap, ...
                          parentPos(3), totalHeight * (1 - bottomFrac) - gap]);

% Draw histogram on both axes
bar(axBot, barCenters, counts, 1, 'FaceColor', [0.3 0.3 0.8], 'EdgeColor', 'k');
ylim(axBot, [0, breakHeight]);
xlabel(axBot, 'FSEOF Slope');
ylabel(axBot, 'Count');

bar(axTop, barCenters, counts, 1, 'FaceColor', [0.3 0.3 0.8], 'EdgeColor', 'k');
ylim(axTop, [breakHeight, yTopMax]);
set(axTop, 'XTick', []);
ylabel(axTop, 'Count');

linkaxes([axTop, axBot], 'x');

% Draw dashed break line
annotation('line', ...
    [parentPos(1), parentPos(1) + parentPos(3)], ...
    [parentPos(2) + totalHeight * bottomFrac + gap, ...
     parentPos(2) + totalHeight * bottomFrac + gap], ...
    'Color', 'k', 'LineStyle', '--');


%% 7. Visualize Gene Targets

fprintf('Creating gene target bar plots... \n')

% Add a 'plotting_category' column to geneTargets for easier filtering
geneTargets.plotting_category = repmat("", size(geneTargets, 1), 1);
geneTargets.plotting_category(contains(geneTargets.category,'Upregulation')) = "Upregulation";
geneTargets.plotting_category(contains(geneTargets.category,'Knockout')) = "Knockout";
% Handle mixed categories: if 'Mixed, Upregulation', treat as Upregulation, etc.
geneTargets.plotting_category(contains(geneTargets.plotting_category, 'Mixed, Upregulation')) = "Upregulation";
geneTargets.plotting_category(contains(geneTargets.plotting_category, 'Mixed, Knockout')) = "Knockout";

% Create the 'max_slope' column from the 'slopes' column.
% This assumes 'slopes' is a cell array where each cell contains a numeric vector of slopes.
% For each gene, we take the maximum slope from its associated vector.
geneTargets.max_slope = zeros(size(geneTargets, 1), 1);
for i = 1:size(geneTargets, 1)
    currentSlopes = geneTargets.slopes{i};
    if ~isempty(currentSlopes)
        geneTargets.max_slope(i) = max(currentSlopes);
    else
        geneTargets.max_slope(i) = 0; % Assign 0 if slopes are empty for a gene, adjust if NaN is preferred
    end
end


% Get unique protein types for coloring.
% Using 'proteinType' as per your correction.
uniqueProteinTypes = unique(geneTargets.proteinType);
numProteinTypes = length(uniqueProteinTypes);

% Define the specific pastel colors
fixedColors = [
    127 201 127; % #7fc97f
    190 174 212; % #beaed4
    253 192 134; % #fdc086
    255 255 153  % #ffff99
] / 255; % Normalize to 0-1 range

% Assign colors from fixedColors to protein types, cycling if needed
proteinTypeColors = zeros(numProteinTypes, 3);
for i = 1:numProteinTypes
    proteinTypeColors(i, :) = fixedColors(mod(i-1, size(fixedColors, 1)) + 1, :);
end

% Create a map for protein type to color for consistent coloring
colorMap = containers.Map(uniqueProteinTypes, num2cell(proteinTypeColors, 2));


%% Plot for Upregulation Targets
upregulationTargets = geneTargets(geneTargets.plotting_category == "Upregulation", :);
% Filter for max_slopes > 10
min_max_slope = 10;
upregulationTargets = upregulationTargets(upregulationTargets.max_slope > min_max_slope, :);
upregulationTargets = sortrows(upregulationTargets, "max_slope", "descend");

if ~isempty(upregulationTargets)
    % hFigUp = figure('Position', [100, 100, 1200, 700]);
    % hAxUp = axes('Parent', hFigUp);
    % hold(hAxUp, 'on');
    figure
    hold on
    % Create colors array for bars
    barColorsUp = zeros(size(upregulationTargets, 1), 3);
    for i = 1:size(upregulationTargets, 1)
        % Use 'proteinType' for coloring
        barColorsUp(i, :) = colorMap(upregulationTargets.proteinType{i});
    end

    % Escape underscores in FinalAnnotation for plotting
    upregulationLabels = cellfun(@(x) strrep(x, '_', '\_'), upregulationTargets.FinalAnnotation, 'UniformOutput', false);

    % Plot horizontal bars. Note: For descending order from top to bottom,
    % we plot the data in reverse order after sorting descending.
    % The CData property needs to be assigned the entire flipped color array.
    bUp = barh(flipud(upregulationLabels), flipud(upregulationTargets.max_slope)); % hAxUp, 
    bUp.FaceColor = 'flat'; % Set face color mode to flat
    bUp.CData = flipud(barColorsUp); % Assign the flipped colors directly

    % Add labels and title (swapped for horizontal bars)
    ylabel('Gene Name');
    xlabel('Maximum Slope');
    title('Upregulation Gene Targets by Maximum Slope (Slope > 10)');
    % Removed YTickLabelInterpreter as it's not supported by your MATLAB version


    % Create a custom legend for protein types
    legendHandlesUp = gobjects(numProteinTypes, 1);
    legendEntriesUp = cell(numProteinTypes, 1);
    for i = 1:numProteinTypes
        legendHandlesUp(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', proteinTypeColors(i,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
        legendEntriesUp{i} = uniqueProteinTypes{i};
    end
    % Position legend in the bottom right within the plotting range
    legend(legendHandlesUp, legendEntriesUp, 'Location', 'SouthEast');

    % grid(hAxUp, 'on');
    % set(hAxUp, 'Box', 'on', 'FontSize', 10);
    % hold(hAxUp, 'off');
    hold off
else
    fprintf('No Upregulation targets with slope > 10 found to plot.\n');
end


%% Plot for Knockout Targets
knockoutTargets = geneTargets(geneTargets.plotting_category == "Knockout", :);
% Filter for max_slopes > 10
knockoutTargets = knockoutTargets(knockoutTargets.max_slope > min_max_slope, :);
knockoutTargets = sortrows(knockoutTargets, "max_slope", "descend");

if ~isempty(knockoutTargets)
    figure
    hold on
    % hFigKO = figure('Position', [100, 100, 1200, 700]);
    % hAxKO = axes('Parent', hFigKO);
    % hold(hAxKO, 'on');

    % Create colors array for bars
    barColorsKO = zeros(size(knockoutTargets, 1), 3);
    for i = 1:size(knockoutTargets, 1)
        % Use 'proteinType' for coloring
        barColorsKO(i, :) = colorMap(knockoutTargets.proteinType{i});
    end

    % Escape underscores in FinalAnnotation for plotting
    knockoutLabels = cellfun(@(x) strrep(x, '_', '\_'), knockoutTargets.FinalAnnotation, 'UniformOutput', false);

    % Plot horizontal bars. Note: For descending order from top to bottom,
    % we plot the data in reverse order after sorting descending.
    % The CData property needs to be assigned the entire flipped color array.
    bKO = barh(flipud(knockoutLabels), flipud(knockoutTargets.max_slope));
    bKO.FaceColor = 'flat'; % Set face color mode to flat
    bKO.CData = flipud(barColorsKO); % Assign the flipped colors directly

    % Add labels and title (swapped for horizontal bars)
    ylabel('Gene Name');
    xlabel('Maximum Slope');
    title('Knockout Gene Targets by Maximum Slope (Slope > 10)');
    % Removed YTickLabelInterpreter as it's not supported by your MATLAB version

    % Create a custom legend for protein types
    legendHandlesKO = gobjects(numProteinTypes, 1);
    legendEntriesKO = cell(numProteinTypes, 1);
    for i = 1:numProteinTypes
        legendHandlesKO(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', proteinTypeColors(i,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
        legendEntriesKO{i} = uniqueProteinTypes{i};
    end
    % Position legend in the bottom right within the plotting range
    legend(legendHandlesKO, legendEntriesKO, 'Location', 'SouthEast');

    % grid(hAxKO, 'on');
    % set(hAxKO, 'Box', 'on', 'FontSize', 10);
    % hold(hAxKO, 'off');
    hold off
else
    fprintf('No Knockout targets with slope > 10 found to plot.\n');
end

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


function [p,S,mu] = polyfitZero(x,y,degree)
% POLYFITZERO Fit polynomial to data, forcing y-intercept to zero.
%   P = POLYFITZERO(X,Y,N) is similar POLYFIT(X,Y,N) except that the
%   y-intercept is forced to zero, i.e. P(N) = 0. In the same way as
%   POLYFIT, the coefficients, P(1:N-1), fit the data Y best in the least-
%   squares sense. You can also use Y = POLYVAL(PZERO,X) to evaluate the
%   polynomial because the output is the same as POLYFIT.
%
%   [P,S,MU] = POLYFITZERO() Return structure, S, similar to POLYFIT for use
%   with POLYVAL to calculate error estimates of predictions with P.
%
%   [P,S,MU] = POLYFITZERO() Scale X by std(X), returns MU = [0,std(X)].
%
%   See also POLYVAL, POLYFIT
%
%   Also see <a href="http://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn">POLYFITN by John D'Errico</a>
%
% Copyright (c) 2013 Mark Mikofski
% Version 1-1, 2013-10-15
%   add delta output
%   center and scale
% Version 1-0, 2011-06-29
%% check args
% X & Y should be numbers
assert(isnumeric(x) && isnumeric(y),'polyfitZero:notNumeric', ...
    'X and Y must be numeric.')
dim = numel(x); % number of elements in X
% DEGREE should be scalar positive number between 1 & 10 inclusive
assert(isnumeric(degree) && isscalar(degree) && degree>0 && degree<=10, ...
    'polyfitZero:degreeOutOfRange', ...
    'DEGREE must be an integer between 1 and 10.')
% DEGREE must be less than number of elements in X & Y
assert(degree<dim && degree==round(degree), ...
    'polyfitZero:DegreeGreaterThanDim', 'DEGREE must be less than numel(X)')
% X & Y should be same size vectors
assert(isvector(x) && isvector(y) && dim==numel(y), ...
    'polyfitZero:vectorMismatch', 'X and Y must be vectors of the same length.')
%% solve
% convert X & Y to column vectors
x = x(:); y = y(:);
% Scale X.
% attribution: this is based on code from POLYFIT by The MathWorks Inc.
if nargout > 2
   mu = [0; std(x)];
   x = (x - mu(1))/mu(2);
end
% using pow() is actually as fast or faster than looping, same # of flops!
z = zeros(dim,degree);
for n = 1:degree
    z(:,n) = x.^(degree-n+1);
end
p = z\y; % solve
p = [p;0]; % set y-intercept to zero
%% error estimates
% attribution: this is based on code from POLYFIT by The MathWorks Inc.
if nargout > 1
    V = [z,ones(dim,1)]; % append constant term for Vandermonde matrix
    % Return upper-triangular factor of QR-decomposition for error estimates
    R = triu(qr(V,0));
    r = y - V*p;
    S.R = R(1:size(R,2),:);
    S.df = max(0,length(y) - (degree+1));
    S.normr = norm(r);
end
p = p'; % polynomial output is row vector by convention
end


