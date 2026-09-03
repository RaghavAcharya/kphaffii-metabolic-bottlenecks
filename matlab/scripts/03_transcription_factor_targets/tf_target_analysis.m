
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
%% 
all_tf_targets = readtable("all_tf_targets.xlsx");
str_tf_targets = all_tf_targets(all_tf_targets.Str_Eng_Target,:);

% Script to create network visualizations from str_tf_targets table
% for each unique environment

% Check if table exists in workspace
if ~exist('str_tf_targets', 'var')
    error('str_tf_targets table not found in workspace');
end

% Get unique environments
environments = unique(str_tf_targets.Environment);
numEnvironments = length(environments);

% Create a figure to arrange all subplots
figure('Position', [100, 100, 1200, 800]);

% Set up color maps
uniqueTFCategories = unique(str_tf_targets.TF_Category);
uniqueGeneCategories = unique(str_tf_targets.Category);

% Create colormap for TF categories (edge colors)
edgeColorMap = lines(length(uniqueTFCategories));
% Store the actual RGB values, not cells
tfCategoryColors = containers.Map(uniqueTFCategories, mat2cell(edgeColorMap, ones(length(uniqueTFCategories), 1), 3));

% Create colormap for gene categories (node colors)
nodeColorMap = jet(length(uniqueGeneCategories));
% Store the actual RGB values, not cells
geneCategoryColors = containers.Map(uniqueGeneCategories, mat2cell(nodeColorMap, ones(length(uniqueGeneCategories), 1), 3));

% Calculate subplot grid dimensions
gridRows = ceil(sqrt(numEnvironments));
gridCols = ceil(numEnvironments/gridRows);

% Process each environment
for envIdx = 1:numEnvironments
    currentEnv = environments{envIdx};
    
    % Filter data for current environment
    envData = str_tf_targets(strcmp(str_tf_targets.Environment, currentEnv), :);
    
    % Create directed graph
    G = digraph();
    
    % Collect all unique nodes (genes and TFs)
    allGenes = unique(envData.Gene);
    allTFs = unique(envData.Transcription_Factor);
    allNodes = unique([allGenes; allTFs]);
    
    % Add nodes
    G = addnode(G, allNodes);
    
    % Add edges
    sources = cell(height(envData), 1);
    targets = cell(height(envData), 1);
    
    for i = 1:height(envData)
        sources{i} = envData.Transcription_Factor{i};
        targets{i} = envData.Gene{i};
    end
    
    G = addedge(G, sources, targets);
    
    % Calculate node properties
    nodeSizes = zeros(length(allNodes), 1);
    nodeColors = cell(length(allNodes), 1);
    
    % Process gene nodes
    for i = 1:length(allGenes)
        geneName = allGenes{i};
        nodeIdx = findnode(G, geneName);
        
        % Get gene category for color
        geneRows = envData(strcmp(envData.Gene, geneName), :);
        if ~isempty(geneRows)
            category = geneRows.Category{1};
            nodeColors{nodeIdx} = geneCategoryColors(category);
            
            % Set size based on Gene_Predicted_Change and environment type
            % Handle different possible data types for Gene_Predicted_Change
            predictedChange = geneRows.Gene_Predicted_Change(1);
            
            % Convert to numeric if it's a string or cell
            if iscell(predictedChange)
                predictedChange = str2double(predictedChange{1});
            elseif isstring(predictedChange)
                predictedChange = str2double(predictedChange);
            elseif istable(predictedChange)
                predictedChange = table2array(predictedChange);
            end
            
            if startsWith(currentEnv, 'HSA')
                % For HSA environments: larger as value increases
                nodeSizes(nodeIdx) = 20 + abs(predictedChange) * 10;
            else
                % For All environments: larger as value decreases
                nodeSizes(nodeIdx) = 20 + (1/max(abs(predictedChange), 0.01)) * 10;
            end
        end
    end
    
    % Process TF nodes
    for i = 1:length(allTFs)
        tfName = allTFs{i};
        nodeIdx = findnode(G, tfName);
        
        % Color TF nodes gray
        nodeColors{nodeIdx} = [0.7 0.7 0.7];
        
        % Size based on number of connected genes
        connectedGenes = length(find(strcmp(envData.Transcription_Factor, tfName)));
        nodeSizes(nodeIdx) = 10 + connectedGenes * 3;
    end
    
    % Create subplot for current environment
    subplot(gridRows, gridCols, envIdx);
    
    % Plot graph
    hPlot = plot(G, 'Layout', 'force');
    
    % Apply node sizes and colors 
    for i = 1:length(allNodes)
        highlight(hPlot, i, 'NodeColor', nodeColors{i}, 'MarkerSize', nodeSizes(i));
    end
    
    % Apply edge colors individually
    for i = 1:height(envData)
        source = envData.Transcription_Factor{i};
        target = envData.Gene{i};
        tfCategory = envData.TF_Category{i};
        
        % Get color from the map - should be a 1x3 RGB array already
        edgeColor = tfCategoryColors(tfCategory);
        
        highlight(hPlot, source, target, 'EdgeColor', edgeColor, 'LineWidth', 1.5);
    end
    
    % Add title and improve layout
    title(['Network for Environment: ', currentEnv]);
    
    % Add node labels
    hPlot.NodeLabel = G.Nodes.Name;
    
    % Adjust layout for better visualization
    layout(hPlot, 'force', 'WeightEffect', 'inverse');
end

% Create legend for TF categories (edge colors)
if numEnvironments < gridRows * gridCols
    subplot(gridRows, gridCols, numEnvironments + 1);
    hold on;
    legendHandles = zeros(length(uniqueTFCategories), 1);
    for i = 1:length(uniqueTFCategories)
        legendHandles(i) = plot([0, 1], [i, i], 'Color', tfCategoryColors(uniqueTFCategories{i}), 'LineWidth', 2);
    end
    hold off;
    axis off;
    legend(legendHandles, uniqueTFCategories, 'Location', 'best', 'Title', 'TF Categories');

    % Create legend for gene categories (node colors)
    if numEnvironments + 1 < gridRows * gridCols
        subplot(gridRows, gridCols, numEnvironments + 2);
        hold on;
        legendHandles = zeros(length(uniqueGeneCategories), 1);
        for i = 1:length(uniqueGeneCategories)
            legendHandles(i) = scatter(1, i, 100, geneCategoryColors(uniqueGeneCategories{i}), 'filled');
        end
        hold off;
        axis off;
        legend(legendHandles, uniqueGeneCategories, 'Location', 'best', 'Title', 'Gene Categories');
    end
end

% Add explanation for node sizes
annotation('textbox', [0.6, 0.05, 0.35, 0.15], 'String', {...
    'Node Size Legend:', ...
    '- Gene nodes: Based on Gene_Predicted_Change', ...
    '  (HSA env: larger = higher value, All env: larger = lower value)', ...
    '- TF nodes: Based on number of connected genes'}, ...
    'EdgeColor', 'none');

% Adjust figure properties
sgtitle('Transcription Factor - Gene Networks by Environment');