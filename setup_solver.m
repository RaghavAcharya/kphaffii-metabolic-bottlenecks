function setup_solver()
% SETUP_SOLVER  Configure solver/toolbox paths for the Paper 1 pipeline.
%   Each analysis script calls this after adding matlab/functions to the path.
%   EDIT the lines below to match your installation, then remove the warning.
%
%   Requirements: MATLAB R2024a, Gurobi, COBRA Toolbox v3.0, RAVEN 2.0.

% --- Gurobi (LP/MILP solver) -------------------------------------------------
% run('C:\gurobi1202\win64\matlab\gurobi_setup');           % Windows example
% run('/opt/gurobi/linux64/matlab/gurobi_setup.m');         % Linux example

% --- COBRA Toolbox v3.0 ------------------------------------------------------
% addpath('/path/to/cobratoolbox'); initCobraToolbox(false);
% changeCobraSolver('gurobi','all');

% --- RAVEN 2.0 (provides importModel/readYAMLmodel etc.) ---------------------
% addpath(genpath('/path/to/RAVEN'));

warning(['setup_solver.m has not been configured. Edit ' ...
    'setup_solver.m to point to your Gurobi, COBRA v3.0, and RAVEN 2.0 installs.']);
end
