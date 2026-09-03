function verify_matlab_dependencies()
% VERIFY_MATLAB_DEPENDENCIES  Authoritative transfer check for the Paper 1 repo.
%   Run in MATLAB from anywhere (it locates the repo root itself). For each entry
%   script it computes the true set of required files via
%   matlab.codetools.requiredFilesAndProducts, then reports:
%     (a) required files that live in the repo and are PRESENT,
%     (b) required files that are MISSING (a real transfer gap),
%     (c) required MathWorks products / external toolboxes (COBRA/RAVEN/Gurobi),
%     (d) bundled functions in matlab/functions that were NOT required by any
%         entry script (safe-to-prune candidates).
%
%   Note: requiredFilesAndProducts resolves what is on the MATLAB path. Add COBRA
%   and RAVEN to your path first so their functions are classified as external
%   products rather than reported as unresolved.

    here     = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(here);
    funcDir  = fullfile(repoRoot,'matlab','functions');
    addpath(genpath(funcDir));

    entry = {
        fullfile(repoRoot,'matlab','scripts','01_HSA_target_identification','HSA_model_Generation_Verification.m')
        fullfile(repoRoot,'matlab','scripts','02_common_target_identification','common_target_identification_validation.m')
        fullfile(repoRoot,'matlab','scripts','03_transcription_factor_targets','tf_target_analysis.m')
        fullfile(repoRoot,'matlab','scripts','04_nutrient_cost','metCostDetermination.m')
        fullfile(repoRoot,'matlab','scripts','04_nutrient_cost','metCostKOBenchmark.m')
    };

    usedRepoFiles = strings(0,1);
    fprintf('== Dependency check ==\n');
    for i = 1:numel(entry)
        s = entry{i};
        [~, nm, ex] = fileparts(s);
        fprintf('\n--- %s%s ---\n', nm, ex);
        if ~isfile(s)
            fprintf('  MISSING ENTRY SCRIPT: %s\n', s); continue;
        end
        try
            [flist, plist] = matlab.codetools.requiredFilesAndProducts(s);
        catch err
            fprintf('  could not analyze (%s). Ensure COBRA/RAVEN are on the path.\n', err.message);
            continue;
        end
        for k = 1:numel(flist)
            f = flist{k};
            if startsWith(f, repoRoot)
                rel = erase(f, [repoRoot filesep]);
                if isfile(f)
                    fprintf('  present : %s\n', rel);
                    usedRepoFiles(end+1,1) = string(f); %#ok<AGROW>
                else
                    fprintf('  MISSING : %s\n', rel);
                end
            end
        end
        fprintf('  external products: %s\n', strjoin({plist.Name}, ', '));
    end

    % prunable: bundled functions never required by any entry script
    bundled = dir(fullfile(funcDir,'**','*.m'));
    fprintf('\n== Bundled functions not required by any entry script (prune candidates) ==\n');
    none = true;
    for i = 1:numel(bundled)
        full = string(fullfile(bundled(i).folder, bundled(i).name));
        if ~any(usedRepoFiles == full)
            fprintf('  %s\n', erase(char(full), [funcDir filesep]));
            none = false;
        end
    end
    if none, fprintf('  (none - every bundled function is required)\n'); end
    fprintf('\nDone. Any line marked MISSING is a genuine transfer gap to fix.\n');
end
