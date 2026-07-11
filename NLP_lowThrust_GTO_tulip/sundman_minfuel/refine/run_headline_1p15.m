function history = run_headline_1p15()
% RUN_HEADLINE_1P15  Prototype demonstration on the certified 1.15x solution.
%
% Prepares a duals-carrying seed from sundman_minfuel_certified.mat, then runs
% the refinement loop and prints the summary table for RESULTS.md.
%
% INPUTS:  none
% OUTPUTS: history - the refine_loop history struct array
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md

here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   % ms_band/setup_paths does NOT add ms_band itself
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();  cd(old);

src  = fullfile(here, '..', 'sundman_minfuel_certified.mat');
seed = fullfile(here, 'seed_1p15.mat');
if ~isfile(seed), prep_refine_seed(src, seed); end

opts = struct('maxRounds', 4, 'tag', 'headline_1p15', 'K', 8, 'maxAdd', 40);
history = refine_loop(seed, opts);

fprintf('\n=== HEADLINE 1.15x SUMMARY ===\n');
fprintf('%-6s %-7s %-4s %-11s %-11s %-7s %-11s\n', ...
        'round', 'nodes', 'sw', 'maxMove', 'dProp(kg)', 'nViol', 'HresMax');
for r = 1:numel(history)
    h = history(r);
    fprintf('%-6d %-7d %-4d %-11.2e %-11.2e %-7d %-11.2e\n', ...
            r-1, h.nNodes, h.switches, h.maxSwitchMove, h.dProp, h.nViol, h.HresMax);
end
end
