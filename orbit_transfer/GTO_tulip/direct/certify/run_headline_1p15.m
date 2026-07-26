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
% One campaign setup (2026-07-26 flatten): direct/setup_paths adds lib, certify,
% viz, cr3bp_common AND the indirect ms_band verifier. This used to cd into
% ms_band to borrow its setup_paths via '../../ms_band' -- a path that stopped
% resolving when ms_band moved to indirect/ on 2026-07-21.
addpath(fullfile(here, '..'));  setup_paths();

src  = fullfile(here, '..', 'lib', 'sundman_minfuel_certified.mat');
cfgH = minfuel_config();
seed = fullfile(cfgH.dirs.root, 'seed_1p15.mat');   % generated cache -> results/
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
