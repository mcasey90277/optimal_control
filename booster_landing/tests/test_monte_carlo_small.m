% TEST_MONTE_CARLO_SMALL  20-run MC sanity: executes, deterministic under
% the seed (two calls identical), success rate positive, outputs shaped.
%
% TASK-7 EVOLUTION: mc.ok requires out.td.landed in addition to miss/vtd/m
% (an arrest or horizon run is a failure regardless of how optimistic its
% vtd looks -- see sim_closed_loop's ADAPTATION note). This test also
% checks the failure-mode breakdown fields are present and add up to Nrun.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);
mc1  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
mc2  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
assert(isequal(size(mc1.land), [20 2]) && numel(mc1.vtd) == 20);
assert(isequal(mc1.land, mc2.land), 'MC not deterministic under seed');
assert(mc1.success_rate > 0, 'zero successes at default dispersions');

% Task-7 evolution: landed-gated ok, failure-mode breakdown.
assert(numel(mc1.landed) == 20 && islogical(mc1.landed));
assert(all(mc1.ok(~mc1.landed) == false), ...
       'ok=true on a non-landed (arrest/horizon) run -- landed gate not applied');
assert(mc1.n_landed + mc1.n_arrest + mc1.n_horizon == 20, ...
       'failure-mode counts do not sum to Nrun');
assert(mc1.n_landed == nnz(mc1.landed), 'n_landed disagrees with mc1.landed');
assert(isequal(mc1.n_landed, mc2.n_landed) && isequal(mc1.n_arrest, mc2.n_arrest) ...
       && isequal(mc1.n_horizon, mc2.n_horizon), 'failure-mode counts not deterministic');

fprintf('test_monte_carlo_small PASS  success=%.0f%%  landed=%d arrest=%d horizon=%d\n', ...
        100*mc1.success_rate, mc1.n_landed, mc1.n_arrest, mc1.n_horizon);
