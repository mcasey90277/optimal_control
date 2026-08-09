% TEST_MONTE_CARLO_SMALL  20-run MC sanity: executes, deterministic under
% the seed (two calls identical), success rate positive, outputs shaped.
%
% TASK-7 EVOLUTION: mc.ok requires out.td.landed in addition to miss/vtd/m
% (an arrest, depletion or horizon run is a failure regardless of how
% optimistic its vtd looks -- see sim_closed_loop's ADAPTATION note). This
% test also checks the failure-mode breakdown fields are present and add
% up to Nrun -- now including n_depleted, split out of n_arrest by the
% 2026-08-09 external code review.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
% NOTE: N=30 is a COARSE smoke grid, chosen only to keep this test fast --
% it is NOT the production grid (P.N=60) and its ~65% success_rate is NOT
% representative of the shipped campaign's performance. Do not read this
% test's success rate as a robustness claim. For the production-grid
% number see process/task-8-report.md (50-run probe at N=P.N) and the
% README's flagship table (200 runs). Corrected 2026-08-09: this note used
% to point at "run_monte_carlo.m's own probe at N=P.N", which does not
% exist -- the probe was a one-off run recorded in the task-8 report, not
% code that ships in that function.
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);
mc1  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
mc2  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
assert(isequal(size(mc1.land), [20 2]) && numel(mc1.vtd) == 20);
assert(isequal(mc1, mc2), 'MC not deterministic under seed (full-struct compare)');
assert(mc1.success_rate > 0, 'zero successes at default dispersions');

% Task-7 evolution: landed-gated ok, failure-mode breakdown.
assert(numel(mc1.landed) == 20 && islogical(mc1.landed));
assert(all(mc1.ok(~mc1.landed) == false), ...
       'ok=true on a non-landed (arrest/horizon) run -- landed gate not applied');
assert(mc1.n_landed + mc1.n_arrest + mc1.n_depleted + mc1.n_horizon == 20, ...
       'failure-mode counts do not sum to Nrun');
assert(mc1.n_landed == nnz(mc1.landed), 'n_landed disagrees with mc1.landed');
% Every stop string must be one of the four documented classes -- the
% counts summing to Nrun would also be satisfied by a miscount that
% happened to cancel (external code review follow-up, 2026-08-09: the
% 'depleted' class existed as an EVENT long before it existed as a LABEL,
% and a test that only checks a sum cannot notice that).
assert(all(ismember(mc1.stop, {'touchdown','arrest','depleted','horizon'})), ...
       'mc.stop contains a class outside {touchdown,arrest,depleted,horizon}');
% (mc1 vs mc2 equality of the breakdown counts is already covered by the
% full-struct isequal(mc1, mc2) assertion above.)

% Bad-input probe: opts.sig validation must reject a typo'd field name and
% a wrong-shaped value, cleanly, before any run executes.
badField = false;
try
    run_monte_carlo(sol, ctrl, P, struct('Nrun', 1, 'sig', struct('thrust_sig', 0.02)));
catch e_
    badField = strcmp(e_.identifier, 'run_monte_carlo:unknownSigField');
    fprintf('  bad-field probe error (expected): %s\n', e_.message);
end
assert(badField, 'unknown opts.sig field was not rejected');

badShape = false;
try
    run_monte_carlo(sol, ctrl, P, struct('Nrun', 1, 'sig', struct('r0', [100 100 50])));
catch e_
    badShape = strcmp(e_.identifier, 'run_monte_carlo:badSigShape');
    fprintf('  bad-shape probe error (expected): %s\n', e_.message);
end
assert(badShape, 'wrong-shaped opts.sig.r0 (1x3) was not rejected');

fprintf('test_monte_carlo_small PASS  success=%.0f%%  landed=%d arrest=%d depleted=%d horizon=%d\n', ...
        100*mc1.success_rate, mc1.n_landed, mc1.n_arrest, mc1.n_depleted, mc1.n_horizon);
