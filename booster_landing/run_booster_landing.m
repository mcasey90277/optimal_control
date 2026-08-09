function R = run_booster_landing(cfg)
% RUN_BOOSTER_LANDING  Front door: full 3-DOF booster-landing campaign.
%
%   solve (colloc + convex) -> certify G1-G5 -> TVLQR -> closed loop ->
%   Monte Carlo -> plots + movie -> results/booster_run.mat
%
% Run me with no arguments for the nominal campaign:
%   /Applications/MATLAB_R2025b.app/bin/matlab -batch ...
%     "cd('~/Desktop/optimal_control/booster_landing'); run_booster_landing"
%
% INPUTS:
%   cfg - (optional) .P params overrides (field-merged onto
%         booster_params()), .doMovie [true], .doMC [true], .Nrun [200],
%         .outdir [results/]
% OUTPUTS:
%   R   - everything: .P .solC .solV .rep .ctrl .out0 .mc .when
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] docs/superpowers/plans/2026-08-08-booster-landing.md
setup_paths;

%% ---------------- ADJUSTABLE PARAMETERS (defaults) ----------------
def = struct('doMovie', true, 'doMC', true, 'Nrun', 200, ...
             'outdir', fullfile(fileparts(mfilename('fullpath')), 'results'));
%% ------------------------------------------------------------------
if nargin < 1, cfg = struct(); end
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end
P = booster_params();
if isfield(cfg, 'P')
    pf = fieldnames(cfg.P);
    for k = 1:numel(pf), P.(pf{k}) = cfg.P.(pf{k}); end
end
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end

fprintf('=== [1/5] Guidance: collocation (N=%d) ===\n', P.N);
R.solC = solve_pdg_colloc(P);
fprintf('    tf=%.3f s  mf=%.2f kg\n', R.solC.tf, R.solC.mf);

fprintf('=== [2/5] Guidance: lossless convexification (golden tf) ===\n');
R.solV = solve_pdg_convex(P);
fprintf('    tf=%.3f s  mf=%.2f kg  gap=%.2e\n', ...
        R.solV.tf, R.solV.mf, R.solV.lossless_gap);

fprintf('=== [3/5] Certification ===\n');
R.rep = certify_pdg(R.solC, R.solV, P);
print_certify_report(R.rep);

fprintf('=== [4/5] Tracking: TVLQR + closed loop%s ===\n', ...
        ternary(cfg.doMC, ' + Monte Carlo', ''));
R.ctrl = tvlqr_design(R.solC, P);
R.out0 = sim_closed_loop(R.solC, R.ctrl, P, struct());
if cfg.doMC
    R.mc = run_monte_carlo(R.solC, R.ctrl, P, struct('Nrun', cfg.Nrun));
end

fprintf('=== [5/5] Products -> %s ===\n', cfg.outdir);
plot_pdg_solution(R.solC, R.solV, fullfile(cfg.outdir, 'pdg_solution.png'));
if cfg.doMC
    plot_footprint(R.mc, P, fullfile(cfg.outdir, 'footprint.png'));
end
if cfg.doMovie
    movie_landing(R.out0, R.solC, P, fullfile(cfg.outdir, 'landing.mp4'));
end
R.P = P;  R.when = datetime('now');
save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

fprintf('\n==================== SUMMARY ====================\n');
fprintf('tf        %.3f s      fuel  %.1f kg (mf %.1f)\n', ...
        R.solC.tf, P.m0 - R.solC.mf, R.solC.mf);
fprintf('gates     %s\n', ternary(R.rep.all_pass, 'ALL PASS', 'FAILURES -- see table'));
fprintf('nom miss  %.2f m  vtd %.2f m/s  (landed=%d, stop=%s)\n', ...
        R.out0.td.miss, R.out0.td.vtd, R.out0.td.landed, R.out0.td.stop);
if cfg.doMC
    fprintf('MC        %d runs, success %.1f%%  (landed %d / arrest %d / horizon %d)\n', ...
            cfg.Nrun, 100*R.mc.success_rate, R.mc.n_landed, R.mc.n_arrest, R.mc.n_horizon);
end
fprintf('design    guidance thrust de-rate etaT=%.2f (ceiling %.2f kN of Tmax=%.2f kN) --\n', ...
        P.etaT, P.etaT*P.Tmax/1e3, P.Tmax/1e3);
fprintf('          adjudicated 2026-08-09 so the tracker keeps thrust headroom against\n');
fprintf('          dispersions (see booster_params.m P.etaT note).\n');
fprintf('          terminal vf=%.2f m/s (not 0) -- v(tf)=0 is singular under Tmin>weight;\n', P.vf(3));
fprintf('          legs absorb touchdown speed up to the %.1f m/s vtd_max gate\n', P.vtd_max);
fprintf('          (adjudicated 2026-08-08, see booster_params.m P.vf note).\n');
fprintf('=================================================\n');
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
