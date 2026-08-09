function R = run_booster_landing(cfg)
% RUN_BOOSTER_LANDING  Front door: full 3-DOF booster-landing campaign.
%
%   solve (colloc + convex) -> certify G1-G5 -> TVLQR -> closed loop ->
%   Monte Carlo -> plots + movie -> results/booster_run.mat
%   -> optional Phase 2 (cfg.phase2): drag-on re-solve (warm-started from
%      the Phase-1 vacuum solution) -> certify (G1/G2/G5 only, no convex
%      twin) -> drag-aware TVLQR + closed loop + wind Monte Carlo ->
%      vacuum-vs-drag comparison + footprint plots
%
% Run me with no arguments for the nominal (Phase 1, vacuum) campaign:
%   /Applications/MATLAB_R2025b.app/bin/matlab -batch ...
%     "cd('~/Desktop/optimal_control/booster_landing'); run_booster_landing"
% Add Phase 2 (atmosphere + drag) with:
%   run_booster_landing(struct('phase2', true))
%
% INPUTS:
%   cfg - (optional) .P params overrides (field-merged onto
%         booster_params()), .doMovie [true], .doMC [true], .Nrun [200],
%         .phase2 [false] (drag-on re-solve stage, see above),
%         .outdir [results/]. A cfg.P override that shadows a DERIVED
%         field's parent (P.Tmax -> P.Tmin=0.40*Tmax, P.g0 -> P.gvec)
%         is automatically re-derived unless the dependent is ALSO
%         overridden explicitly -- see the re-derive block below.
% OUTPUTS:
%   R   - everything: .P .solC .solV .rep .ctrl .out0 .mc .when, and, when
%         cfg.phase2 is true: .solD .repD .ctrlD .Pd (drag-on collocation
%         solution, its certify_pdg report, its TVLQR design, and the
%         drag-on params -- the same product KINDS as .solC/.rep/.ctrl/.P,
%         for the drag-on problem instead of vacuum), plus .mcD (wind
%         Monte Carlo, same shape as .mc) ONLY when cfg.doMC is also true
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md --
%       source of truth for every adjudicated number/date quoted below
%       (etaT, vf, G3 gate); this file only echoes them for a cold
%       console reader, never re-adjudicates.
%   [2] docs/superpowers/plans/2026-08-08-booster-landing.md
setup_paths;

%% ---------------- ADJUSTABLE PARAMETERS (defaults) ----------------
def = struct('doMovie', true, 'doMC', true, 'Nrun', 200, 'phase2', false, ...
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
    % Re-derive dependents a flat field-overwrite would otherwise leave
    % STALE (fix report, 2026-08-09 review): booster_params.m derives
    % P.Tmin = 0.40*P.Tmax and P.gvec = [0;0;-g0] from their parents. A
    % caller who overrides only the parent (e.g. cfg.P.Tmax) would
    % silently get an inconsistent (Tmax, Tmin) pair that certify_pdg
    % would certify without complaint. Re-derive whenever the parent was
    % overridden and the dependent was NOT also explicitly overridden
    % (an explicit dependent override always wins).
    if isfield(cfg.P, 'Tmax') && ~isfield(cfg.P, 'Tmin')
        P.Tmin = 0.40 * P.Tmax;
        fprintf('    [cfg.P] Tmax overridden -> re-derived Tmin = 0.40*Tmax = %.1f N\n', P.Tmin);
    end
    if isfield(cfg.P, 'g0') && ~isfield(cfg.P, 'gvec')
        P.gvec = [0; 0; -P.g0];
        fprintf('    [cfg.P] g0 overridden -> re-derived gvec = [0;0;-g0] = [%.4f; %.4f; %.4f]\n', P.gvec);
    end
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

% CHECKPOINT SAVE (fix report, 2026-08-09 review): solve + certify + MC
% is the expensive ~5-min part of the campaign. Save it BEFORE the viz
% stage, not only after -- a movie_landing throw (the H.264
% divisible-by-16 frame-size path is this campaign's known bite point,
% see that function's own header) used to discard the entire run with
% nothing on disk. The final save below re-writes the same file with
% R.when refreshed and any viz-stage products folded in.
R.P = P;  R.when = datetime('now');
fprintf('    [checkpoint] solve+certify+MC products -> %s\n', ...
        fullfile(cfg.outdir, 'booster_run.mat'));
save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

fprintf('=== [5/5] Products -> %s ===\n', cfg.outdir);
try
    plot_pdg_solution(R.solC, R.solV, fullfile(cfg.outdir, 'pdg_solution.png'));
    if cfg.doMC
        plot_footprint(R.mc, P, fullfile(cfg.outdir, 'footprint.png'));
    end
    if cfg.doMovie
        movie_landing(R.out0, R.solC, P, fullfile(cfg.outdir, 'landing.mp4'));
    end
catch vizErr
    warning('run_booster_landing:vizFailed', ...
        ['viz stage failed (%s) -- solve/certify/MC products are already ' ...
         'checkpointed at %s; continuing with whatever plots/movie ' ...
         'completed before the failure.'], ...
        vizErr.message, fullfile(cfg.outdir, 'booster_run.mat'));
end
R.when = datetime('now');
save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

%% ---------------- Phase 2: atmosphere + drag (opt-in) ----------------
% Drag-on re-solve, warm-started from the just-computed Phase-1 vacuum
% solution R.solC (P.drag.on is opt-in throughout the campaign -- see
% lib/pdg_dynamics.m -- so Phase 1 above ran entirely unmodified, in
% vacuum, before this block ever executes). No convex twin exists for the
% drag-on problem (lossless convexification is only exact in vacuum -- see
% certify_pdg.m's solV doc), so certify_pdg(R.solD, [], Pd) skips G3/G4;
% G5's primer-alignment sub-check is also LOOSENED (not removed -- kept as
% a real, 10-deg-threshold detector) under drag -- see that file's "G5
% primer LOOSENING under drag" note.
if cfg.phase2
    fprintf('=== [P2] Drag-on re-solve (warm-started) ===\n');
    Pd = P;  Pd.drag.on = true;
    R.solD = solve_pdg_colloc(Pd, struct('init', R.solC));
    fprintf('    tf=%.3f s  mf=%.2f kg\n', R.solD.tf, R.solD.mf);
    R.repD = certify_pdg(R.solD, [], Pd);
    print_certify_report(R.repD);
    ctrlD   = tvlqr_design(R.solD, Pd);
    % Saved for parity with Phase 1's R.ctrl/R.P (task-11 close-out review
    % minor): no flagship re-run needed for this to take effect, it
    % populates on the next full run.
    R.ctrlD = ctrlD;  R.Pd = Pd;
    % cfg.doMC GUARD (task-11 close-out review, Important-1): the drag
    % Monte Carlo + its footprint plot are the expensive, dispersion-only
    % products -- same thing cfg.doMC already gates in Phase 1 above. This
    % branch used to run both unconditionally regardless of cfg.doMC, the
    % canonical advertised-but-ignored-cfg-option bug (see
    % test_run_front_door.m's own header for why that bug class gets a
    % dedicated regression check).
    if cfg.doMC
        R.mcD = run_monte_carlo(R.solD, ctrlD, Pd, struct('Nrun', cfg.Nrun));
    end

    % CHECKPOINT SAVE (same discipline as Phase 1, task-9/10 fix report):
    % solve+certify+MC is the expensive part; save it before viz can throw.
    R.when = datetime('now');
    fprintf('    [checkpoint] Phase-2 solve+certify+MC products -> %s\n', ...
            fullfile(cfg.outdir, 'booster_run.mat'));
    save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

    try
        plot_vacuum_vs_drag(R.solC, R.solD, fullfile(cfg.outdir, 'phase2_vac_vs_drag.png'));
        if cfg.doMC
            plot_footprint(R.mcD, Pd, fullfile(cfg.outdir, 'phase2_footprint.png'));
        end
    catch vizErrD
        warning('run_booster_landing:phase2VizFailed', ...
            ['Phase-2 viz stage failed (%s) -- solve/certify/MC products are ' ...
             'already checkpointed at %s.'], ...
            vizErrD.message, fullfile(cfg.outdir, 'booster_run.mat'));
    end
    R.when = datetime('now');
    save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

    if cfg.doMC
        fprintf('P2: fuel vac %.1f kg -> drag %.1f kg; MC(wind) %.1f%%\n', ...
                P.m0 - R.solC.mf, Pd.m0 - R.solD.mf, 100*R.mcD.success_rate);
    else
        fprintf('P2: fuel vac %.1f kg -> drag %.1f kg (MC skipped, cfg.doMC=false)\n', ...
                P.m0 - R.solC.mf, Pd.m0 - R.solD.mf);
    end
end

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
if cfg.phase2
    fprintf('phase2    tf %.3f s (d%+.3f)  fuel %.1f kg (dfuel %+.1f)  gates %s\n', ...
            R.solD.tf, R.solD.tf - R.solC.tf, P.m0 - R.solD.mf, ...
            (P.m0 - R.solC.mf) - (P.m0 - R.solD.mf), ...
            ternary(R.repD.all_pass, 'ALL PASS', 'FAILURES -- see table'));
    if cfg.doMC
        fprintf('          MC(wind) %d runs, success %.1f%%  (landed %d / arrest %d / horizon %d)\n', ...
                cfg.Nrun, 100*R.mcD.success_rate, R.mcD.n_landed, R.mcD.n_arrest, R.mcD.n_horizon);
    else
        fprintf('          MC(wind) skipped (cfg.doMC=false)\n');
    end
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
