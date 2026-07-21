function results = run_tf_sweep(tfFactors, maxIter, epsSched)
% RUN_TF_SWEEP  Map the minimum-fuel Delta-V vs transfer-time Pareto front.
%
% Fixed-endpoint front: the GTO departure and the tulip rendezvous point are
% held FIXED; only the transfer time t_f (imposed through the carried time
% state, t(tau_f)=t_f) is varied. For each t_f we run the energy->fuel homotopy
% and record the sharp bang-bang point.
%
% Continuation is done on the SMOOTH (energy, eps=1) solution, never on the
% bang-bang one: warm-starting a bang-bang solve from a bang-bang solution at a
% shifted t_f forces coast arcs to become burns (an active-set change that
% wedges IPOPT). Instead we carry the eps=1 energy solution across t_f (smooth,
% so continuation is robust) and re-sharpen eps:1->0 at each t_f. The tau-domain
% and mesh are held fixed; only the warm-start time state is rescaled by
% t_f^new/t_f^prev.
%
% Expected shape: monotone-decreasing dV(t_f), steep near t_f^min (where
% min-fuel -> min-time, ~4.4665 km/s, 0 switches) and flattening for large t_f.
%
% INPUTS:
%   tfFactors - t_f / t_f^min values to solve [default below]
%   maxIter   - IPOPT max iters per solve [default 1200]
%
% OUTPUTS:
%   results - struct array: .factor .tf_days .dV .prop_kg .switches .edge
%             .defect .success  (saved to tf_sweep_results.mat + a PNG)

here = fileparts(mfilename('fullpath'));  addpath(here);
if nargin<1||isempty(tfFactors), tfFactors=[1.00 1.05 1.10 1.15 1.23 1.33 1.45 1.60 1.80]; end
if nargin<2||isempty(maxIter),   maxIter=1200; end
% per-t_f re-sharpen schedule. Default is coarse (fast shape preview, ~96% edge,
% dV ~1.5% above optimum); pass a finer schedule for the accurate front, e.g.
% [1 0.6 0.3 0.15 0.08 0.04 0.02 0.01 0.005 0.002] (edge ~99.5%, ~25 switches).
if nargin<3||isempty(epsSched),  epsSched=[1 0.4 0.15 0.05 0.02]; end
pSund = 1.5;

p = cr3bp_lt_params(0.025, 15, 2100);

% energy-seed warm start (collocation-feasible time-mesh solution at 1.15x)
S = load(fullfile(here, 'minfuel_from_energy_seed.mat'));
rv0 = S.rv0;  rvf = S.rvf;
[sigma, Xseed, Useed, tauf0] = sundman_seed_map(S.nlp.X, S.nlp.U, S.tf, S.sigma, ...
                                                pSund, p.muStar, rv0, rvf);
tfAnchor = S.tf;   tfMin = tfAnchor/1.15;

facs  = unique(sort(tfFactors(:).'));
downF = flip(facs(facs <= 1.15+1e-9));    % anchor first, then descending
upF   = facs(facs > 1.15+1e-9);           % ascending above the anchor

results = struct('factor',{},'tf_days',{},'dV',{},'prop_kg',{},'switches',{}, ...
                 'edge',{},'defect',{},'success',{},'tf',{},'X',{},'U',{}, ...
                 'lamDef',{},'primerAlignDeg',{},'lamMassEnd',{});

% ---- down pass (from the energy seed at anchor, descending) ----
savePath = fullfile(here,'tf_sweep_results.mat');
Xe = Xseed;  Ue = Useed;  tfPrev = tfAnchor;
for f = downF
    [rec, Xe, Ue, tfPrev] = solve_tf(f, tfMin, tfPrev, Xe, Ue, sigma, tauf0, ...
                                     rv0, rvf, pSund, p, epsSched, maxIter);
    results(end+1) = rec; %#ok<AGROW>
    save(savePath, 'results', 'tfMin');            % incremental
end
% ---- up pass (re-seed from the energy seed at anchor, ascending) ----
Xe = Xseed;  Ue = Useed;  tfPrev = tfAnchor;
for f = upF
    [rec, Xe, Ue, tfPrev] = solve_tf(f, tfMin, tfPrev, Xe, Ue, sigma, tauf0, ...
                                     rv0, rvf, pSund, p, epsSched, maxIter);
    results(end+1) = rec; %#ok<AGROW>
    save(savePath, 'results', 'tfMin');            % incremental
end

[~,ord] = sort([results.factor]);  results = results(ord);
save(fullfile(here,'tf_sweep_results.mat'), 'results', 'tfMin');

fprintf('\n=== dV-TIME FRONT (factor, days, dV, prop, switches, primer-align, ok) ===\n');
for k = 1:numel(results)
    r = results(k);
    fprintf('  %.2f  %6.2f d  %7.4f km/s  %7.4f kg  %3d sw  primer=%.2f deg  %d\n', ...
            r.factor, r.tf_days, r.dV, r.prop_kg, r.switches, r.primerAlignDeg, r.success);
end
plot_front(results, tfMin, p.tStar, here);
end

% -------------------------------------------------------------------------
function [rec, Xe, Ue, tfOut] = solve_tf(f, tfMin, tfPrev, Xe, Ue, sigma, tauf0, ...
                                         rv0, rvf, pSund, p, epsSched, maxIter)
% Solve at t_f=f*tfMin by re-sharpening the (continued) energy warm start
% (Xe,Ue). Returns the front record and the new eps=1 energy solution for the
% next t_f. The tau-domain is fixed; the warm-start time state is rescaled.
tf = f*tfMin;
Xe(8,:) = Xe(8,:) * (tf/tfPrev);
fprintf('\n==== t_f = %.3f (%.3fx = %.2f d) ====\n', tf, f, tf*p.tStar/86400);
X0 = Xe;  U0 = Ue;  best = [];  energy = [];
for ie = 1:numel(epsSched)
    e = epsSched(ie);
    % The first (largest-eps, ~energy) solve is a genuine continuation move to
    % the new t_f -> use the LOOSE warm start so IPOPT can explore; the
    % subsequent sharpening solves re-solve AT a near-bang-bang point -> TIGHT.
    tight = (ie > 1);
    o = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                               X0, U0, tauf0, pSund, maxIter, e, tight);
    ok = o.success && o.maxDefect < 1e-6;
    if ie==1 && ok, energy = o; end          % eps=1 solution for continuation
    if ok, X0 = o.X;  U0 = o.U;  best = o; end
end
if isempty(best),  best = o;  end
if isempty(energy), energy = struct('X',Xe,'U',Ue); end   % fall back: keep prior
Xe = energy.X;  Ue = energy.U;  tfOut = tf;
dV = p.c*log(1/best.mf)*p.lStar/p.tStar;
adv = best.success && best.maxDefect < 1e-6;
rec = struct('factor',f, 'tf_days',tf*p.tStar/86400, 'dV',dV, ...
             'prop_kg',p.m0kg*(1-best.mf), 'switches',best.switches, ...
             'edge',best.edge, 'defect',best.maxDefect, 'success',adv, ...
             'tf',tf, 'X',best.X, 'U',best.U, ...          % full state+control
             'lamDef',best.lamDef, 'primerAlignDeg',best.primerAlignDeg, ...
             'lamMassEnd',best.lamMassEnd);                % discrete costates + PMP check
fprintf('  -> dV=%.4f km/s  prop=%.4f kg  switches=%d  edge=%.1f%%  defect=%.2g  %s\n', ...
        dV, p.m0kg*(1-best.mf), best.switches, 100*best.edge, best.maxDefect, ...
        ternary(adv,'OK','(loose)'));
end

% -------------------------------------------------------------------------
function plot_front(results, tfMin, tStar, here)
ok = [results.success]==1;
d  = [results.tf_days];  v = [results.dV];  sw = [results.switches];
fig = figure('Color','w','Position',[100 100 760 560],'Visible','off');
try
    theme(fig,'light');
catch
end
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(ax1, d(ok), v(ok), '-o', 'Color',[0.60 0.10 0.10], 'MarkerFaceColor',[0.60 0.10 0.10], 'LineWidth',1.8);
if any(~ok), plot(ax1, d(~ok), v(~ok), 'x', 'Color',[0.5 0.5 0.5], 'MarkerSize',9); end
plot(ax1, tfMin*tStar/86400, 4.4665, 'ks', 'MarkerFaceColor','k', 'MarkerSize',7);
xlabel(ax1,'transfer time t_f (days)'); ylabel(ax1,'\DeltaV (km/s)');
title(ax1,'Minimum-fuel \DeltaV--time Pareto front (GTO \rightarrow tulip, fixed endpoints)');
text(ax1, tfMin*tStar/86400, 4.4665, '  min-time (4.4665, 0 switches)', 'FontSize',9, 'Color',[0.2 0.2 0.2]);
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
plot(ax2, d(ok), sw(ok), '-s', 'Color',[0.15 0.35 0.75], 'MarkerFaceColor',[0.15 0.35 0.75], 'LineWidth',1.6);
xlabel(ax2,'transfer time t_f (days)'); ylabel(ax2,'bang-bang switches');
title(ax2,'switch count vs transfer time');
exportgraphics(fig, fullfile(here,'tf_dv_front.png'), 'Resolution',150);
close(fig);
fprintf('WROTE %s\n', fullfile(here,'tf_dv_front.png'));
end

function y = ternary(c,a,b)
if c
    y = a;
else
    y = b;
end
end
