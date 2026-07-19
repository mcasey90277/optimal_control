% TEST_MEE_SEED  Seed sanity: monotone P growth, mass linearity in time,
% small trapezoid-stencil defect, stopP event capture, ballistic invariance.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
p = kepler_lt_params(10, 1500, 2000);

% (a)+(b) transverse constant-thrust seed: P monotone increasing, dL0 exact,
% mass consistent with mdot = -(Tm/c)*thr (linear in time, not L).
% NOTE nRev=3 (not the paper's optimal-solution 7.5): an always-on,
% full-throttle, fixed-RTN-direction burn is an open-loop (not a feedback)
% control law, so nothing stops it at GEO. Diagnosed empirically (see
% task-2-report.md): both 'transverse' and 'tangential' fixed-direction
% full-throttle burns overshoot P=1 around rev~3.1 and destabilize into a
% coordinate-singular (e>1) regime around rev~4.0-4.1, well short of the
% paper's 7.5-rev OPTIMAL bang-bang solution (which reaches GEO in 7.5 revs
% only because it switches thrust off along most of each orbit). nRev=3 stays
% safely inside the smooth, non-escaping window while still exercising every
% behavior the check cares about (monotone P growth well past its initial
% value, exact dL0, mass/time consistency, small stencil defect).
opts = struct('thr', 1, 'betaMode', 'transverse', 'N', 190, 'nRev', 3);
[sg, X0, U0, dL0, inf1] = mee_seed(p, opts);

assert(all(isfinite(sg(:))) && all(isfinite(X0(:))) && all(isfinite(U0(:))) ...
    && isfinite(dL0), 'seed outputs must be finite');
assert(isequal(size(sg), [191 1]), 'sigma size');
assert(isequal(size(X0), [7 191]), 'X0 size');
assert(isequal(size(U0), [4 191]), 'U0 size');
assert(all(diff(sg) > 0) && sg(1) == 0 && sg(end) == 1, 'sigma must be uniform 0->1');

assert(all(diff(X0(1,:)) > 0), 'P must strictly increase under transverse thrust');
assert(abs(dL0 - 2*pi*opts.nRev) < 1e-9, 'dL0 must match 2*pi*nRev to 1e-9');

assert(all(diff(X0(6,:)) < 0), 'mass must strictly decrease under thrust');
mdrop_expected = (p.Tmax/p.c) * opts.thr * inf1.tEnd;
mdrop_actual   = 1 - X0(6,end);
assert(abs(mdrop_actual - mdrop_expected) < 1e-6*mdrop_expected, ...
    'mass drop inconsistent with mdot = -(Tm/c)*thr: got %.6e vs %.6e', ...
    mdrop_actual, mdrop_expected);

% (c) trapezoid-stencil defect on the solver's sigma-grid (mirrors
% seed_stencil_defect in test_seed.m). NOTE the dL0 factor: lt_mee_rhs
% returns d/dL, the grid is sigma, and dL/dsigma = dL0.
dmax = seed_stencil_defect(sg, X0, U0, dL0, p);
assert(dmax < 1e-2, 'seed defect too big: %.2e', dmax);

% (d) stopP variant: integrate to GEO (P=1 in these nondim units, since
% LU_km = GEO radius), assert the event fired and info.nRev is finite.
opts2 = struct('thr', 1, 'betaMode', 'transverse', 'N', 100, 'stopP', 1.0);
[sg2, X02, ~, dL02, inf2] = mee_seed(p, opts2);
assert(abs(X02(1,end) - 1.0) < 1e-6, 'stopP event did not land on P=1 to 1e-6');
assert(isfinite(inf2.nRev) && inf2.nRev > 0, 'info.nRev must be finite and positive');
assert(all(diff(sg2) > 0), 'sigma2 must be strictly increasing');
assert(dL02 > 0, 'dL0 must be positive');

% (e) ballistic sanity: thr=0 -> elements frozen, tEnd ~= 2 Kepler periods.
opts3 = struct('thr', 0, 'betaMode', 'transverse', 'N', 50, 'nRev', 2);
[~, X03, ~, ~, inf3] = mee_seed(p, opts3);
elemSpan = max(X03(1:5,:), [], 2) - min(X03(1:5,:), [], 2);
assert(all(elemSpan < 1e-8), 'elements must stay frozen under zero thrust');
P0 = X03(1,1);  e0 = 0.75;
a0 = P0 / (1 - e0^2);
Tperiod = 2*pi*sqrt(a0^3 / p.mu);
assert(abs(inf3.tEnd - 2*Tperiod) < 0.01*2*Tperiod, ...
    'ballistic tEnd not ~2 Kepler periods: got %.6f vs %.6f', inf3.tEnd, 2*Tperiod);

% (f) [check 6] tangential betaMode: unit RTN velocity-direction steering,
% previously exercised by zero tests despite being implemented (Task 4
% depends on it). Same N/thr/nRev as check 2 (N=190, thr=1, nRev=3) so the
% comparison sits in the same stable, non-escaping window documented above
% (both 'transverse' and 'tangential' overshoot P=1 around rev~3.1).
opts4 = struct('thr', 1, 'betaMode', 'tangential', 'N', 190, 'nRev', 3);
[sg4, X04, U04, dL04, inf4] = mee_seed(p, opts4);

assert(all(isfinite(sg4(:))) && all(isfinite(X04(:))) && all(isfinite(U04(:))) ...
    && isfinite(dL04), 'tangential seed outputs must be finite');

betaNorm4 = sqrt(sum(U04(1:3,:).^2, 1));
assert(all(abs(betaNorm4 - 1) < 1e-12), 'tangential beta must be unit-norm to 1e-12');
assert(all(U04(3,:) == 0), 'tangential beta normal component must be exactly 0');
assert(all(U04(2,:) > 0), 'tangential beta transverse component must be positive (prograde)');
assert(abs(U04(1,1)) < 1e-9, ...
    'radial component at apogee start node (L=pi, e=0.75) must be ~0: got %.3e', U04(1,1));

assert(all(diff(X04(1,:)) > 0), 'P must strictly increase under tangential thrust');
mdrop4_expected = (p.Tmax/p.c) * opts4.thr * inf4.tEnd;
mdrop4_actual   = 1 - X04(6,end);
assert(abs(mdrop4_actual - mdrop4_expected) < 1e-6*mdrop4_expected, ...
    'tangential mass drop inconsistent with mdot = -(Tm/c)*thr: got %.6e vs %.6e', ...
    mdrop4_actual, mdrop4_expected);

fprintf('test_mee_seed: ALL PASS (nRev1=%.4f mEnd1=%.6f, stopP nRev=%.4f, ballistic tEnd=%.4f vs 2T=%.4f, tangential nRev=%.4f mEnd=%.6f)\n', ...
    inf1.nRev, inf1.mEnd, inf2.nRev, inf3.tEnd, 2*Tperiod, inf4.nRev, inf4.mEnd);

function dmax = seed_stencil_defect(sg, X0, U0, dL0, p)
% Trapezoid defect of the seed under the solver's sigma-stencil (numeric
% mirror of test_seed.m's seed_stencil_defect, adapted to the L-domain MEE
% RHS: dXdsigma = dL0 * dXdL).
N = numel(sg) - 1;
f = zeros(7, N+1);
for k = 1:N+1
    Lk = pi + sg(k)*dL0;
    par_k = p;  par_k.L = Lk;
    dXdL = lt_mee_rhs(X0(:,k), U0(:,k), par_k);
    f(:,k) = dL0 * dXdL;
end
dmax = 0;
for k = 1:N
    d = X0(:,k+1) - X0(:,k) - ((sg(k+1)-sg(k))/2)*(f(:,k)+f(:,k+1));
    dmax = max(dmax, max(abs(d)));
end
end
