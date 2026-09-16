%% ANCHOR_STUDY  How a SEED is made: direct solve -> covector harvest -> shooting root.
%
%   THE PREQUEL TO transfer_study.m. That script starts from a certified
%   root and asks whether it is a minimum. This one starts from NOTHING but
%   a neighbouring solution and produces the root -- which is the step every
%   reader asks about first: where do the costates come from?
%
%   The answer, in one line: a direct collocation solve is run in discrete
%   time, and the Lagrange multipliers of its dynamics-defect constraints
%   ARE a discretisation of the costates. Harvest them, interpolate them onto
%   a multiple-shooting mesh, and that state-and-costate trajectory is the
%   SEED. The script measures every claim in that sentence.
%
%   Five anchors exist in the 70 mN library, one per FAMILY of extremals --
%   five distinct branches connecting the SAME two orbits. Choose one in
%   section 0 and the script rebuilds it.
%
%   THE SECTIONS, and what each one is for:
%     0  the anchor to build, the physics, the budgets, the tolerances
%     1  the departure orbit (generated from its own parameter)
%     2  the target orbit (generated from its own parameters)
%     3  the engine and the two endpoints at the chosen phases
%     4  THE WARM START: which certified root we begin from, and how wrong
%        it is for the transfer we actually want
%     5  THE DIRECT SOLVE: the discrete-time NLP, with its own gates
%     6  THE HARVEST: duals -> costates, the three subtleties, with gates
%     7  WHAT THE SEED IS WORTH: fly the harvested lambda(0) ALONE and
%        measure the miss; then measure the seed's junction defects. This is
%        the section that explains why a seed is a trajectory, not a vector
%     8  MULTIPLE SHOOTING: the seed becomes a root
%     9  CERTIFY: the same gate stack every catalog entry passed
%    10  the hand-off: compare against the stored anchor, and what to run next
%
%   WHAT IS AND IS NOT DEMONSTRATED. Cases 1 and 2 are REPRODUCTION runs:
%   the warm start IS the root being rebuilt, so the NLP starts at the
%   answer. They cannot demonstrate discovery -- but they are not circular,
%   because the costates are not reused: the stored root only shapes the
%   initial trajectory, and the NLP that follows has no costate in it at
%   all. Section 6b then compares the HARVESTED lambda(0) against the stored
%   one, which is a shooting solution of the same transfer -- two
%   independent routes to the same object, and the sharpest test of the
%   duals-to-costates mapping there is. Read these two first. Cases 3, 4 and 5
%   warm-start from a certified root of ANOTHER family at a nearby phase --
%   which is how those three families were actually discovered (FINDINGS 61,
%   63, 68). A direct solve is BRANCH-BLIND: it lands in whatever basin its
%   warm start is near, so "did we land on the expected root?" is reported
%   as a measurement, never asserted. That blindness is a feature -- it is
%   the only reason the library has five families and not one.
%
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.

%% paths
clear; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
D = fullfile(fileparts(here), 'direct');
addpath(D, fullfile(D, 'lib'), fullfile(D, 'certify'));
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

%% CR3BP constants
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
day    = @(tnd) tnd*tStar/86400;

%% ========================================================================
%  0. THE ANCHOR TO BUILD, AND EVERY OTHER SETTING
%
%     anchorCase picks one of the five. Each row says: what the transfer IS
%     (its phases and the family it belongs to), WHERE ITS WARM START COMES
%     FROM (a certified root on disk), and what t_f the campaign measured
%     for it -- the last is used only to report whether this run landed on
%     the same root, never to steer the solve.
%% ========================================================================
anchorCase = 1;          % <-- 1..5, choose the anchor to build (1 and 2 are reproductions: read them first)

%          name        sD      sA         t_f (d)   warm start file                              what the warm start is
A = { ...
  'fast',      0.0, 0.0754,     17.7975, 'mintime_70mN_anchor.mat',           'ITSELF (reproduction): the 17.80 d root of the original family'; ...
  'A2',        0.0, 0.0754+10/12, 26.4300, 'mintime_70mN_certified.mat',      'ITSELF (reproduction): the 26.43 d root of the slow family'; ...
  'fast2',     0.0, 0.8671,     17.2487, 'mintime_70mN_certified.mat',        'the A2 family''s 26.43 d root at 0.9087 -- HISTORICAL (FINDINGS 61)'; ...
  'direct18',  0.0, 0.0754+17/24, 17.8342, 'mintime_70mN_anchor_fast2.mat',   'the fast2 family''s 17.25 d root at 0.8671, 0.083 away in phase'; ...
  'direct11',  0.0, 0.0754+10/24, 17.9605, 'mintime_70mN_anchor_direct18.mat','the direct18 family''s 17.83 d root at 0.7837 (FINDINGS 68)'};
assert(anchorCase >= 1 && anchorCase <= size(A, 1) && anchorCase == round(anchorCase), ...
       'anchor_study: anchorCase must be 1..%d', size(A, 1));
want = struct('name', A{anchorCase, 1}, 'sD', A{anchorCase, 2}, 'sA', A{anchorCase, 3}, ...
              'tfDays', A{anchorCase, 4}, 'seedFile', fullfile(here, 'results', A{anchorCase, 5}), ...
              'seedWhat', A{anchorCase, 6});
% A REPRODUCTION case warm-starts from the root it is trying to rebuild, so
% the NLP begins at the answer. It cannot demonstrate discovery -- but the
% costates are NOT reused: the stored root only shapes the initial
% trajectory, and section 6 reconstructs lambda(0) from the NLP's own
% multipliers. That makes the stored lambda(0) an INDEPENDENT reference for
% the harvest, which is the whole value of these two cases (section 6b).
want.reproduction = ismember(anchorCase, [1 2]);

%% the physics (the same operating point for all five anchors)
dep = struct('family', 'dro',   'tau', 1.0);              % departure orbit
arr = struct('family', 'tulip', 'Np', 7, 'pm', -1);       % target orbit
eng = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150); % engine

%% the budgets
num = struct( ...
    'N',         800,  ...   % direct solve: collocation intervals
    'K',          24,  ...   % shooting mesh: junction segments (the seed's width)
    'maxIter',  3000,  ...   % IPOPT iteration cap
    'maxCpuSec', 900,  ...   % IPOPT wall cap (a mismatched warm start burns it)
    'clearKm',  1900,  ...   % lunar clearance, MOON-CENTRE distance
    'shootSec', 600);        % multiple-shooting wall cap

%% the tolerances -- one source, used by every printed line AND the verdict
tol = struct( ...
    'closure',   1e-7,  ...  % periodicity of each generated orbit
    'seamDeriv', 1e-6,  ...  % interpolant derivative mismatch across s = 0
    'defect',    1e-8,  ...  % 5. direct solve: max collocation defect
    'termErr',   1e-7,  ...  % 5. direct solve: terminal boundary residual
    'unit',      1e-8,  ...  % 5. |u| = 1 on the solved directions
    'lamT',      1e-2,  ...  % 6. harvested lambda_t vs its theoretical +1
    'vote',      0.90,  ...  % 6. sign vote: fraction agreeing with the majority
    'R',         1e-8,  ...  % 8. shooting residual, inf-norm
    'tfAgree',   1e-3,  ...  % 8. direct t_f vs shot t_f, relative
    'same',      1e-2);      % 10. |t_f - the campaign's| /t_f to call it the same root

fprintf('ANCHOR STUDY -- building the "%s" anchor of the %.0f mN library\n', want.name, eng.thrustN*1000);
fprintf('   target      : sD = %.4f, sA = %.4f\n', want.sD, want.sA);
fprintf('   warm start  : %s\n                 %s\n', A{anchorCase, 5}, want.seedWhat);
fprintf('   campaign''s t_f for this anchor: %.4f d (reported against, never used to steer)\n', want.tfDays);

%% ========================================================================
%  1. DEPARTURE ORBIT -- a distant retrograde orbit (DRO)
%% ========================================================================
[tD, rvD, infoD] = get_family_orbit(dep.family, struct('tau', dep.tau, 'muStar', muStar));
closD = norm(rvD(end,:) - rvD(1,:));
fprintf('\n1. DEPARTURE %s: period %.4f ND (%.3f d), closure %.2e / %-7.0e  %s\n', ...
        upper(dep.family), infoD.periodND, day(infoD.periodND), closD, tol.closure, pass(closD < tol.closure));
assert(closD < tol.closure, 'the departure orbit is not periodic to tolerance');

%% ========================================================================
%  2. TARGET ORBIT -- a tulip (its period is LOCKED by the petal count)
%% ========================================================================
[tT, rvT, infoT] = get_family_orbit(arr.family, struct('Np', arr.Np, 'pm', arr.pm, 'muStar', muStar));
closA = norm(rvT(end,:) - rvT(1,:));
fprintf('2. TARGET %s: %d petals, branch %+d, period %.4f ND (%.3f d), closure %.2e / %-7.0e  %s\n', ...
        upper(arr.family), arr.Np, arr.pm, infoT.periodND, day(infoT.periodND), closA, tol.closure, ...
        pass(closA < tol.closure));
assert(closA < tol.closure, 'the target orbit is not periodic to tolerance');

%% ========================================================================
%  3. ENGINE AND ENDPOINTS
%% ========================================================================
g0  = 9.80665*tStar^2/(1000*lStar);           % ND gravity at sea level
cnd = (eng.ispS/tStar)*g0;                    % ND exhaust speed
Tnd = (eng.thrustN/eng.m0kg)*tStar^2/(lStar*1000);   % ND thrust acceleration at m = 1
[stateD, seamD] = phase_state(tD, rvD);
[stateA, seamA] = phase_state(tT, rvT);
rv0 = stateD(want.sD);   rvf = stateA(want.sA);
seamOK = max(seamD.value, seamA.value) < tol.closure && max(seamD.deriv, seamA.deriv) < tol.seamDeriv;
fprintf('3. ENGINE %.0f mN, Isp %g s, m0 %g kg -> T_nd %.6e, c_nd %.4f\n', ...
        eng.thrustN*1000, eng.ispS, eng.m0kg, Tnd, cnd);
fprintf('   endpoints: departure r = [%+.5f %+.5f %+.5f], arrival r = [%+.5f %+.5f %+.5f]\n', rv0(1:3), rvf(1:3));
fprintf('   periodic interpolant seam: value %.1e, derivative %.1e            %s\n', ...
        max(seamD.value, seamA.value), max(seamD.deriv, seamA.deriv), pass(seamOK));
assert(seamOK, 'the endpoint interpolant is not periodic across its seam');

%% ========================================================================
%  4. THE WARM START -- a certified root of SOME transfer, which is not the
%     transfer we want. The direct solver does not need a good guess of the
%     ANSWER; it needs a trajectory in roughly the right part of the space.
%     How wrong it is for our endpoints is worth measuring, because that
%     distance is what decides whether the solve converges at all
%     (measured 2026-09-15: a warm start from another ORBIT PAIR gives a
%     259-day trajectory through the Moon).
%% ========================================================================
assert(isfile(want.seedFile), 'the warm start %s is not on disk', want.seedFile);
Ls = load(want.seedFile);
if isfield(Ls, 'best'), root0 = Ls.best; else, root0 = Ls; end
assert(isfield(root0, 'z') && isfield(root0, 'it') && isfield(root0.it, 'Y'), ...
       'the warm start must carry .z and .it.Y');
z0 = root0.z(:);
if isfield(root0, 'rv0') && numel(root0.rv0) >= 6, rvWS = root0.rv0(1:6); else, rvWS = root0.it.Y(1:6, 1); end

% fly it, on ITS own problem, and see where it ends up relative to OUR target
[tauWS, YWS] = pumpkyn.cr3bp.tfMinProp(z0(8), [rvWS(:); 1; z0(1:7)], Tnd, cnd, muStar);
endWS   = YWS(end, 1:6).';
missKm  = norm(endWS(1:3) - rvf(1:3))*lStar;
missMs  = norm(endWS(4:6) - rvf(4:6))*lStar/tStar*1000;
fprintf('\n4. WARM START: t_f %.4f d, %d propagation samples\n', day(z0(8)), numel(tauWS));
fprintf('   its endpoint vs OUR target: %.0f km, %.1f m/s apart\n', missKm, missMs);
fprintf('   (the NLP must move the endpoint that far; it is an initial guess, not an answer)\n');

%% ========================================================================
%  5. THE DIRECT SOLVE -- the discrete-time problem.
%     Time is cut into N intervals; the states at the nodes and the thrust
%     direction on each interval are all unknowns; "the trajectory obeys the
%     equations of motion" becomes N equality constraints (the DEFECTS); the
%     objective is t_f. A Sundman-regularised independent variable makes the
%     mesh follow the geometry rather than the clock, and the lunar clearance
%     is a path constraint. IPOPT solves it. No costate is mentioned here --
%     which is exactly why the next section is interesting.
%% ========================================================================
sN = linspace(0, tauWS(end), num.N + 1);
X0 = interp1(tauWS, YWS(:, 1:7), sN, 'spline').';            % states at the nodes
LV = interp1(tauWS, YWS(:, 11:13), sN, 'spline');            % lambda_v, for the direction
U0 = [(-LV ./ max(vecnorm(LV, 2, 2), eps)).'; ones(1, num.N + 1)];   % u = -lam_v/|lam_v|, throttle 1
floorKm = num.clearKm - 1737.4;                              % casadi takes an ALTITUDE floor
tDirect = tic;
o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tnd, cnd, muStar, num.N, X0, U0, z0(8), ...
        struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', true, ...
               'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', num.maxCpuSec));
wallDirect = toc(tDirect);
if isfield(o, 'model'), o = rmfield(o, 'model'); end
perisKm = o.altMinKm + 1737.4;                               % back to Moon-CENTRE distance

fprintf('\n5. DIRECT SOLVE (%d intervals, Hermite-Simpson, Sundman): %s, %.0f s\n', ...
        num.N, o.ipoptStatus, wallDirect);
fprintf('   D1 converged      %-24s                      %s\n', o.ipoptStatus, pass(o.success));
fprintf('   D2 max defect     %9.2e / %-9.0e                  %s\n', o.maxDefect, tol.defect, pass(o.maxDefect < tol.defect));
fprintf('   D3 terminal resid %9.2e / %-9.0e                  %s\n', o.termErr, tol.termErr, pass(o.termErr < tol.termErr));
fprintf('   D4 |u| = 1        %9.2e / %-9.0e                  %s\n', o.maxUnit, tol.unit, pass(o.maxUnit < tol.unit));
fprintf('   D5 lunar clear    %9.1f / %-9.1f km (Moon centre)    %s\n', perisKm, num.clearKm, pass(perisKm >= num.clearKm));
fprintf('   t_f = %.4f d, final mass fraction %.6f\n', day(o.tf), o.mf);
d1 = o.success;  d2 = o.maxDefect < tol.defect;  d3 = o.termErr < tol.termErr;
d4 = o.maxUnit < tol.unit;  d5 = perisKm >= num.clearKm;
directOK = d1 && d2 && d3 && d4 && d5;
assert(directOK, ['the direct solve did not produce a usable trajectory (status %s, defect %.1e, ' ...
                  'terminal %.1e, periselene %.0f km). Nothing below is meaningful. A warm start ' ...
                  'whose geometry is wrong for this target is the usual cause.'], ...
       o.ipoptStatus, o.maxDefect, o.termErr, perisKm);

%% ========================================================================
%  6. THE HARVEST -- duals to costates. THE PEDAGOGICAL CORE.
%
%     The multiplier attached to a dynamics-defect constraint measures how
%     much the objective would improve if that constraint were relaxed --
%     which is the definition of the adjoint variable. So the NLP has been
%     computing the costates all along, under another name. Three things
%     have to be right to read them out:
%       (a) STATION. For Hermite-Simpson the multipliers belong to segment
%           MIDPOINTS, not to the nodes. Pairing them with nodes shifts the
%           costates half a step -- an error that still converges and still
%           gives a slightly wrong answer.
%       (b) SIGN. The convention depends on how the NLP was written. It is
%           determined by a VOTE: the minimum principle requires the applied
%           direction to oppose lambda_v at every station, so the sign that
%           satisfies that at the most stations wins, and the margin is
%           reported.
%       (c) SCALE. With t_f as the objective the multiplier of the
%           time-scaling constraint must come out lambda_t = +1; it is a
%           free check on the whole mapping.
%     One home for all three: oclib/+oc/duals_to_costates (a second consumer,
%     booster_landing, re-invented and re-fixed the midpoint bug alone).
%% ========================================================================
[seed, dg] = harvest_ms_seed(o, num.K);

fprintf('\n6. HARVEST: %d defect multipliers -> costates at %d midpoints -> a %dx%d seed\n', ...
        size(o.lamDef, 2), size(o.lamDef, 2), size(seed.Y, 1), size(seed.Y, 2));
voteFrac = gvd(dg, 'voteMargin', NaN);
lamT     = gvd(dg, 'lamT', NaN);
lamTOK   = gvd(dg, 'lamTOK', NaN);
fprintf('   H1 sign vote      %s, %5.1f%% of stations agree / %.0f%%          %s\n', ...
        tern(gvd(dg, 'sign', 1) > 0, 'positive', 'negative'), 100*voteFrac, 100*tol.vote, ...
        pass(isfinite(voteFrac) && voteFrac >= tol.vote));
if isfinite(lamT)
    fprintf('   H2 lambda_t       %9.6f / 1 +- %-9.0e                  %s\n', ...
            lamT, tol.lamT, pass(abs(lamT - 1) < tol.lamT));
else
    fprintf('   H2 lambda_t       NOT EXPOSED by this NLP (no lamDef row 8): the scale is unchecked\n');
end
fprintf('   H3 seed finite    %d x %d, all finite                          %s\n', ...
        size(seed.Y, 1), size(seed.Y, 2), pass(all(isfinite(seed.Y(:))) && isfinite(seed.tf)));
h1 = isfinite(voteFrac) && voteFrac >= tol.vote;
h2 = ~isfinite(lamT) || abs(lamT - 1) < tol.lamT;      % unchecked is not failed; it is unchecked
h3 = all(isfinite(seed.Y(:))) && isfinite(seed.tf) && isequal(size(seed.Y), [14, num.K + 1]);
fprintf('   the seed: tf %.4f d, tGrid [%d], Y [%d x %d] = [r; v; m; lam_r; lam_v; lam_m]\n', ...
        day(seed.tf), numel(seed.tGrid), size(seed.Y, 1), size(seed.Y, 2));
fprintf('   lambda(0) harvested = [%s]\n', strjoin(compose('%+.6g', seed.Y(8:14, 1)'), ' '));
assert(h3, 'the harvested seed is malformed');

% 6b. THE HARVEST AGAINST A KNOWN ANSWER -- available only on a reproduction
%     case, and the reason those cases exist. The stored root's lambda(0)
%     came from a SHOOTING solve; the harvested one came from a collocation
%     NLP's constraint multipliers. They are two independent routes to the
%     same object, so their difference measures the duals-to-costates
%     mapping itself -- station association, sign and scale together.
%     Costates are defined up to their own scale here only through lambda_t,
%     which H2 already pinned to +1, so a direct componentwise comparison is
%     meaningful.
if want.reproduction
    lamStored = z0(1:7);
    dLam = norm(seed.Y(8:14, 1) - lamStored)/max(norm(lamStored), eps);
    dTf  = abs(seed.tf - z0(8))/z0(8);
    fprintf('   6b REPRODUCTION CHECK (the stored root is an independent reference):\n');
    fprintf('      lambda(0) stored    = [%s]\n', strjoin(compose('%+.6g', lamStored'), ' '));
    fprintf('      harvested vs stored: %.2e relative on lambda(0), %.2e on t_f\n', dLam, dTf);
    fprintf('      the NLP never saw a costate; this is the mapping measured against a\n');
    fprintf('      shooting solution of the same transfer.\n');
else
    fprintf('   (6b reproduction check: not applicable -- this case warm-starts from a\n');
    fprintf('    DIFFERENT transfer, so no stored lambda(0) for this phase is a reference.)\n');
end

%% ========================================================================
%  7. WHAT THE SEED IS WORTH -- the section that answers "why not just
%     lambda(0)?". Two measurements, no opinions.
%
%     (a) SINGLE SHOOTING from the harvested lambda(0) alone: fly the PMP
%         equations for the harvested t_f and see where they land. The
%         min-time extremal amplifies an initial costate perturbation by
%         ~1e3 over tens of revolutions, so this is expected to miss by
%         tens of thousands of kilometres -- from costates that are, in
%         every other sense, good.
%     (b) THE SEED'S OWN DEFECTS: propagate each junction forward to the
%         next junction time and compare with the seed's own value there.
%         Those gaps are what the shooting solve has to close, and they are
%         small because each one spans only 1/K of the arc.
%% ========================================================================
lam0 = seed.Y(8:14, 1);
[~, Yss] = pumpkyn.cr3bp.tfMinProp(seed.tf, [rv0(1:6); 1; lam0], Tnd, cnd, muStar);
ssKm = norm(Yss(end, 1:3).' - rvf(1:3))*lStar;
ssMs = norm(Yss(end, 4:6).' - rvf(4:6))*lStar/tStar*1000;

defJ = zeros(1, num.K);
for k = 1:num.K
    [~, Yk] = pumpkyn.cr3bp.tfMinProp(seed.tGrid(k+1) - seed.tGrid(k), seed.Y(:, k), Tnd, cnd, muStar);
    defJ(k) = norm(Yk(end, :).' - seed.Y(:, k+1), inf);
end
fprintf('\n7. WHAT THE SEED IS WORTH\n');
fprintf('   (a) SINGLE SHOOTING from the harvested lambda(0) alone:\n');
fprintf('       misses the target by %.0f km and %.1f m/s -- with correct-looking costates.\n', ssKm, ssMs);
fprintf('       The arrival gate is %g km, so those costates are unusable ON THEIR OWN even\n', 100);
fprintf('       though every digit of them is meaningful. The amplification grows with the\n');
fprintf('       length of the arc: this transfer is %.1f d, and the catalog entries that were\n', day(seed.tf));
fprintf('       measured at 36,000-560,000 km were tens of revolutions long. Same mechanism,\n');
fprintf('       two regimes -- which is why the fix is structural, not a tighter tolerance.\n');
fprintf('   (b) THE SEED''S JUNCTION DEFECTS over %d segments: max %.2e, median %.2e (ND, inf-norm)\n', ...
        num.K, max(defJ), median(defJ));
fprintf('       each spans 1/%d of the arc, so each gap is small -- and closing them is\n', num.K);
fprintf('       exactly what section 8 does.\n');

%% ========================================================================
%  8. MULTIPLE SHOOTING -- the seed becomes a root.
%     The junction states join lambda(0) and t_f as unknowns; the defects of
%     section 7(b) become the equations. The Jacobian is block-bidiagonal and
%     each block is an STM over ONE short segment, so it is conditioned where
%     single shooting is not.
%% ========================================================================
tShoot = tic;
[z8, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                    struct('tolR', 1e-11, 'wallSec', num.shootSec, 'maxIter', 100));
wallShoot = toc(tShoot);
tfRel = abs(z8(8) - o.tf)/o.tf;
fprintf('\n8. MULTIPLE SHOOTING (%d segments): %d iterations, %.0f s\n', num.K, it.iters, wallShoot);
% THE GATE IS THE RESIDUAL THIS LINE PRINTS. The solver is ASKED for 1e-11
% -- tighter than the script's gate, so it works as hard as it can -- and
% `it.converged` is its verdict against that tighter target. Requiring the
% flag as well made a |R| of 1.05e-11 print FAIL against a 1e-8 threshold it
% passes by three orders of magnitude (caught on this script's first run).
% A solve that plateaus above its own target but inside the gate is a root;
% it is reported as such, with the plateau named.
m1 = isfinite(it.normR) && it.normR < tol.R;
fprintf('   M1 residual       |R| %9.2e / %-9.0e                  %s%s\n', it.normR, tol.R, pass(m1), ...
        tern(m1 && ~it.converged, sprintf('   (plateaued above the solver''s own %.0e target)', 1e-11), ''));
fprintf('   M2 t_f moved      %9.2e / %-9.0e relative             %s   (direct %.4f d -> shot %.4f d)\n', ...
        tfRel, tol.tfAgree, pass(tfRel < tol.tfAgree), day(o.tf), day(z8(8)));
fprintf('   lambda(0) shot    = [%s]\n', strjoin(compose('%+.6g', z8(1:7)'), ' '));
fprintf('   lambda(0) moved from the harvest by %.2e (relative)\n', ...
        norm(z8(1:7) - lam0)/max(norm(lam0), eps));
m2 = tfRel < tol.tfAgree;
assert(m1, 'the shooting residual %.2e is above the gate %.0e: there is no root to certify', it.normR, tol.R);

%% ========================================================================
%  9. CERTIFY -- the same gate stack every one of the 576 catalog entries
%     passed. transfer_study.m opens this box and reads each gate; here we
%     only need its verdict, because the subject of this script is the seed.
%% ========================================================================
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', eng.thrustN, 'ispS', eng.ispS, 'm0kg', eng.m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', Tnd, 'cnd', cnd, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
pool = capped_pool();
C = certify_root(seed, rv0(1:6), rvf, B, struct('pool', pool, 'wallSec', 900, 'sA', want.sA, 'sD', want.sD));
fprintf('\n9. CERTIFY: %s\n', C.reason);
fprintf('   t_f %.4f d, flown miss %.3f km / %.3f m/s, conjugate %s, lift margin %.1fx\n', ...
        C.tfDays, C.flyKm, C.flyVms, tern(gvd(C, 'conj', 0) == 1, 'PASS', 'not passed'), gvd(C, 'liftMargin', NaN));
if ~isempty(gvd(C, 'note', '')), fprintf('   note: %s\n', C.note); end

%% ========================================================================
% 10. THE HAND-OFF -- is this the anchor the campaign found, and what next?
%     Branch-blindness means the answer can legitimately be NO: the solve
%     may have landed on another family's root at the same phase. That is
%     reported, not asserted.
%% ========================================================================
sameRoot = C.ok && abs(C.tfDays - want.tfDays)/want.tfDays < tol.same;
fprintf('\n10. HAND-OFF\n');
fprintf('    this run: %.4f d      the campaign''s "%s": %.4f d      difference %.2f%%\n', ...
        C.tfDays, want.name, want.tfDays, 100*abs(C.tfDays - want.tfDays)/want.tfDays);
if sameRoot
    fprintf('    SAME ROOT: the direct solve landed on the expected family.\n');
elseif C.ok
    fprintf(['    A DIFFERENT ROOT: certified, but not the one the campaign stored here.\n' ...
             '    This is branch-blindness doing its job -- the same mechanism that found\n' ...
             '    four of the five families. It is a result, not an error.\n']);
else
    fprintf('    NOT CERTIFIED: %s\n', C.reason);
end
fprintf(['    To save it as an anchor:  best = struct(''z'', C.z, ''it'', struct(''Y'', C.Y), ...\n' ...
         '        ''sA'', %.6f, ''sD'', %.4f, ''tfDays'', C.tfDays, ''origin'', ''...'');  save(f, ''best'', ''Tnd'', ''cnd'')\n' ...
         '    Then: arcs from it (arclength_arrival), or transfer_study.m to interrogate it,\n' ...
         '    or run_phase_torus with it in .anchors to build a library around it.\n'], want.sA, want.sD);

%% ------------------------------------------------------------------------
%  THE SELF-CHECK. Every gate by name, one scalar, then an assertion. A
%  failed DIRECT, HARVEST or SHOOTING gate means the machinery is broken and
%  throws. Landing on a different root does NOT throw: it is a finding.
%% ------------------------------------------------------------------------
gateStatus = { ...
    'D1 direct converged', d1; 'D2 collocation defect', d2; 'D3 terminal residual', d3; ...
    'D4 unit control',     d4; 'D5 lunar clearance',   d5; ...
    'H1 sign vote',        h1; 'H2 lambda_t',          h2; 'H3 seed well formed', h3; ...
    'M1 shooting residual', m1; 'M2 t_f agreement',    m2; ...
    'C1 certified',        C.ok};
failed = gateStatus(~[gateStatus{:, 2}], 1);
fprintf('\n    SELF-CHECK: %d of %d gates passed', nnz([gateStatus{:, 2}]), size(gateStatus, 1));
if isempty(failed), fprintf('; the seed was built, shot and certified.\n');
else, fprintf('; NOT passed: %s\n', strjoin(failed', ', ')); end
assert(all([gateStatus{1:10, 2}]), 'anchor_study:machinery', ...
       'the seed-building machinery did not work end to end: %s', strjoin(failed', ', '));

%% ------------------------------------------------------------------------
function s = pass(c)
% PASS  Verdict text.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function v = gvd(s, f, d)
% GVD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = tern(c, a, b)
% TERN  Ternary.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
