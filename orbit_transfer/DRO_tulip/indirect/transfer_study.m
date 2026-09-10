%% TRANSFER_STUDY  One minimum-time low-thrust transfer, with the scaffolding showing.
%
%   Edit the parameter blocks, press Run. Every step is done HERE rather than
%   behind a front door, so the machinery is visible:
%
%     1  generate the DEPARTURE orbit      (family + parameters -> a real orbit)
%     2  generate the TARGET orbit         (same)
%     3  engine, phases, endpoints         (nondimensionalisation spelled out)
%     4  solve                             (seed -> multiple shooting -> costates)
%     5  independent verification          (a second solver must not move them)
%     6  NECESSARY conditions              (Pontryagin, one at a time)
%     7  SUFFICIENCY hypotheses            (Bonnard-Caillau-Trelat, one at a time)
%     8  interactive 3D plot               (drag to rotate)
%
%   Sections 6 and 7 are separate because the theory separates them: the
%   first-order conditions are HYPOTHESES of the sufficiency theorem, not
%   consequences of it, and the conjugate test is computed ALONG the extremal
%   -- off one, its determinant means nothing. On this problem 12 of 14
%   candidates satisfied every first-order condition and were then refuted.
%
%   Section 5 is not an optimality condition at all. It guards against OUR
%   solver: a bug in our shooting could give a self-consistent answer to the
%   wrong problem, and only a second implementation catches that.
%
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.

clear; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
day    = @(tnd) tnd*tStar/86400;

%% ========================================================================
%  1. DEPARTURE ORBIT -- a distant retrograde orbit (DRO)
%     Parameter: the PERIOD tau (ND). Larger tau = larger orbit.
%       tau 0.5 -> 2.22 d      tau 1.0 -> 4.43 d      tau 2.0 -> 8.87 d
%% ========================================================================
dep.family = 'dro';
dep.tau    = 1.0;                            % <-- change the departure orbit here

[tD, rvD, infoD] = get_family_orbit(dep.family, struct('tau', dep.tau, 'muStar', muStar));
closureD = norm(rvD(end,:) - rvD(1,:));
fprintf('1. DEPARTURE %s: period %.4f ND (%.3f d), %d samples\n', ...
        upper(dep.family), infoD.periodND, day(infoD.periodND), numel(tD));
fprintf('   max radius from the barycentre %.4f ND (%.0f km)\n', ...
        max(vecnorm(rvD(:,1:3),2,2)), max(vecnorm(rvD(:,1:3),2,2))*lStar);
fprintf('   periodicity closure |x(T) - x(0)| = %.2e   %s\n', closureD, ...
        pass(closureD < 1e-8));

%% ========================================================================
%  2. TARGET ORBIT -- a tulip
%     Parameters: the PETAL COUNT Np and the branch pm (+1 is the z-mirror
%     of -1). The period is NOT free: it is locked to 2*pi*(Np-2)/(Np-1).
%       Np 3 -> 13.93 d    Np 5 -> 20.89 d    Np 7 -> 23.21 d    Np 9 -> 24.37 d
%% ========================================================================
arr.family = 'tulip';
arr.Np     = 7;                              % <-- change the target orbit here
arr.pm     = -1;                             %     branch: -1 or +1

[tT, rvT, infoT] = get_family_orbit(arr.family, struct('Np', arr.Np, 'pm', arr.pm, 'muStar', muStar));
closureT = norm(rvT(end,:) - rvT(1,:));
fprintf('\n2. TARGET %s: %d petals, branch %+d, period %.4f ND (%.3f d) LOCKED by Np\n', ...
        upper(arr.family), arr.Np, arr.pm, infoT.periodND, day(infoT.periodND));
fprintf('   out-of-plane extent max|z| = %.4f ND (%.0f km)\n', ...
        max(abs(rvT(:,3))), max(abs(rvT(:,3)))*lStar);
fprintf('   periodicity closure |x(T) - x(0)| = %.2e   %s\n', closureT, ...
        pass(closureT < 1e-7));

%% ========================================================================
%  3. ENGINE, PHASES, ENDPOINTS
%     A phase is a fraction of its orbit's period. The target state is read
%     from a PERIODIC SPLINE through the propagated orbit -- that interpolant,
%     not the raw table and not the dynamics, is what the boundary condition
%     matches.
%% ========================================================================
thrustN = 0.070;      ispS = 900;      m0kg = 150;
sD      = 0.0;                                % departure phase, fraction of tau_D
sA      = 0.0754;                             % arrival phase,  fraction of tau_A

g0  = 9.80665*tStar^2/(1000*lStar);           % ND gravity at sea level
cnd = (ispS/tStar)*g0;                        % ND exhaust speed
Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);    % ND thrust acceleration at m = 1

ppD = spline(tD(:).', rvD.');   ppA = spline(tT(:).', rvT.');
stateD = @(s) ppval(ppD, mod(s,1)*tD(end));
stateA = @(s) ppval(ppA, mod(s,1)*tT(end));
rv0 = stateD(sD);   rvf = stateA(sA);

fprintf('\n3. ENGINE %.0f mN, Isp %g s, m0 %g kg  ->  T_nd = %.6e, c_nd = %.4f\n', ...
        thrustN*1000, ispS, m0kg, Tnd, cnd);
fprintf('   departure phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sD, rv0(1:3));
fprintf('   arrival   phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sA, rvf(1:3));

%% ========================================================================
%  4. SOLVE -- seed, then multiple shooting on the PMP boundary-value problem
%     Unknowns: the seven initial costates and the final time. The seed comes
%     from the certified library when this operating point is in it; a new
%     orbit pair needs a continuation from a solved one (run_dro_tulip).
%% ========================================================================
lib = dro_tulip_library();
match = find(abs([lib.sD] - mod(sD,1)) < 1e-6 & abs([lib.sA] - mod(sA,1)) < 1e-6, 1);
assert(~isempty(match) && abs(dep.tau - 1) < 1e-12 && arr.Np == 7, ...
    ['no seed for this operating point. The library covers the tau = 1 DRO to the\n' ...
     '7-petal tulip at 70 mN. For another orbit pair, walk to it with\n' ...
     '  T = run_dro_tulip(sD, sA, opts)   and then study T here.']);

K = size(lib(match).Y, 2);
seed = struct('tf', lib(match).z(8), 'tGrid', linspace(0, lib(match).z(8), K+1), ...
              'Y', [lib(match).Y, lib(match).Y(:,end)]);
seed.Y(1:7,1) = [rv0(1:6); 1];   seed.Y(8:14,1) = lib(match).z(1:7);

[z8, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                    struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));

[tu, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], Tnd, cnd, muStar);
mf   = Y(end,7);
dV   = cnd*log(1/mf)*lStar/tStar;
fprintf('\n4. SOLVED: %d Newton iterations, |R| = %.2e, converged = %d\n', ...
        it.iters, it.normR, it.converged);
fprintf('   t_f = %.6f ND (%.4f d)   Delta-V = %.4f km/s   propellant %.2f kg\n', ...
        z8(8), day(z8(8)), dV, m0kg*(1-mf));
fprintf('   lambda_0 = [%s]\n', strjoin(compose('%+.6g', z8(1:7)'), ' '));

%% ========================================================================
%  5. INDEPENDENT VERIFICATION -- a second solver must not move the costates
%% ========================================================================
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', Tnd, 'cnd', cnd, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
V = verify_with_pumpkyn(struct('z', z8, 'sD', sD, 'sA', sA), B);

%% ========================================================================
%  6. NECESSARY CONDITIONS (Pontryagin, first order) -- one at a time
%% ========================================================================
fprintf('\n6. NECESSARY CONDITIONS\n');

% N1  the boundary-value residual: costate equations, terminal matching,
%     transversality. This IS the statement that the first variation vanishes.
[~, chk] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, struct('assembleOnly', true));
R1 = chk.residual([z8(1:7); reshape(it.Y(:,2:end), [], 1); z8(8)]);
fprintf('   N1 BVP residual  |R|_inf = %.2e                       %s\n', ...
        norm(R1, inf), pass(norm(R1, inf) < 1e-8));

% N2  the Hamiltonian. Autonomous problem, free final time  =>  H == 0.
%     H = 1 + lam_r.v + lam_v.(g + (T/m) alpha) - lam_m T/c,  alpha = -lam_v/|lam_v|
lamr = Y(:, 8:10);   lamv = Y(:, 11:13);   lamm = Y(:, 14);   mm = Y(:, 7);
nlv  = vecnorm(lamv, 2, 2);
Hval = zeros(size(tu));
for k = 1:numel(tu)
    F = mintime_rhs_point(Y(k,:).', Tnd, cnd, muStar);     % [xdot; lamdot]
    Hval(k) = 1 + Y(k,8:14)*F(1:7);
end
fprintf('   N2 Hamiltonian   max|H| along the arc = %.2e         %s\n', ...
        max(abs(Hval)), pass(max(abs(Hval)) < 1e-6));

% N3  the flight actually reaches the target, in position AND velocity
missKm  = norm(Y(end,1:3) - rvf(1:3)')*lStar;
missVms = norm(Y(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
fprintf('   N3 arrival       %.4f km, %.4f m/s                   %s\n', ...
        missKm, missVms, pass(missKm < 100 && missVms < 10));

% N4  transversality on the free mass: lam_m(t_f) = 0
fprintf('   N4 transversality lam_m(t_f) = %+.3e                  %s\n', ...
        Y(end,14), pass(abs(Y(end,14)) < 1e-6));

% N5  the control obeys the minimum principle: alpha = -lam_v/|lam_v|, |alpha| = 1
alpha = -lamv ./ max(nlv, realmin);
fprintf('   N5 control law   max||alpha| - 1| = %.2e               %s\n', ...
        max(abs(vecnorm(alpha,2,2) - 1)), pass(max(abs(vecnorm(alpha,2,2)-1)) < 1e-12));

% N6  the second solver agreed (section 5), which is not a theory condition
fprintf('   N6 independent   |dz| = %.2e                          %s\n', ...
        V.dz, pass(~V.moved));

necessary = norm(R1,inf) < 1e-8 && max(abs(Hval)) < 1e-6 && missKm < 100 && ...
            missVms < 10 && abs(Y(end,14)) < 1e-6 && ~V.moved;

%% ========================================================================
%  7. SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat) -- one at a time
%     With section 6, these give a strict STRONG local minimizer among
%     trajectories with the same endpoints.
%% ========================================================================
fprintf('\n7. SUFFICIENCY HYPOTHESES\n');

% S1  strengthened Legendre. For a direction on the unit sphere the second
%     derivative of H in the control, restricted to that sphere, is
%     (T/m)|lam_v| I -- positive definite exactly when |lam_v| > 0.
[minLamV, iLV] = min(nlv);
fprintf('   S1 Legendre      min|lam_v| = %.4e at t/t_f = %.3f   %s\n', ...
        minLamV, tu(iLV)/tu(end), pass(minLamV > 0));

% S2  all-burn really is the PMP control: the switching function stays positive
Qmt = nlv./mm + lamm/cnd;
fprintf('   S2 switching     min Q = %.4e                        %s\n', ...
        min(Qmt), pass(min(Qmt) > 0));

% S3  normality: no abnormal lift of the SAME trajectory (dim S = 1)
gates = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, muStar, struct());
fprintf('   S3 normality     dim S = %d (1 = no abnormal lift)        %s\n', ...
        gates.dimS, pass(gates.dimS == 1));

% S4  no conjugate time in (0, t_f], by the free-time quotiented Jacobi test
cj = it.conj;
fprintf('   S4 conjugate     %s, %d crossing(s), min|det| = %.2e       %s\n', ...
        cj.verdict, gv(cj, 'nCrossings'), min(abs(cj.detScaled)), pass(cj.pass == 1));

sufficient = minLamV > 0 && min(Qmt) > 0 && gates.dimS == 1 && cj.pass == 1;

fprintf('\n   VERDICT: ');
if necessary && sufficient
    fprintf(['strict strong local minimizer among trajectories with the same\n' ...
             '            endpoints -- numerically certified at the sampled times.\n']);
elseif ~necessary
    fprintf(['NOT an extremal to tolerance. The second-order test is meaningless\n' ...
             '            off an extremal, so no minimality is claimed.\n']);
else
    fprintf('an extremal, but a sufficiency hypothesis fails: no minimality claimed.\n');
end

% cross-check the explicit numbers against the shared instrument, so this
% script cannot quietly become a second, unverified implementation
assert(abs(minLamV - gates.minLamV) < 1e-6*max(gates.minLamV,1) && ...
       abs(min(Qmt) - gates.minQmt) < 1e-6*max(gates.minQmt,1), ...
       'the inline gates disagree with mintime_hypothesis_gates');

%% ========================================================================
%  8. INTERACTIVE 3D PLOT  (drag to rotate, scroll to zoom)
%% ========================================================================
T = struct('z', z8, 'sD', sD, 'sA', sA, 'tfDays', day(z8(8)), 'dvKms', dV, ...
           'mfKg', m0kg*(1-mf));
P = plot_transfer_3d(T, B);
fprintf('\n8. Figure %d is rotatable.\n', P.fig.Number);

%% ------------------------------------------------------------------------
function s = pass(c)
% PASS  Verdict text.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function v = gv(s, f)
% GV  Field or NaN.  INPUTS: s; f.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = NaN; end
end
