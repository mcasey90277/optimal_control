function out = casadi_lt_2body(sigma, X0, U0, tauf0, term, opts)
% CASADI_LT_2BODY  Sundman-collocated 2-body low-thrust NLP (CasADi+IPOPT).
%
% Trapezoidal collocation in the Sundman variable tau with a cScale slack state:
%     dt/dtau = cScale * kappa(r),   kappa = ||r||^pSund,   dcScale/dtau = 0,
% tau_f held FIXED (= tauf0) so the KKT stays banded (Betts' sparse free-time
% trick; a free scalar tau_f makes a dense KKT column). State per node
% x = [r(3); v(3); m; t; cScale] (9), control u = [alpha(3); s] (4), cone-
% eliminated thrust = s*Tmax*alpha/m with ||alpha||=1, s in [0,1].
%
% Modes:
%   'mintime' - s == 1 (all-burn restriction; optimal for this transfer),
%               objective min t(tau_f); t_f found via cScale.
%   'fixedtf' - constraint t(tau_f) = opts.tfTarget; Bertrand-Epenoy objective
%               J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (dt = cScale*kappa dtau)
%               eps=1 energy (smooth), eps=0 fuel (bang-bang).
%
% Terminal (term from GEO_TERMINAL): 'fixed' pins X(1:6,end) = term.rvf;
% 'manifold' poses the 5 insertion constraints symbolically.
%
% INPUTS:  sigma [(N+1)x1] 0->1;  X0 [8|9 x N+1] warm start (cScale row appended
%          as 1s if absent);  U0 [4xN+1];  tauf0 [scalar];  term [struct];
%          opts .par .mode .eps .tfTarget .rv0 .maxIter .warmTight .printLevel
% OUTPUTS: out - see header table in PLAN.md Task 7 (X,U,success,maxDefect,
%          lamDef, primerAlignDeg, m_f_kg, dV_kms, switches, edge, ...)
%
% REFERENCES:
%   [1] NLP_lowThrust_GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m (parent).
%   [2] NLP_lowThrust_GTO_tulip/elfo/casadi_energy_freetf.m (cScale pattern).
%   [3] DESIGN.md secs 2-4.
cp = getenv('CASADI_PATH');
if isempty(cp), cp = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cp);
par = opts.par;
d = @(f,v) getdef(opts, f, v);
mode      = d('mode', 'fixedtf');
epsv      = d('eps', 0);
tfTarget  = d('tfTarget', []);
maxIter   = d('maxIter', 1500);
warmTight = d('warmTight', false);
printLvl  = d('printLevel', 0);

N    = numel(sigma) - 1;
dtau = diff(sigma(:)).' * tauf0;                    % [1xN]
if size(X0,1) == 8, X0 = [X0; ones(1, N+1)]; end

opti = casadi.Opti();
X = opti.variable(9, N+1);
U = opti.variable(4, N+1);
m = X(7,:);  t = X(8,:);  cS = X(9,:);  al = U(1:3,:);  s = U(4,:);

% node dynamics f{k} = dX/dtau and the clock row kapAll (for the objective)
f = cell(1, N+1);  kapCell = cell(1, N+1);
for k = 1:N+1
    rk  = X(1:3,k);
    rn2 = rk(1)^2 + rk(2)^2 + rk(3)^2 + 1e-12;
    kap = rn2^(par.pSund/2);
    fd  = lt2b_rhs_time(X(1:8,k), U(:,k), par);     % d/dt of [r v m t]
    f{k} = [cS(k)*kap*fd; 0];                       % d/dtau; cScale constant
    kapCell{k} = kap;
end
kapAll = [kapCell{:}];

% collocation defects (KEEP HANDLES for the duals)
conDef = cell(1, N);
for k = 1:N
    conDef{k} = X(:,k+1) - X(:,k) - (dtau(k)/2)*(f{k} + f{k+1}) == 0;
    opti.subject_to(conDef{k});
end

% control cone + throttle (NEVER chain a<=x<=b -- MATLAB gotcha)
for k = 1:N+1
    opti.subject_to(al(1,k)^2 + al(2,k)^2 + al(3,k)^2 == 1);
end
if strcmp(mode, 'mintime')
    opti.subject_to(s == 1);
else
    opti.subject_to(s(:) >= 0);  opti.subject_to(s(:) <= 1);
end

% generous boxes (review lesson: bounds only block divergence)
% NB: CasADi Opti requires inequality expressions to be a vector or square,
% so non-square/row-vector slices are flattened with (:) -- matches the
% (:) convention already used for box bounds in sundman_minfuel/
% casadi_minfuel_sundman.m (the parent file); still two separate calls per
% bound, never chained.
Xr = X(1:3,:);  Xv = X(4:6,:);
opti.subject_to(Xr(:) >= -5);      opti.subject_to(Xr(:) <= 5);
opti.subject_to(Xv(:) >= -8);      opti.subject_to(Xv(:) <= 8);
opti.subject_to(m(:) >= 0.3);      opti.subject_to(m(:) <= 1.001);
% t upper bound: scale with the problem's own timescale, not a flat
% constant (Task 7c fix-round, 2026-07-17 -- see casadi_lt_mee.m's matching
% comment for the incident this guards against: a fixed t<=300 ceiling
% became a structural infeasibility once tfTarget exceeded it, hidden
% behind machine-precision maxDefect/termErr). Verified a no-op for every
% banked Cartesian config in results/ (max observed tf = 66.67 ND, all
% well under 300); kept here purely for consistency/future-proofing since
% this file shares the same "generous box" pattern casadi_lt_mee.m copied
% from.
if strcmp(mode, 'fixedtf')
    tUB = max(300, 2*tfTarget);
else
    tUB = max(300, 3*X0(8,end));   % mintime: 3x the warm-start seed's t-span
end
opti.subject_to(t(:) >= 0);        opti.subject_to(t(:) <= tUB);
opti.subject_to(cS(:) >= 0.05);    opti.subject_to(cS(:) <= 20);
opti.subject_to(al(:) >= -1.01);   opti.subject_to(al(:) <= 1.01);

% boundary conditions
opti.subject_to(X(1:6,1) == opts.rv0(:));
opti.subject_to(m(1) == 1);
opti.subject_to(t(1) == 0);
switch term.type
    case 'fixed'
        opti.subject_to(X(1:6,end) == term.rvf(:));
    case 'manifold'
        re = X(1:3,end);  ve = X(4:6,end);  a = term.aGeo;
        opti.subject_to(re(3) == 0);
        opti.subject_to(ve(3) == 0);
        opti.subject_to(re(1)^2 + re(2)^2 + re(3)^2 == a^2);
        opti.subject_to(ve(1)^2 + ve(2)^2 + ve(3)^2 == par.mu/a);
        opti.subject_to(re(1)*ve(1) + re(2)*ve(2) + re(3)*ve(3) == 0);
        % PROGRADE GUARD (2026-07-17 triage, both reviewers): the 5 residuals
        % above admit the retrograde GEO circle as a legitimate branch (h_z<0)
        % -- a solve that converges onto it is silently wrong, not a
        % convergence failure. Exclude it without constraining the converged
        % prograde value: in canonical units (mu=1, GEO circular r=v=1) the
        % prograde solution has h_z=+1 exactly, retrograde h_z=-1, so
        % h_min=0.1 is inactive at the solution and cleanly excludes h_z<0.
        h_min = 0.1;
        opti.subject_to(re(1)*ve(2) - re(2)*ve(1) >= h_min);
end

% objective + t_f handling
if strcmp(mode, 'mintime')
    opti.minimize(t(end));
else
    assert(~isempty(tfTarget), 'fixedtf mode requires opts.tfTarget');
    opti.subject_to(t(end) == tfTarget);
    w = cS .* kapAll .* (s - epsv*(s.*(1 - s)));    % homotopy integrand * clock
    opti.minimize(sum((dtau/2) .* (w(1:N) + w(2:N+1))));
end

% warm start + IPOPT
opti.set_initial(X, X0);
opti.set_initial(U, U0);
ip = struct('max_iter', maxIter, 'tol', 1e-9, 'constr_viol_tol', 1e-10, ...
            'print_level', printLvl, 'mu_strategy', 'adaptive', ...
            'linear_solver', 'mumps');
if warmTight
    ip.mu_strategy = 'monotone';  ip.mu_init = 1e-4;
    ip.warm_start_init_point = 'yes';
    ip.warm_start_bound_push = 1e-9;  ip.warm_start_mult_bound_push = 1e-9;
end
opti.solver('ipopt', struct('print_time', printLvl > 0), ip);
success = true;
try
    sol = opti.solve();
catch
    sol = opti.debug;  success = false;
end
st = opti.stats();
status = st.return_status;
success = success && any(strcmp(status, {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));

% extraction + numeric re-check of the defects
Xs = sol.value(X);  Us = sol.value(U);
dmax = 0;  fn = zeros(9, N+1);
for k = 1:N+1
    rn = norm(Xs(1:3,k));
    fn(:,k) = [Xs(9,k) * rn^par.pSund * lt2b_rhs_time(Xs(1:8,k), Us(:,k), par); 0];
end
for k = 1:N
    dk = Xs(:,k+1) - Xs(:,k) - (dtau(k)/2)*(fn(:,k) + fn(:,k+1));
    dmax = max(dmax, max(abs(dk)));
end
lamDef = nan(9, N);
try
    for k = 1:N, lamDef(:,k) = sol.value(opti.dual(conDef{k})); end
catch
end
ss = Us(4,:);
burn = ss > 0.5;
% primer alignment on burn nodes (global costate sign resolved by best fit)
lamV = lamDef(4:6, :);
angs = @(sgn) mean(arrayfun(@(k) real(acosd(max(-1,min(1, ...
        dot(Us(1:3,k), sgn*(-lamV(:,min(k,N)))) / max(norm(lamV(:,min(k,N))),1e-30))))), ...
        find(burn)));
primer = min(angs(1), angs(-1));
switch term.type
    case 'fixed',    termErr = norm(Xs(1:6,end) - term.rvf(:));
    case 'manifold', termErr = max(abs(term.resid(Xs(1:6,end))));
end
mf = Xs(7,end);
out = struct('X', Xs, 'U', Us, 'tauf0', tauf0, 'success', success, ...
    'ipoptStatus', status, 'maxDefect', dmax, ...
    'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), 'termErr', termErr, ...
    'mf', mf, 'm_f_kg', par.m0kg*mf, 'dV_kms', par.c*log(1/mf)*par.VU_kms, ...
    'tf', Xs(8,end), 'switches', sum(abs(diff(burn))), ...
    'edge', mean(ss > 0.95 | ss < 0.05), 'lamDef', lamDef, ...
    'primerAlignDeg', primer, 'lamMassEnd', lamDef(7,end));
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
% GETDEF  Optional-field default (mirrors campaign helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
