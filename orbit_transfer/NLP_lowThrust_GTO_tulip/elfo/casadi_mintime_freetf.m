function out = casadi_mintime_freetf(sigma, rv0, rvf, Tmax, cEx, muStar, X0, U0, tauf0, opts)
% CASADI_MINTIME_FREETF  Hard all-burn (s==1) free-final-time min-TIME Sundman
% collocation with a two-primary clock (CasADi+IPOPT).
%
% Sibling of CASADI_ENERGY_FREETF: identical two-primary Sundman clock and cScale
% free-t_f slack, but throttle is HARD-PINNED s==1 (control is steering alpha
% only, ||alpha||=1) and the objective is min t(tau_f). No epsilon, no tfTarget,
% no throttle integrand -- the objective IS min-time, so IPOPT drives t_f down
% until the rendezvous BC can no longer be met all-burn. That floor is tfMin.
%
% State  x = [r(3); v(3); m; t; cScale]  (9).  Control u = alpha (3), ||alpha||=1.
% Always-burn: thrust accel = (Tmax/m)*alpha,  mdot = -Tmax/cEx.
%
% INPUTS:
%   sigma   - normalized independent-variable nodes [(N+1)x1], 0 -> 1
%   rv0,rvf - initial / target position-velocity (ND) [1x6]
%   Tmax    - max thrust (ND) [scalar]
%   cEx     - exhaust velocity (ND) [scalar]
%   muStar  - CR3BP mass ratio [scalar]
%   X0      - warm-start states [8x(N+1)] ([r;v;m;t]) or [9x(N+1)] ([...;cScale])
%   U0      - warm-start control [3x(N+1)] (alpha) or [4x(N+1)] ([alpha;s]); a
%             4th throttle row is dropped (s==1 is enforced here)
%   tauf0   - fixed total regularized length [scalar]
%   opts    - (optional) struct: .pSund[1.5] .qSund[4] .moonZone[0.15]
%             .cBox[[0.10 8]] .c0[1] .tfCapMult[4] .maxIter[3000] .warmTight[false]
%
% OUTPUTS:
%   out - struct: .X [9x(N+1)] .U [3x(N+1)] .tauf .tf(=tfMin) .cScale .mf
%         .maxDefect .maxUnit .minR1 .tMonotone (logical) .lamDef [9xN] .lamAll
%         .primerAlignDeg .lamMassEnd .success .ipoptStatus
%
% REFERENCES:
%   [1] casadi_energy_freetf.m (the energy/fuel sibling this mirrors).
%   [2] Betts, "Practical Methods for Optimal Control," SIAM (2010) -- sparse
%       free-final-time via a constant slack state.

if nargin < 10 || isempty(opts), opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
pSund    = gd('pSund',    1.5);
qSund    = gd('qSund',    4);
moonZone = gd('moonZone', 0.15);
cBox     = gd('cBox',     [0.10 8]);
c0       = gd('c0',       1);
tfCapMult= gd('tfCapMult',4);
maxIter  = gd('maxIter',  3000);
warmTight= gd('warmTight',false);

cpath = getenv('CASADI_PATH');
if isempty(cpath), cpath = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cpath);
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
dsig  = diff(sigma).';

% --- warm start: 9-row state, 3-row control ---------------------------------
assert(size(X0,1) >= 8, 'X0 must have >=8 rows ([r;v;m;t]); got %d', size(X0,1));
if size(X0,1) == 8
    X0 = [X0; c0*ones(1, size(X0,2))];
end
assert(size(U0,1) >= 3, 'U0 must have >=3 rows (alpha); got %d', size(U0,1));
U0 = U0(1:3,:);                               % drop any throttle row (s==1)
tf_ws = X0(8,end);
tfCap = tfCapMult * max(tf_ws, eps);

% --- symbolic always-burn two-primary Sundman dynamics dX/dtau --------------
x = MX.sym('x', 9);  u = MX.sym('u', 3);
r = x(1:3);  v = x(4:6);  m = x(7);  cs = x(9);   al = u;
dd = [r(1)+muStar; r(2); r(3)];               % vector from Earth
rr = [r(1)-1+muStar; r(2); r(3)];             % vector from Moon
d2 = dd.'*dd + 1e-12;   e2 = rr.'*rr + 1e-12;
r1 = sqrt(d2);  r2 = sqrt(e2);
d3 = d2^1.5;    r3 = e2^1.5;
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;   % full gravity (muGain=1)
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (Tmax/m)*al;                % s == 1
mdot  = -(Tmax/cEx);                          % s == 1
if moonZone > 0
    kappa = ( r1^(-qSund) + (r2/moonZone)^(-qSund) )^(-pSund/qSund);
else
    kappa = r1^pSund;
end
fdyn  = Function('f', {x,u}, {[ cs*kappa*[v; accel; mdot; 1]; 0 ]});  % 9x1
Fmap  = fdyn.map(nN);

% --- NLP --------------------------------------------------------------------
opti = Opti();
X    = opti.variable(9, nN);
U    = opti.variable(3, nN);
tauf = tauf0;                                 % FIXED (t_f floats via cScale)
F    = Fmap(X, U);

% trapezoidal defects in sigma: dX/dsigma = tauf * dX/dtau
D = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,9,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);

% unit steering direction
opti.subject_to((sum(U(1:3,:).^2, 1) - 1).' == 0);

% bounds (explicit two-sided); cScale in its box, t in [0, tfCap]
lbX = repmat([-3;-3;-3;-12;-12;-12;0.3;0;    cBox(1)], 1, nN);
ubX = repmat([ 3; 3; 3; 12; 12; 12;1.0;tfCap;cBox(2)], 1, nN);
opti.subject_to(X(:) >= lbX(:));   opti.subject_to(X(:) <= ubX(:));
lbU = repmat([-1.1;-1.1;-1.1], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1], 1, nN);
opti.subject_to(U(:) >= lbU(:));   opti.subject_to(U(:) <= ubU(:));

% boundary conditions; cScale free
opti.subject_to(X(1:6,1) == rv0(:));   opti.subject_to(X(7,1) == 1);   opti.subject_to(X(8,1) == 0);
opti.subject_to(X(1:6,nN) == rvf(:));

% objective: MIN-TIME
opti.minimize(X(8,nN));
opti.set_initial(X, X0);
opti.set_initial(U, U0);

p = struct;
p.print_time      = true;
p.ipopt.max_iter  = maxIter;
p.ipopt.tol       = 1e-7;
p.ipopt.constr_viol_tol = 1e-7;
p.ipopt.acceptable_tol  = 1e-5;
p.ipopt.acceptable_iter = 15;
p.ipopt.nlp_scaling_method = 'gradient-based';
p.ipopt.linear_solver = 'mumps';
p.ipopt.print_level = 5;
if warmTight
    p.ipopt.mu_strategy                 = 'monotone';
    p.ipopt.warm_start_init_point       = 'yes';
    p.ipopt.mu_init                     = 1e-4;
    p.ipopt.warm_start_bound_push       = 1e-9;
    p.ipopt.warm_start_bound_frac       = 1e-9;
    p.ipopt.warm_start_slack_bound_push = 1e-9;
    p.ipopt.warm_start_slack_bound_frac = 1e-9;
    p.ipopt.warm_start_mult_bound_push  = 1e-9;
else
    p.ipopt.mu_strategy           = 'monotone';
    p.ipopt.warm_start_init_point = 'yes';
    p.ipopt.mu_init               = 0.1;
end
opti.solver('ipopt', p);

success = true;  status = 'solved';
try
    sol = opti.solve();
    Xs = sol.value(X);  Us = sol.value(U);
    lamAll = full(sol.value(opti.lam_g));
    status = char(opti.return_status());
catch solveErr
    Xs = opti.debug.value(X);  Us = opti.debug.value(U);
    try
        lamAll = full(opti.debug.value(opti.lam_g));
    catch
        lamAll = [];
    end
    success = false;  status = solveErr.message;
end

% metrics
Fs = full(Fmap(Xs, Us));
Dd = Xs(:,2:end) - Xs(:,1:end-1) - tauf*(repmat(dsig,9,1)/2).*(Fs(:,1:end-1) + Fs(:,2:end));

% --- KKT multipliers -> discrete costates + PMP primer (all burn) -----------
lamDef = [];  primerAlignDeg = NaN;  lamMassEnd = NaN;
if numel(lamAll) >= 9*N
    lamDef = reshape(lamAll(1:9*N), 9, N);
    lamV   = lamDef(4:6, :);
    primer = -lamV ./ max(sqrt(sum(lamV.^2,1)), 1e-12);
    aMid   = 0.5*(Us(1:3,1:end-1) + Us(1:3,2:end));
    aMid   = aMid ./ max(sqrt(sum(aMid.^2,1)), 1e-12);
    cang   = sum(primer.*aMid, 1);            % every interval is a burn (s==1)
    if mean(cang) < 0, cang = -cang; end
    primerAlignDeg = mean(acosd(min(max(cang,-1),1)));
    lamMassEnd = lamDef(7,end);
end

% physicality diagnostics
r1N   = sqrt((Xs(1,:)+muStar).^2 + Xs(2,:).^2 + Xs(3,:).^2);
minR1 = min(r1N);
tMonotone = all(diff(Xs(8,:)) > 0);

out = struct('X', Xs, 'U', Us, 'tauf', tauf, 'tf', Xs(8,end), ...
             'cScale', Xs(9,end), 'mf', Xs(7,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'minR1', minR1, 'tMonotone', tMonotone, ...
             'lamDef', lamDef, 'lamAll', lamAll, ...
             'primerAlignDeg', primerAlignDeg, 'lamMassEnd', lamMassEnd, ...
             'success', success, 'ipoptStatus', status);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
