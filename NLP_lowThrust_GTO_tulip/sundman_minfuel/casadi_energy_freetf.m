function out = casadi_energy_freetf(sigma, rv0, rvf, Tmax, cEx, muStar, X0, U0, tauf0, opts)
% CASADI_ENERGY_FREETF  Free-final-time min-ENERGY Sundman collocation with a
% two-primary clock and Moon-gravity homotopy (CasADi+IPOPT).
%
% Purpose: manufacture a GTO->ELFO min-ENERGY solution (the homotopy root the
% PSR min-fuel pipeline consumes). It is a redesign of CASADI_MINFUEL_SUNDMAN
% along three axes recommended by an external design review (GPT-5.6-terra +
% Gemini 3.1 Pro, 2026-07-13), each targeting a distinct wall that stalled the
% fixed-t_f, single-primary, linear-target-homotopy route at s=0.45:
%
%   (1) TWO-PRIMARY Sundman clock.  The original kappa = r1^pSund uses only the
%       EARTH distance, so as the terminal moves into the Moon's well the mesh
%       starves exactly where the ~mu/r2^3 lunar Hessian spikes.  Here
%           kappa = ( r1^(-q) + (r2/D)^(-q) )^(-p/q)          [moonZone D > 0]
%       a smooth min-distance clock: kappa -> r1^p near Earth, -> (r2/D)^p near
%       the Moon, redistributing nodes into the lunar-capture arc.  q sets the
%       transition sharpness; D (moonZone) sets the crossover radius (~lunar SOI,
%       0.15 ND).  moonZone <= 0 recovers the original single-primary r1^p clock.
%
%   (2) FREE physical t_f with a BANDED KKT.  A free scalar tau_f couples every
%       defect into one dense KKT column (OOM at large N) -- which is why the
%       original FIXES tau_f and pins t(tau_f)=t_f, hitting a "can't-reach-
%       terminal" wall when the target moves Moon-ward.  Instead we float t_f via
%       a slack STATE cScale (Betts' sparse free-time trick):
%           dt/dtau = cScale * kappa,   dcScale/dtau = 0.
%       cScale is one number tied across nodes by LOCAL continuity constraints,
%       so the Jacobian stays banded.  t_f = t(tau_f) is now a free outcome
%       (loose box), so every intermediate target is reachable at some t_f.
%       tau_f itself remains FIXED (= tauf0).
%
%   (3) MOON-GRAVITY homotopy (driver-facing).  The linear Cartesian target walk
%       forced dynamically-inconsistent intermediate rendezvous states; instead
%       hold rvf FIXED and continue muGain: 0 -> 1 scaling ONLY the Moon-gravity
%       acceleration term  -muGain*muStar*rr/r3.  muStar in the Coriolis /
%       centrifugal / Earth terms and in the frame is UNCHANGED (scaling it would
%       move the barycenter and rotate the frame under the boundary conditions).
%       muGain=0 is a well-less near-2-body transfer that converges easily;
%       muGain=1 is the true problem.
%
% State  x = [r(3); v(3); m; t; cScale]  (9).   Control u = [alpha(3); s]  (4).
% Cone-eliminated thrust = s*Tmax*alpha/m, ||alpha||=1, s in [0,1].
% Objective (Bertrand-Epenoy, physical-time measure dt = cScale*kappa dtau):
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt  [+ tfWeight * t(tau_f) if t_f free]
%   eps=1 -> Int[s^2]dt (ENERGY).
%
% t_f MODE.  Min-ENERGY is only well-posed at a PINNED t_f (unpinned, the energy
% optimum drifts t_f: longer time -> thinner thrust -> lower Int[s^2]dt, so IPOPT
% wanders off the warm start).  So by default pass opts.tfTarget: the constraint
% t(tau_f)=tfTarget is added and the slack state cScale floats to satisfy it --
% a clean single-DOF way to hold a fixed transfer time under the (changing)
% two-primary clock.  Leave tfTarget empty ONLY for a genuinely free-t_f problem
% (e.g. min-time), in which case tfWeight>0 is needed to keep it well-posed.
%
% INPUTS:
%   sigma   - normalized independent-variable nodes [(N+1)x1], 0 -> 1
%   rv0,rvf - initial / target position-velocity (ND) [1x6]
%   Tmax    - max thrust (ND) [scalar]
%   cEx     - exhaust velocity (ND) [scalar]   (Tmax/cEx = mdot per unit throttle)
%   muStar  - CR3BP mass ratio [scalar]
%   X0      - warm-start states.  [8x(N+1)] ([r;v;m;t]) or [9x(N+1)] ([...;cScale]).
%             If 8 rows, a cScale row = opts.c0 (default 1) is appended, which
%             reproduces the original dt/dtau=kappa timing of the warm start.
%   U0      - warm-start controls [4x(N+1)] ([alpha;s])
%   tauf0   - fixed total regularized length [scalar]
%   opts    - (optional) struct of knobs (all have defaults):
%       .muGain   gravity-homotopy scale on the Moon term, [0,1]   [1]
%       .pSund    Sundman power p                                  [1.5]
%       .qSund    two-primary transition sharpness q               [4]
%       .moonZone crossover radius D (ND); <=0 -> single-primary   [0.15]
%       .epsilon  homotopy: 1=energy, 0=fuel                       [1]
%       .tfTarget pin t(tau_f)=tfTarget (well-posed energy); []=free[warm-start t_f]
%       .tfWeight penalty rho on t_f, used ONLY when tfTarget=[]    [0]
%       .cBox     [lo hi] bounds on cScale (time-scale box)        [0.2 5]
%       .c0       cScale seed when appending to an 8-row X0        [1]
%       .tfCapMult t_f upper bound = tfCapMult * (warm-start t_f)  [4]
%       .maxIter  IPOPT cap                                        [3000]
%       .warmTight tight warm start (sharpen) vs loose (continue)  [false]
%
% OUTPUTS:
%   out - struct: .X [9x(N+1)] .U [4x(N+1)] .tauf .mf .tf .cScale
%         .maxDefect .maxUnit .switches .edge .success .ipoptStatus
%         .lamDef [9xN] discrete costates .lamAll .primerAlignDeg .lamMassEnd
%
% REFERENCES:
%   [1] Bertrand & Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," Optim. Control Appl. Methods 23 (2002).
%   [2] Betts, "Practical Methods for Optimal Control...," SIAM (2010) -- sparse
%       free-final-time via a constant slack state.
%   [3] Sundman regularization; two-primary min-distance clock.
%   [4] casadi_minfuel_sundman.m (the fixed-t_f single-primary predecessor).

if nargin < 10 || isempty(opts), opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
muGain   = gd('muGain',   1);
pSund    = gd('pSund',    1.5);
qSund    = gd('qSund',    4);
moonZone = gd('moonZone', 0.15);
epsilon  = gd('epsilon',  1);      % 1 = energy
tfTarget = gd('tfTarget', []);     % [] -> free t_f; else pin t(tau_f)=tfTarget
tfWeight = gd('tfWeight', 0);
cBox     = gd('cBox',     [0.2 5]);
c0       = gd('c0',       1);
tfCapMult= gd('tfCapMult',4);
maxIter  = gd('maxIter',  3000);
warmTight= gd('warmTight',false);  % continuation moves want LOOSE (see below)

cpath = getenv('CASADI_PATH');
if isempty(cpath), cpath = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cpath);
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
dsig  = diff(sigma).';                       % 1 x N

% --- warm start: ensure a 9-row [r;v;m;t;cScale] state guess ----------------
assert(size(X0,1) >= 8, 'X0 must have >=8 rows ([r;v;m;t]); got %d', size(X0,1));
if size(X0,1) == 8
    X0 = [X0; c0*ones(1, size(X0,2))];       % append constant cScale row
end
tf_ws = X0(8,end);                            % warm-start transfer time
tfCap = tfCapMult * max([tf_ws; tfTarget(:); eps]);

% --- symbolic free-time two-primary Sundman dynamics dX/dtau ----------------
x = MX.sym('x', 9);  u = MX.sym('u', 4);
r = x(1:3);  v = x(4:6);  m = x(7);  cs = x(9);   al = u(1:3);  s = u(4);
dd = [r(1)+muStar; r(2); r(3)];              % vector from Earth
rr = [r(1)-1+muStar; r(2); r(3)];            % vector from Moon
d2 = dd.'*dd + 1e-12;   e2 = rr.'*rr + 1e-12;
r1 = sqrt(d2);  r2 = sqrt(e2);
d3 = d2^1.5;    r3 = e2^1.5;                  % guarded denoms
% Moon-gravity homotopy: scale ONLY the Moon acceleration term (frame fixed).
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muGain*muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (s*Tmax/m)*al;
mdot  = -(Tmax/cEx)*s;
% two-primary min-distance clock (moonZone<=0 -> original single-primary r1^p)
if moonZone > 0
    kappa = ( r1^(-qSund) + (r2/moonZone)^(-qSund) )^(-pSund/qSund);
else
    kappa = r1^pSund;
end
% free physical time via the constant slack state cs:  dt/dtau = cs*kappa
fdyn  = Function('f', {x,u}, {[ cs*kappa*[v; accel; mdot; 1]; 0 ]});  % dX/dtau (9x1)
Fmap  = fdyn.map(nN);
% energy->fuel integrands on the physical-time measure dt = cs*kappa dtau
gint  = Function('g', {x,u}, {[s*cs*kappa; s*(1-s)*cs*kappa]});
Gmap  = gint.map(nN);

% --- NLP --------------------------------------------------------------------
opti = Opti();
X    = opti.variable(9, nN);
U    = opti.variable(4, nN);
tauf = tauf0;                                 % FIXED (t_f floats via cScale)
F    = Fmap(X, U);                            % 9 x nN, = dX/dtau

% trapezoidal defects in sigma: dX/dsigma = tauf * dX/dtau
D = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,9,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);

% unit thrust direction
opti.subject_to((sum(U(1:3,:).^2, 1) - 1).' == 0);

% bounds (explicit two-sided); cScale in its box, t in [0, tfCap]
lbX = repmat([-3;-3;-3;-12;-12;-12;0.3;0;    cBox(1)], 1, nN);
ubX = repmat([ 3; 3; 3; 12; 12; 12;1.0;tfCap;cBox(2)], 1, nN);
opti.subject_to(X(:) >= lbX(:));   opti.subject_to(X(:) <= ubX(:));
lbU = repmat([-1.1;-1.1;-1.1;0], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1;1], 1, nN);
opti.subject_to(U(:) >= lbU(:));   opti.subject_to(U(:) <= ubU(:));

% boundary conditions; cScale free (pinned by tfTarget below, if given)
opti.subject_to(X(1:6,1) == rv0(:));   opti.subject_to(X(7,1) == 1);   opti.subject_to(X(8,1) == 0);
opti.subject_to(X(1:6,nN) == rvf(:));
if ~isempty(tfTarget)
    opti.subject_to(X(8,nN) == tfTarget);     % PIN t_f (well-posed energy); cScale floats
end

% objective: J = Int[s]dt - eps*Int[s(1-s)]dt  (+ tfWeight*t_f only if t_f free)
G    = Gmap(X, U);                            % 2 x nN
IntF = tauf * sum((dsig/2).*(G(1,1:end-1) + G(1,2:end)));
IntS = tauf * sum((dsig/2).*(G(2,1:end-1) + G(2,2:end)));
Jobj = IntF - epsilon*IntS;
if isempty(tfTarget), Jobj = Jobj + tfWeight*X(8,nN); end
opti.minimize(Jobj);
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
    % TIGHT: re-solve AT a converged point (sharpen). Hug bounds, small barrier.
    p.ipopt.mu_strategy                 = 'monotone';
    p.ipopt.warm_start_init_point       = 'yes';
    p.ipopt.mu_init                     = 1e-4;
    p.ipopt.warm_start_bound_push       = 1e-9;
    p.ipopt.warm_start_bound_frac       = 1e-9;
    p.ipopt.warm_start_slack_bound_push = 1e-9;
    p.ipopt.warm_start_slack_bound_frac = 1e-9;
    p.ipopt.warm_start_mult_bound_push  = 1e-9;
else
    % LOOSE: a genuine continuation move (a muGain step, or a t_f-floating
    % re-solve). Monotone barrier, default bound push, larger initial barrier so
    % IPOPT has room to move (tight pinning makes inf_du blow up on a real step).
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
ss = Us(4,:);

% --- KKT multipliers -> discrete costates + PMP primer check (9-state) -------
lamDef = [];  primerAlignDeg = NaN;  lamMassEnd = NaN;
if numel(lamAll) >= 9*N
    lamDef = reshape(lamAll(1:9*N), 9, N);            % [9 x N] discrete costates
    lamV   = lamDef(4:6, :);
    primer = -lamV ./ max(sqrt(sum(lamV.^2,1)), 1e-12);
    aMid   = 0.5*(Us(1:3,1:end-1) + Us(1:3,2:end));
    aMid   = aMid ./ max(sqrt(sum(aMid.^2,1)), 1e-12);
    burn   = (Us(4,1:end-1) > 0.5) & (Us(4,2:end) > 0.5);
    if any(burn)
        cang = sum(primer(:,burn).*aMid(:,burn), 1);
        if mean(cang) < 0, cang = -cang; end
        primerAlignDeg = mean(acosd(min(max(cang,-1),1)));
    end
    lamMassEnd = lamDef(7,end);
end

out = struct('X', Xs, 'U', Us, 'tauf', tauf, 'mf', Xs(7,end), ...
             'tf', Xs(8,end), 'cScale', Xs(9,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'switches', sum(abs(diff(ss > 0.5))), ...
             'edge', mean(ss > 0.95 | ss < 0.05), ...
             'lamDef', lamDef, 'lamAll', lamAll, ...
             'primerAlignDeg', primerAlignDeg, 'lamMassEnd', lamMassEnd, ...
             'success', success, 'ipoptStatus', status);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
