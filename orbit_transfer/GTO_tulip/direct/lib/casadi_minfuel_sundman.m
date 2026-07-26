function out = casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsilon, warmTight, opts)
% CASADI_MINFUEL_SUNDMAN  Sundman-regularized min-fuel collocation (CasADi+IPOPT).
%
% Path A: fix the near-perigee ill-conditioning that stalled the plain
% collocation by changing the independent variable from time t to a
% regularizing variable tau, with the Sundman relation
%     dt/dtau = kappa(r) = r1^pSund,     r1 = ||r - r_Earth||.
% Every state ODE is multiplied by kappa, so the near-perigee gravity Hessian
% terms (~1/r^3) that blew up IPOPT's exact Hessian become r1^(pSund-3)
% (bounded for pSund ~ 3, mild for pSund ~ 1.5). A UNIFORM mesh in tau also
% concentrates nodes near perigee in TIME automatically. Time t is carried as
% an 8th state; the fixed transfer time is the terminal constraint
% t(tau_f) = tf. The total regularized length tau_f is held FIXED (= tauf0, the
% warm-start value), NOT a decision variable: a free scalar tau_f multiplies
% every collocation defect, producing one dense KKT column -> catastrophic
% MUMPS fill-in / OOM at large N. Fixing tau_f and enforcing t(tau_f)=tf on the
% carried time state instead keeps the Jacobian sparse; the trajectory adjusts
% so that int(kappa dtau) = tf.
%
% The objective is the Bertrand-Epenoy energy->fuel homotopy in epsilon:
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (physical-time measure dt=kappa dtau)
%   eps=1 -> Int[s^2]dt (energy, strictly convex, smooth ramp)
%   eps=0 -> Int[s]dt   (fuel, linear -> bang-bang; equals propellant up to a
%                        positive constant, since m(tf)=1-(Tmax/c)Int[s]dt).
% Sweep eps 1->0, warm-starting each solve from the last (see run_sundman_*).
%
% State  x = [r(3); v(3); m; t]  (8).   Control u = [alpha(3); s]  (4).
% Cone-eliminated: thrust = s*Tmax*alpha/m, ||alpha|| = 1, s in [0,1].
%
% INPUTS:
%   sigma   - normalized independent-variable nodes [(N+1)x1], 0 -> 1
%   tf      - fixed transfer TIME (ND) [scalar]
%   rv0,rvf - initial / target position-velocity (ND) [1x6]
%   Tmax,c,muStar - dynamics constants [scalars]
%   X0      - warm-start states [8x(N+1)] ([r;v;m;t])
%   U0      - warm-start controls [4x(N+1)] ([alpha;s])
%   tauf0   - fixed total regularized length [scalar]
%   pSund   - Sundman power [scalar, default 1.5]
%   maxIter - IPOPT max iterations [scalar, default 3000]
%   epsilon - homotopy parameter in [0,1]: 0=fuel, 1=energy [scalar, default 0]
%   warmTight - true (default): tight warm start for re-solving AT a
%           near-bang-bang solution (homotopy sharpening); false: loose
%           (adaptive barrier, default bound_push) for a genuine continuation
%           move such as an energy re-solve at a shifted t_f [logical]
%   opts    - (optional) struct: .vBox position/velocity... see below
%           .vBox - velocity box half-width, ND [scalar, default 12]
%           .rBox - position box half-width, ND [scalar, default 3]
%           .returnModel - (default false) true -> ADDITIONALLY attach
%           out.model = struct('opti',opti,'creg',creg) with the live solved
%           CasADi Opti object and a constraint registry creg (struct array,
%           fields .label[char] .rows[1xk] row range into opti.g) recording
%           the 'defect', 'betaNorm', 'thrLo', 'thrHi' constraint groups, for
%           the generic FOC/KKT gate (verify_common/foc_check.m). Purely
%           additive: with the flag absent/false, X/U are byte-identical and
%           out.model is absent (Task 8, 2026-07-25).
%           Omit or pass [] / struct() for the nominal (byte-identical) bounds.
%
% OUTPUTS:
%   out - struct: .X [8x(N+1)] .U [4x(N+1)] .tauf .mf .maxDefect .maxUnit
%         .switches .edge (bang-bang node fraction) .success .ipoptStatus
%         .lamDef [8xN] discrete costates (defect-constraint KKT multipliers,
%           [lam_r;lam_v;lam_m;lam_t] per interval, up to a positive mesh-weight
%           scaling and a global sign), .lamAll (full stacked g-multiplier),
%         .primerAlignDeg (mean angle between the NLP thrust direction and the
%           costate primer -lam_v/||lam_v|| on burn arcs; ~0 certifies PMP),
%         .lamMassEnd (terminal mass-costate proxy; ~0 is the transversality)
%         .boundSat - struct('minSlack',s,'worst',label,'hit',logical): the
%           tightest nonphysical-box slack at INTERIOR nodes (BCs pin the
%           endpoints by construction); .hit true warns the box may be
%           binding and should be widened via opts before trusting the result
%         .regHistory [1xnIter or []] - IPOPT's per-iteration Hessian
%           regularization delta_w (st.iterations.regularization_size), or []
%           if the CasADi build lacks it; interpret via
%           verify_common/foc_ipopt_inertia.m (same field name/shape as
%           PSR/lib/casadi_minfuel_sundman.m and psr_ipopt_certify.m expect)
%         .model (opts.returnModel only) - struct('opti',opti,'creg',creg)
%
% REFERENCES:
%   [1] Bertrand & Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," Optim. Control Appl. Methods 23 (2002).
%   [2] Sundman regularization of the two/three-body problem; e.g. dt = r dtau.
%   [3] Andersson et al., "CasADi," Math. Prog. Comp. 11 (2019); Wachter &
%       Biegler (IPOPT), Math. Prog. 106 (2006).
%   [4] earth_elliptic_to_geo/direct/core/casadi_lt_mee.m (returnModel/creg
%       registry pattern this mirrors); GTO_tulip/direct/lib/
%       casadi_minfuel_sundman.m (regHistory capture pattern this mirrors).

if nargin < 11 || isempty(pSund),  pSund  = 1.5;  end
if nargin < 12 || isempty(maxIter), maxIter = 3000; end
if nargin < 13 || isempty(epsilon), epsilon = 0;   end   % 0=fuel, 1=energy
if nargin < 14 || isempty(warmTight), warmTight = true; end  % see IPOPT opts
if nargin < 15 || isempty(opts), opts = struct(); end
vBox = 12;  if isfield(opts,'vBox') && ~isempty(opts.vBox), vBox = opts.vBox; end
rBox = 3;   if isfield(opts,'rBox') && ~isempty(opts.rBox), rBox = opts.rBox; end
returnModel = false;
if isfield(opts,'returnModel') && ~isempty(opts.returnModel), returnModel = opts.returnModel; end
cpath = getenv('CASADI_PATH');
if isempty(cpath), cpath = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cpath);
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
dsig  = diff(sigma).';                       % 1 x N

% --- symbolic Sundman-regularized dynamics dX/dtau = kappa * f -------------
x = MX.sym('x', 8);  u = MX.sym('u', 4);
r = x(1:3);  v = x(4:6);  m = x(7);  al = u(1:3);  s = u(4);
dd = [r(1)+muStar; r(2); r(3)];              % vector from Earth
rr = [r(1)-1+muStar; r(2); r(3)];            % vector from Moon
r1 = sqrt(dd.'*dd + 1e-12);                  % Earth distance (guarded)
d3 = (dd.'*dd + 1e-12)^1.5;  r3 = (rr.'*rr + 1e-12)^1.5;   % guarded denoms
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (s*Tmax/m)*al;
mdot  = -(Tmax/c)*s;
kappa = r1^pSund;
fdyn  = Function('f', {x,u}, {kappa*[v; accel; mdot; 1]});   % dX/dtau (8x1)
Fmap  = fdyn.map(nN);
% energy->fuel homotopy integrands (Bertrand-Epenoy), on the physical-time
% measure dt = kappa dtau:  q_fuel = s*kappa,  q_smooth = s(1-s)*kappa.
gint  = Function('g', {x,u}, {[s*kappa; s*(1-s)*kappa]});
Gmap  = gint.map(nN);

% --- NLP ------------------------------------------------------------------
opti = Opti();
X    = opti.variable(8, nN);
U    = opti.variable(4, nN);
% SOSC/FOC-gate registry (Task 8, additive only): records the row range each
% subject_to group occupies in opti.g, purely for the generic FOC/KKT gate
% (verify_common/foc_check.m) -- zero effect on the solve itself (bracketing
% reads size(opti.g,1), never writes it).
creg = struct('label',{},'rows',{});
% tau-length is FIXED (from the warm start), NOT a decision variable: a free
% scalar tau_f couples to every defect -> a dense KKT column -> catastrophic
% MUMPS fill-in / OOM at large N. Fixed transfer time is still enforced
% exactly by the t-state terminal condition t(tau_end) = tf below; the
% trajectory adjusts so that int(kappa dtau) = tf.
tauf = tauf0;
F    = Fmap(X, U);                           % 8 x nN, = dX/dtau

% trapezoidal defects in sigma: dX/dsigma = tauf * dX/dtau
r0 = size(opti.g,1)+1;
D = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,8,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);
if returnModel, creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1)); end

% unit-direction
r0 = size(opti.g,1)+1;
opti.subject_to((sum(U(1:3,:).^2, 1) - 1).' == 0);
if returnModel, creg(end+1) = struct('label','betaNorm','rows',r0:size(opti.g,1)); end

% bounds (explicit two-sided)
lbX = repmat([-rBox;-rBox;-rBox;-vBox;-vBox;-vBox;0.3;0], 1, nN);
ubX = repmat([ rBox; rBox; rBox; vBox; vBox; vBox;1.0; 2*tf], 1, nN);
opti.subject_to(X(:) >= lbX(:));   opti.subject_to(X(:) <= ubX(:));
lbU = repmat([-1.1;-1.1;-1.1;0], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1;1], 1, nN);
r0 = size(opti.g,1)+1;
opti.subject_to(U(:) >= lbU(:));
% M2 fix (final-review wave): U(:) is a SINGLE vectorized box constraint over
% all 4 control rows (3 direction + 1 throttle), stacked column-major, so the
% raw r0:size(opti.g,1) range covers all 4 rows/column, not just the throttle
% row. Constraint content/order is untouched (still one subject_to over the
% full U(:) box); only the creg METADATA is narrowed here to the throttle
% sub-range -- row 4 of every 4-row column block, i.e. rows 4:4:end of this
% range -- so foc_check's throttle-row bookkeeping (Sd/switching-function
% sign check) doesn't get contaminated by the direction-vector bound rows.
if returnModel
    allRowsLo = r0:size(opti.g,1);
    creg(end+1) = struct('label','thrLo','rows',allRowsLo(4:4:end));
end
r0 = size(opti.g,1)+1;
opti.subject_to(U(:) <= ubU(:));
if returnModel
    allRowsHi = r0:size(opti.g,1);
    creg(end+1) = struct('label','thrHi','rows',allRowsHi(4:4:end));
end

% boundary conditions (fixed transfer TIME via the t-state)
opti.subject_to(X(1:6,1) == rv0(:));   opti.subject_to(X(7,1) == 1);   opti.subject_to(X(8,1) == 0);
opti.subject_to(X(1:6,nN) == rvf(:));  opti.subject_to(X(8,nN) == tf);

% objective: J(eps) = Int[s]dt - eps*Int[s(1-s)]dt  (trapezoid in tau)
%   eps=0 -> Int s dt   = fuel (linear in s -> bang-bang)
%   eps=1 -> Int s^2 dt = energy (strictly convex -> smooth ramp, no restoration)
G    = Gmap(X, U);                            % 2 x nN  [q_fuel; q_smooth]
IntF = tauf * sum((dsig/2).*(G(1,1:end-1) + G(1,2:end)));
IntS = tauf * sum((dsig/2).*(G(2,1:end-1) + G(2,2:end)));
opti.minimize(IntF - epsilon*IntS);
opti.set_initial(X, X0);
opti.set_initial(U, U0);

% IPOPT options: single source in cr3bp_common (Tier-0 extraction 2026-07-26).
% These ~20 assignments were byte-identical across this file,
% casadi_energy_freetf and casadi_mintime_freetf; the helper's header explains
% the two warm-start regimes and why the earth campaign and PSR/lib are
% deliberately NOT sharing it. Gate: cr3bp_common/tests/test_cr3bp_ipopt_opts.m
% asserts the helper reproduces the former inline struct exactly.
p = cr3bp_ipopt_opts(maxIter, warmTight);

opti.solver('ipopt', p);

success = true;  status = 'solved';  regHistory = [];
try
    sol = opti.solve();
    Xs = sol.value(X);  Us = sol.value(U);
    lamAll = full(sol.value(opti.lam_g));
    status = char(opti.return_status());
    % IPOPT per-iteration Hessian regularization (delta_w). At a genuine local
    % min IPOPT's inertia-controlled linear solver adds ZERO regularization at
    % convergence -- the reduced Hessian is PD without correction. Captured
    % here (ported from PSR/lib/casadi_minfuel_sundman.m, which had this but
    % this file's fork had dropped it) so foc_ipopt_inertia / psr_ipopt_certify
    % can read the native (well-scaled) 2nd-order verdict.
    try
        st = sol.stats();
        if isfield(st,'iterations') && isfield(st.iterations,'regularization_size')
            regHistory = st.iterations.regularization_size(:).';
        end
    catch  %#ok<CTCH>
    end
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
Dd = Xs(:,2:end) - Xs(:,1:end-1) - tauf*(repmat(dsig,8,1)/2).*(Fs(:,1:end-1) + Fs(:,2:end));
ss = Us(4,:);

% --- KKT multipliers -> discrete costates + PMP primer-vector check ---------
% The duals of the dynamics-defect constraints ARE the discrete costates
% [lam_r; lam_v; lam_m; lam_t] (one per interval, up to a positive mesh-weight
% scaling and a global sign convention). The primer condition
% alpha* = -lam_v/||lam_v|| is scale-invariant, so comparing it to the NLP
% thrust direction on burn arcs is an independent optimality certificate.
lamDef = [];  primerAlignDeg = NaN;  lamMassEnd = NaN;
if numel(lamAll) >= 8*N
    lamDef = reshape(lamAll(1:8*N), 8, N);          % [8 x N] discrete costates
    lamV   = lamDef(4:6, :);                          % velocity-costate proxy
    primer = -lamV ./ max(sqrt(sum(lamV.^2,1)), 1e-12);
    aMid   = 0.5*(Us(1:3,1:end-1) + Us(1:3,2:end));   % node dirs -> interval mids
    aMid   = aMid ./ max(sqrt(sum(aMid.^2,1)), 1e-12);
    burn   = (Us(4,1:end-1) > 0.5) & (Us(4,2:end) > 0.5);
    if any(burn)
        cang = sum(primer(:,burn).*aMid(:,burn), 1);
        if mean(cang) < 0, cang = -cang; end          % absorb global costate sign
        primerAlignDeg = mean(acosd(min(max(cang,-1),1)));
    end
    lamMassEnd = lamDef(7,end);                        % mass costate ~0 (transversality)
end

% Bound-saturation diagnostic (2026-07-21 triage C4; output-only). Nonphysical
% boxes checked at INTERIOR nodes (BCs pin the endpoints by construction).
Xi = Xs(:,2:end-1);
slk = [ rBox - max(abs(Xi(1:3,:)),[],'all');            % position box
        vBox - max(abs(Xi(4:6,:)),[],'all');            % velocity box
        min(Xi(7,:),[],'all') - 0.3;                    % mass lower
        1.0 - max(Xi(7,:),[],'all') ];                  % mass upper
lbl = {'rBox','vBox','massLo','massHi'};
[minSlack, iw] = min(slk);
boundSat = struct('minSlack', minSlack, 'worst', lbl{iw}, 'hit', minSlack < 1e-4);
if boundSat.hit
    warning('casadi_minfuel_sundman:boundSaturation', ...
        'nonphysical box ''%s'' within %.2g of binding -- widen via opts before trusting', ...
        lbl{iw}, max(minSlack,0));
end

out = struct('X', Xs, 'U', Us, 'tauf', tauf, 'mf', Xs(7,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'switches', sum(abs(diff(ss > 0.5))), ...
             'edge', mean(ss > 0.95 | ss < 0.05), ...
             'lamDef', lamDef, 'lamAll', lamAll, ...
             'primerAlignDeg', primerAlignDeg, 'lamMassEnd', lamMassEnd, ...
             'boundSat', boundSat, ...
             'success', success, 'ipoptStatus', status, 'regHistory', regHistory);

if returnModel
    out.model = struct('opti', opti, 'creg', creg);
end
end
