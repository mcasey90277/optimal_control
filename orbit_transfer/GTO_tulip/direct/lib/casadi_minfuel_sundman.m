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
% !! FIXED tau_f IS NOT EQUIVALENT TO THE FIXED-TIME PROBLEM (external review
% 2026-09-06, E1; measured 09-06/07, GTO_tulip/direct/certify/
% probe_e1_free_tauf.m). Every physical trajectory has its own regularized
% length tau_f = int_0^tf dt/kappa, so fixing tau_f = tauf0 ADDS the
% isoperimetric constraint int dt/kappa = tauf0 inherited from the seed, and
% solutions of the fixed-tau_f NLP are extremals of that restricted problem.
% On the campaign's dV-t_f front every stored row carried the 1.15x seed's
% tauf0 = 151.68 whatever its t_f, and releasing tau_f at the same pinned t_f
% moved m_f by +3e-5 .. +9.4e-3 (up to 142 g of propellant, cScale 0.98-1.05).
% First-order PMP checks cannot detect this (the control minimization is
% unaffected; only the costate ODE gains a -grad(kappa)*(K/kappa) term).
%
% THE FIX IS IN THIS FILE: opts.freeTauf = true (ported 2026-09-07 from
% casadi_energy_freetf, Betts' sparse free-time trick). A constant slack STATE
% cScale (9th row) multiplies the clock,
%     dt/dtau = cScale * kappa,   dcScale/dtau = 0,
% so the regularized length is effectively cScale*tauf0 while the KKT matrix
% stays banded (a free scalar tau_f would couple to every defect -> one dense
% column -> MUMPS fill-in). t(tau_f) = tf is still pinned exactly, so this IS
% the fixed-time problem. Off (default) the solver is byte-identical to the
% fixed-tau_f engine every stored artifact was built with. Outputs keep the
% 8-row X/lamDef contract in either mode; the slack row is exposed as
% out.cScale (and out.X9/out.lamDef9 for the 9-state FOC gate, manifest
% 'tulip_free'). A free-tau_f solution with slack c is EXACTLY a fixed-tau_f
% KKT point at tauf0' = c*tauf0 (rows 1:8 of the defects and the objective
% coincide), so out.tauf = c*tauf0 can be stored as the effective tauf0 of an
% 8-state artifact and re-solved by either mode.
%
% The objective is the Bertrand-Epenoy energy->fuel homotopy in epsilon:
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (physical-time measure dt=kappa dtau)
%   eps=1 -> Int[s^2]dt (energy, strictly convex, smooth ramp)
%   eps=0 -> Int[s]dt   (fuel, linear -> bang-bang; equals propellant up to a
%                        positive constant, since m(tf)=1-(Tmax/c)Int[s]dt).
% Sweep eps 1->0, warm-starting each solve from the last (see minfuel_at_tf,
% the per-t_f driver, and sundman_homotopy).
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
%           (monotone barrier with larger mu_init, default bound_push -- see
%           cr3bp_ipopt_opts) for a genuine continuation
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
%           .freeTauf - (default false) true -> release the regularized
%           length through the cScale slack state (see the header block);
%           t(tau_f) = tf stays pinned. .cBox [lo hi] box on cScale
%           [default 0.2 5]; .c0 cScale seed when X0 has 8 rows [default 1].
%           X0 may carry 9 rows (a previous free-tau_f solution) in this mode.
%           .lamG0 - initial constraint multipliers for opti.lam_g (the
%           .lamAll of a previous solve of the SAME NLP: same mode, same
%           mesh). cr3bp_ipopt_opts sets warm_start_init_point = 'yes', but
%           without this CasADi hands IPOPT ZERO multipliers, so a "sit-still"
%           re-solve at a converged point begins with a dual jump (inf_du
%           ~1e4 at iteration 1) and can wander for hundreds of iterations or
%           leave the basin (7 of 17 free front rows did, 2026-09-07). With
%           lamG0 the re-solve is a genuine warm start. Ignored with a
%           warning if its length does not match size(opti.g,1).
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
%           PSR/lib/casadi_minfuel_sundman.m and ipopt_certify.m expect)
%         .model (opts.returnModel only) - struct('opti',opti,'creg',creg);
%           in freeTauf mode also .model.manifest = 'tulip_free' (nx = 9)
%         .cScale (freeTauf only) - the slack value; .tauf = cScale*tauf0 is
%           then the EFFECTIVE regularized length; .X9 [9x(N+1)] and
%           .lamDef9 [9xN] carry the full free-time primal/costates
%
% REFERENCES:
%   [1] Bertrand & Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," Optim. Control Appl. Methods 23 (2002).
%   [2] Sundman regularization of the two/three-body problem; e.g. dt = r dtau.
%   [3] Andersson et al., "CasADi," Math. Prog. Comp. 11 (2019); Wachter &
%       Biegler (IPOPT), Math. Prog. 106 (2006).
%   [4] earth_elliptic_to_geo/direct/core/casadi_lt_mee.m (returnModel/creg
%       registry pattern this mirrors); the regHistory capture was ported from
%       the former PSR/lib copy of this file (dissolved 2026-07-26).
%   [5] GTO_ELFO/direct/elfo/casadi_energy_freetf.m -- the free-tau_f (cScale)
%       sibling the freeTauf branch was ported from; doc/reviews/
%       direct_core_chain_gpt6astra_review_2026-09-06.md (E1) and
%       results/e1_freetauf/ for why it matters.
%   [6] Betts, "Practical Methods for Optimal Control...," SIAM (2010) --
%       sparse free-final-time via a constant slack state.

if nargin < 11 || isempty(pSund),  pSund  = 1.5;  end
if nargin < 12 || isempty(maxIter), maxIter = 3000; end
if nargin < 13 || isempty(epsilon), epsilon = 0;   end   % 0=fuel, 1=energy
if nargin < 14 || isempty(warmTight), warmTight = true; end  % see IPOPT opts
if nargin < 15 || isempty(opts), opts = struct(); end
vBox = 12;  if isfield(opts,'vBox') && ~isempty(opts.vBox), vBox = opts.vBox; end
rBox = 3;   if isfield(opts,'rBox') && ~isempty(opts.rBox), rBox = opts.rBox; end
returnModel = false;
if isfield(opts,'returnModel') && ~isempty(opts.returnModel), returnModel = opts.returnModel; end
freeTauf = false;
if isfield(opts,'freeTauf') && ~isempty(opts.freeTauf), freeTauf = logical(opts.freeTauf); end
cBox = [0.2 5];  if isfield(opts,'cBox') && ~isempty(opts.cBox), cBox = opts.cBox; end
c0   = 1;        if isfield(opts,'c0')   && ~isempty(opts.c0),   c0   = opts.c0;   end
lamG0 = [];      if isfield(opts,'lamG0') && ~isempty(opts.lamG0), lamG0 = opts.lamG0(:); end
cpath = getenv('CASADI_PATH');
if isempty(cpath), cpath = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cpath);
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
dsig  = diff(sigma).';                       % 1 x N

% --- symbolic Sundman-regularized dynamics dX/dtau = kappa * f -------------
% nx = 8 ([r;v;m;t]) fixed tau_f; nx = 9 ([r;v;m;t;cScale]) free tau_f, where
% the constant slack cs multiplies the clock: dX/dtau = cs*kappa*f, dcs/dtau=0.
if freeTauf, nx = 9; else, nx = 8; end
x = MX.sym('x', nx);  u = MX.sym('u', 4);
r = x(1:3);  v = x(4:6);  m = x(7);  al = u(1:3);  s = u(4);
if freeTauf, cs = x(9); else, cs = 1; end
dd = [r(1)+muStar; r(2); r(3)];              % vector from Earth
rr = [r(1)-1+muStar; r(2); r(3)];            % vector from Moon
r1 = sqrt(dd.'*dd + 1e-12);                  % Earth distance (guarded)
d3 = (dd.'*dd + 1e-12)^1.5;  r3 = (rr.'*rr + 1e-12)^1.5;   % guarded denoms
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (s*Tmax/m)*al;
mdot  = -(Tmax/c)*s;
kappa = r1^pSund;
if freeTauf
    fdyn = Function('f', {x,u}, {[cs*kappa*[v; accel; mdot; 1]; 0]});   % dX/dtau (9x1)
else
    fdyn = Function('f', {x,u}, {kappa*[v; accel; mdot; 1]});           % dX/dtau (8x1)
end
Fmap  = fdyn.map(nN);
% energy->fuel homotopy integrands (Bertrand-Epenoy), on the physical-time
% measure dt = cs*kappa dtau:  q_fuel = s*cs*kappa,  q_smooth = s(1-s)*cs*kappa.
gint  = Function('g', {x,u}, {[s*cs*kappa; s*(1-s)*cs*kappa]});
Gmap  = gint.map(nN);

% --- warm start: nx-row state guess ----------------------------------------
assert(size(X0,1) >= 8, 'casadi_minfuel_sundman: X0 must have >= 8 rows ([r;v;m;t]); got %d', size(X0,1));
if freeTauf && size(X0,1) == 8
    X0 = [X0; c0*ones(1, size(X0,2))];           % append the constant slack row
elseif ~freeTauf && size(X0,1) > 8
    X0 = X0(1:8,:);                              % a free-tau_f seed into the fixed engine
end

% --- NLP ------------------------------------------------------------------
opti = Opti();
X    = opti.variable(nx, nN);
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
F    = Fmap(X, U);                           % nx x nN, = dX/dtau

% trapezoidal defects in sigma: dX/dsigma = tauf * dX/dtau
r0 = size(opti.g,1)+1;
D = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,nx,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);
if returnModel, creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1)); end

% unit-direction
r0 = size(opti.g,1)+1;
opti.subject_to((sum(U(1:3,:).^2, 1) - 1).' == 0);
if returnModel, creg(end+1) = struct('label','betaNorm','rows',r0:size(opti.g,1)); end

% bounds (explicit two-sided); cScale in its own box in freeTauf mode
lbX = repmat([-rBox;-rBox;-rBox;-vBox;-vBox;-vBox;0.3;0], 1, nN);
ubX = repmat([ rBox; rBox; rBox; vBox; vBox; vBox;1.0; 2*tf], 1, nN);
if freeTauf
    lbX = [lbX; cBox(1)*ones(1,nN)];  ubX = [ubX; cBox(2)*ones(1,nN)];
end
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
if ~isempty(lamG0)
    if numel(lamG0) == size(opti.g,1)
        opti.set_initial(opti.lam_g, lamG0);      % genuine dual warm start (see header)
    else
        warning('casadi_minfuel_sundman:lamG0size', ...
            'opts.lamG0 has %d entries, opti.g has %d rows -- multipliers ignored', ...
            numel(lamG0), size(opti.g,1));
    end
end

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
    % this file's fork had dropped it) so foc_ipopt_inertia / ipopt_certify
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
Dd = Xs(:,2:end) - Xs(:,1:end-1) - tauf*(repmat(dsig,nx,1)/2).*(Fs(:,1:end-1) + Fs(:,2:end));
ss = Us(4,:);

% --- KKT multipliers -> discrete costates + PMP primer-vector check ---------
% The duals of the dynamics-defect constraints ARE the discrete costates
% [lam_r; lam_v; lam_m; lam_t] (one per interval, up to a positive mesh-weight
% scaling and a global sign convention). The primer condition
% alpha* = -lam_v/||lam_v|| is scale-invariant, so comparing it to the NLP
% thrust direction on burn arcs is an independent optimality certificate.
lamDef = [];  primerAlignDeg = NaN;  lamMassEnd = NaN;  lamDef9 = [];
if numel(lamAll) >= nx*N
    lamDef = reshape(lamAll(1:nx*N), nx, N);        % [nx x N] discrete costates
    if freeTauf, lamDef9 = lamDef;  lamDef = lamDef(1:8,:); end   % 8-row contract
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
if freeTauf                                             % cScale (time-scale) box
    slk(end+1) = min(Xs(9,end)-cBox(1), cBox(2)-Xs(9,end));
    lbl{end+1} = 'cBox';
end
[minSlack, iw] = min(slk);
boundSat = struct('minSlack', minSlack, 'worst', lbl{iw}, 'hit', minSlack < 1e-4);
if boundSat.hit
    warning('casadi_minfuel_sundman:boundSaturation', ...
        'nonphysical box ''%s'' within %.2g of binding -- widen via opts before trusting', ...
        lbl{iw}, max(minSlack,0));
end

% 8-row X contract in both modes; the free-tau_f slack is a scalar (constant
% state) and out.tauf becomes the EFFECTIVE regularized length cScale*tauf0.
X9 = [];  cScale = 1;
if freeTauf, X9 = Xs;  cScale = Xs(9,end);  Xs = Xs(1:8,:);  tauf = cScale*tauf0; end
out = struct('X', Xs, 'U', Us, 'tauf', tauf, 'mf', Xs(7,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'switches', sum(abs(diff(ss > 0.5))), ...
             'edge', mean(ss > 0.95 | ss < 0.05), ...
             'lamDef', lamDef, 'lamAll', lamAll, ...
             'primerAlignDeg', primerAlignDeg, 'lamMassEnd', lamMassEnd, ...
             'boundSat', boundSat, ...
             'success', success, 'ipoptStatus', status, 'regHistory', regHistory);

if freeTauf
    out.cScale = cScale;  out.X9 = X9;  out.lamDef9 = lamDef9;  out.taufSeed = tauf0;
end
if returnModel
    out.model = struct('opti', opti, 'creg', creg);
    if freeTauf, out.model.manifest = 'tulip_free'; else, out.model.manifest = 'tulip'; end
end
end
