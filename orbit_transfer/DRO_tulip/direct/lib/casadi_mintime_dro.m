function out = casadi_mintime_dro(rv0, rvf, Tmax, c, muStar, N, X0, U0, tf0, opts)
% CASADI_MINTIME_DRO  Direct-collocation minimum-time CR3BP transfer, free t_f.
%
% The direct twin of pumpkyn.cr3bp.tfMin. Same problem, same dynamics, same
% endpoints -- solved by transcription instead of by shooting, so the two can be
% compared. That comparison is the point: DRO->tulip is the only problem in the
% catalog where an independently converged INDIRECT solution already exists, so
% it is the one place the direct method's dual-to-costate mapping can be checked
% against a second opinion rather than against itself.
%
% NO SUNDMAN REGULARIZATION -- originally justified by the SHAPE of the transfer
% (1.18 revolutions about the Moon, closest approach ~4650 km) rather than by
% measurement, which was the wrong way to decide it. Measuring the continuous
% error afterwards (certify/dro_residual) showed the second-order scheme does not
% resolve this transfer at any density tried.
%
% The remedy applied was the fourth-order scheme below (opts.scheme), and it was
% sufficient: 3.3 m worst-interval POSITION error at N = 1600. A Sundman change
% of variable and a periselene-concentrated mesh remain untried, and are now
% OPTIONAL rather than necessary. Do not restore the old claim: see FINDINGS.md.
%
% FREE FINAL TIME BY NORMALIZED TIME, WITH t_f LIFTED. The independent variable
% is s in [0,1] and dX/ds = t_f * f(X,U). A single scalar t_f couples to every
% defect, giving one DENSE KKT column -- and MUMPS does not survive it: the
% first version of this file used a scalar and MATLAB died with a fatal error
% immediately after IPOPT printed the Jacobian structure, at only N = 400. That
% is the same failure the tulip solver avoids by fixing tau_f outright, and the
% same one the earth CR3BP campaign cured with liftDL.
%
% The cure here is that lift: t_f is replicated as a per-node variable with
% local continuity T(k+1) = T(k), so each defect touches only its own two nodes
% and the arrowhead becomes banded. The solution is unchanged -- the continuity
% constraints force every T(k) equal -- but the sparsity is what MUMPS needs.
%
% THE THROTTLE IS LEFT FREE. For minimum time the optimum is thrust at its upper
% bound, and pinning it would be simpler -- but leaving it free and CHECKING
% that it saturates is a genuine test of the formulation. Note the reference
% implementation also permits a switch: tfMinEoM sets u = 0 when
% S = -||lambda_v||*c/m - lambda_m > 0. In the converged solution lambda_m runs
% 5.50 -> 0 while staying positive, so S < 0 and u = 1 throughout.
%
% INPUTS:
%   rv0, rvf - departure / arrival states [1x6, ND rotating barycentric]
%   Tmax     - ND thrust acceleration at unit mass fraction [scalar]
%   c        - ND exhaust velocity [scalar]
%   muStar   - mass ratio [scalar]
%   N        - number of collocation intervals [scalar]
%   X0       - warm-start states [7x(N+1)] = [r;v;m]; [] for an internal guess
%   U0       - warm-start controls [4x(N+1)] = [alpha;thr]; [] for internal
%   tf0      - warm-start final time [scalar]
%   opts     - struct (optional):
%              .maxIter [3000]  .printLevel [0]  .returnModel [false]
%              .thrLock [false] true pins thr == 1 (the classical min-time form)
%              .scheme  ['trapezoid'] or 'hermite-simpson'. Trapezoid is second
%                       order and is the historical default, kept so the record
%                       in FINDINGS.md still reproduces. Hermite-Simpson is
%                       fourth order at the same NODE count -- but NOT the same
%                       constraint count: the separated form below adds 7N
%                       interpolation equalities and 4N midpoint controls with
%                       their own norm and bound constraints. It is what the
%                       accuracy gate (certify/certify_dro_mintime) requires.
%              .minAltKm [] minimum altitude above the LUNAR SURFACE, km. Empty
%                       (default) leaves the problem UNCONSTRAINED, which
%                       reproduces the original ill-posed sweep. Any finite
%                       value adds the path constraint below.
%              .lStarKm [389703.264829278] characteristic length, for the km
%                       conversion
%              .rMoonKm [1737.4] lunar radius
%
% OUTPUTS:
%   out - struct: .X [7x(N+1)] .U [4x(N+1)] .tf .s (the grid) .success
%         .scheme, and .Um [4xN] midpoint controls (hermite-simpson only, else [])
%         .minAltKm (the constraint applied, NaN if none)
%         .altMinKm (the ACHIEVED minimum altitude at the nodes)
%         .altActive (true if the constraint is within 1 km of binding)
%         .ipoptStatus .maxDefect .maxUnit .termErr .thrMin (lowest throttle --
%         the saturation check) .mf .lamDef [7xN] defect multipliers
%         .lamBCf [6x1] TERMINAL boundary multipliers -- the Hager terminal
%         covector, i.e. lambda(t_f) up to sign; .lamBC0 [7x1] initial ones.
%         All three require opts.returnModel = true.
%         .tfSpread (max deviation across the lifted copies -- should be ~0)
%         .model (opts.returnModel only)
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMinEoM (the dynamics and control law mirrored here).
%   [2] orbit_transfer/doc/transfer_problem_space.md.

if nargin < 10, opts = struct(); end
g = @(f,v) local_default(opts, f, v);
maxIter  = g('maxIter', 3000);
prnt     = g('printLevel', 0);
retModel = g('returnModel', false);
thrLock  = g('thrLock', false);
scheme   = lower(g('scheme', 'trapezoid'));
minAltKm = g('minAltKm', []);
lStarKm  = g('lStarKm', 389703.264829278);
rMoonKm  = g('rMoonKm', 1737.4);

import casadi.*
nN = N + 1;
s  = linspace(0, 1, nN);
ds = diff(s);

% --- dynamics, mirroring tfMinEoM ------------------------------------------
x = MX.sym('x',7);  u = MX.sym('u',4);
r = x(1:3);  v = x(4:6);  m = x(7);  al = u(1:3);  th = u(4);
dd = sqrt((r(1)+muStar)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+muStar)^2 + r(2)^2 + r(3)^2 + 1e-12);
d3 = dd^3;  r3 = rr^3;
gr = [r(1) - (1-muStar)*(r(1)+muStar)/d3 - muStar*(r(1)-1+muStar)/r3;
      r(2) - (1-muStar)*r(2)/d3          - muStar*r(2)/r3;
           - (1-muStar)*r(3)/d3          - muStar*r(3)/r3];
hv = [2*v(2); -2*v(1); 0];
acc = gr + hv + (th*Tmax/m)*al;
fdyn = Function('f', {x,u}, {[v; acc; -(Tmax/c)*th]});

% --- NLP --------------------------------------------------------------------
opti = Opti();
X  = opti.variable(7, nN);
U  = opti.variable(4, nN);
TF = opti.variable(1, nN);            % lifted final time, one copy per node
creg = struct('label',{},'rows',{});

F = fdyn.map(nN);
Fv = F(X, U);
Um = [];  Xm = [];
switch scheme
case 'trapezoid'
    % Second order. dX/ds = tf*f, integrated by the trapezoid rule. The control
    % is implicitly linear across each interval.
    D = X(:,2:end) - X(:,1:end-1) ...
        - repmat(TF(1:end-1),7,1).*(repmat(ds,7,1)/2).*(Fv(:,1:end-1) + Fv(:,2:end));

case 'hermite-simpson'
    % FOURTH order, in the SEPARATED form: the interior midpoint state Xm is a
    % decision variable constrained to the Hermite interpolant, rather than an
    % expression substituted into the dynamics.
    %
    %   interpolation:  Xm - (Xk + Xk1)/2 - (h/8)(fk - fk1) = 0
    %   defect:         Xk1 - Xk - (h/6)(fk + 4 f(Xm,Um) + fk1) = 0
    %
    % WHY SEPARATED AND NOT COMPRESSED. The compressed form -- substituting the
    % interpolant directly into f -- is algebraically identical and carries half
    % the equations, and it was tried first. It solves at N = 400 and MATLAB
    % dies SILENTLY at N = 800: nesting Xm(X,TF) inside f means the Lagrangian
    % Hessian needs the chain rule through the interpolant, the expression graph
    % deepens with every node, and MUMPS does not survive it. That is the same
    % failure, and the same cure, as the lift of t_f above and as liftDL in the
    % earth CR3BP campaign: keep every constraint SHALLOW and let extra
    % variables carry the coupling. Do not "simplify" this back.
    %
    % The midpoint CONTROL is likewise a genuine variable, not the average of
    % its neighbours. Averaging would tie the control to a piecewise-linear
    % family and cost the scheme its fourth order; and here the control lives on
    % a sphere, where the average of two unit vectors is not a unit vector. Um
    % carries its own norm and throttle constraints below.
    hND = repmat(TF(1:end-1),7,1).*repmat(ds,7,1);       % 7 x N, the ND step
    Um  = opti.variable(4, N);
    Xm  = opti.variable(7, N);
    Fm  = fdyn.map(N);
    Fmv = Fm(Xm, Um);
    rI = size(opti.g,1)+1;
    % reshape needs EXPLICIT dimensions here -- CasADi has no MX overload for
    % the MATLAB '[]' placeholder.
    HSI = Xm - 0.5*(X(:,1:end-1) + X(:,2:end)) - (hND/8).*(Fv(:,1:end-1) - Fv(:,2:end));
    opti.subject_to(reshape(HSI, 7*N, 1) == 0);
    if retModel, creg(end+1) = struct('label','hsInterp','rows',rI:size(opti.g,1)); end
    D = X(:,2:end) - X(:,1:end-1) - (hND/6).*(Fv(:,1:end-1) + 4*Fmv + Fv(:,2:end));

otherwise
    error('casadi_mintime_dro:scheme', ...
        'unknown scheme ''%s''; use ''trapezoid'' or ''hermite-simpson''.', scheme);
end
r0 = size(opti.g,1)+1;
opti.subject_to(D(:) == 0);
if retModel, creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1)); end

r0 = size(opti.g,1)+1;
opti.subject_to((TF(2:end) - TF(1:end-1)).' == 0);     % the lift's continuity
if retModel, creg(end+1) = struct('label','tfCont','rows',r0:size(opti.g,1)); end

r0 = size(opti.g,1)+1;
opti.subject_to((sum(U(1:3,:).^2,1) - 1).' == 0);
if retModel, creg(end+1) = struct('label','betaNorm','rows',r0:size(opti.g,1)); end

if thrLock
    opti.subject_to(U(4,:).' == 1);
else
    opti.subject_to(U(4,:).' >= 0);
    opti.subject_to(U(4,:).' <= 1);
end

if ~isempty(Um)
    % The midpoint control gets the same treatment as the node controls --
    % otherwise it is free to leave the unit sphere and the interpolated
    % direction stops meaning anything.
    r0 = size(opti.g,1)+1;
    opti.subject_to((sum(Um(1:3,:).^2,1) - 1).' == 0);
    if retModel, creg(end+1) = struct('label','betaNormMid','rows',r0:size(opti.g,1)); end
    if thrLock
        opti.subject_to(Um(4,:).' == 1);
    else
        opti.subject_to(Um(4,:).' >= 0);
        opti.subject_to(Um(4,:).' <= 1);
    end
end

r0 = size(opti.g,1)+1;
opti.subject_to(X(1:6,1) == rv0(:));
opti.subject_to(X(7,1) == 1);
if retModel, creg(end+1) = struct('label','bc0','rows',r0:size(opti.g,1)); end

% The TERMINAL boundary rows are registered separately because their multipliers
% are the Hager terminal covector: stationarity of the Lagrangian with respect to
% X(:,end) reads  nu_N'*dD_N/dX_end + nu_psi'*dpsi/dX_end = 0, and since
% dD_N/dX_end = I + O(h), the terminal costate is -nu_psi to leading order. That
% is a cleaner reading of lambda(t_f) than extrapolating the interval
% multipliers, and it is what foc_check uses elsewhere in this repository.
r0 = size(opti.g,1)+1;
opti.subject_to(X(1:6,end) == rvf(:));
if retModel, creg(end+1) = struct('label','bcf','rows',r0:size(opti.g,1)); end
opti.subject_to(TF(:) >= 1e-3);   % non-strict: an NLP cannot impose '>'
opti.subject_to(X(7,:).' >= 0.05);         % mass stays physical (non-strict)

% --- minimum-altitude path constraint ---------------------------------------
% WITHOUT THIS THE PROBLEM IS ILL-POSED. Nothing else in the formulation stops
% the trajectory approaching the Moon, and a deeper flyby buys a stronger
% gravity assist, so inf t_f is approached as the trajectory grazes the lunar
% surface -- 'the minimum time' is a family parameterized by flyby altitude, not
% a number. Measured: an unconstrained solve at N = 800 returned t_f BELOW the
% indirect reference by flying a 442 km pass against the reference's 4674 km.
%
% Imposed SQUARED, which avoids a square root whose derivative blows up as the
% radius shrinks -- precisely where the constraint is active.
%
% NOTE WHAT THIS DOES AND DOES NOT GUARANTEE. Like every constraint in a
% collocation transcription it binds AT THE NODES ONLY. Between nodes the
% trajectory is never evaluated, so a mesh too coarse to resolve the periselene
% can still permit an excursion below the floor. Verify with a high-accuracy
% propagation of the reconstructed control before trusting a tight case.
rhoMinND = NaN;
if ~isempty(minAltKm)
    rhoMinND = (rMoonKm + minAltKm) / lStarKm;
    dMoon2 = sum((X(1:3,:) - [1-muStar; 0; 0]).^2, 1).';
    r0 = size(opti.g,1)+1;
    opti.subject_to(dMoon2 >= rhoMinND^2);
    % Under Hermite-Simpson the interior midpoint state is a computable
    % expression, so the floor can be enforced there too. That halves the
    % unguarded span without adding a variable -- it does not eliminate it.
    if ~isempty(Xm)
        dMid2 = sum((Xm(1:3,:) - [1-muStar; 0; 0]).^2, 1).';
        opti.subject_to(dMid2 >= rhoMinND^2);
    end
    if retModel, creg(end+1) = struct('label','minAlt','rows',r0:size(opti.g,1)); end
end

opti.minimize(TF(1));

% --- warm start -------------------------------------------------------------
% Seed the state and the control INDEPENDENTLY. The earlier version built both
% inside 'if isempty(X0)', so a caller supplying U0 but not X0 had its U0
% silently discarded, and a caller supplying X0 but not U0 reached
% set_initial(U,[]) and failed.
if isempty(X0)
    % Crude: linear in state, mass ramped. Deliberately NOT the indirect answer
    % -- if this converges to the same t_f, the two methods agree independently
    % rather than by construction.
    X0 = zeros(7,nN);
    for k = 1:6, X0(k,:) = linspace(rv0(k), rvf(k), nN); end
    X0(7,:) = linspace(1, 0.92, nN);
end
if isempty(U0)
    ch = rvf(1:3) - rv0(1:3);  ch = ch/max(norm(ch),eps);
    if ~(norm(ch) > 0.5), ch = [1;0;0]; end      % degenerate endpoints
    U0 = [repmat(ch(:),1,nN); ones(1,nN)];
end
if isempty(tf0), tf0 = 4.0; end   % also covers a caller who supplies X0 but not tf0
opti.set_initial(X, X0);
opti.set_initial(U, U0);
opti.set_initial(TF, tf0*ones(1,nN));
if ~isempty(Xm)
    % Seed Xm with the actual HERMITE INTERPOLANT of the seed, not the plain
    % midpoint average. The average violates the interpolation constraint, so
    % IPOPT begins in restoration and walks out of the seed's basin before it
    % starts optimizing -- measured: seeded from the indirect reference at
    % t_f = 4.01524, the average-seeded solve converged to 4.68089 (+16.6%)
    % while the compressed form, whose Xm was consistent by construction,
    % stayed at 4.01734. A warm start must be feasible for the constraints that
    % define it, or it is not a warm start.
    F0 = full(F(X0, U0));            % F is fdyn.map(nN); no MATLAB-level loop
    h0 = tf0*repmat(ds,7,1);
    opti.set_initial(Xm, 0.5*(X0(:,1:end-1) + X0(:,2:end)) ...
                         + (h0/8).*(F0(:,1:end-1) - F0(:,2:end)));
end
if ~isempty(Um)
    Um0 = 0.5*(U0(:,1:end-1) + U0(:,2:end));
    nrm = vecnorm(Um0(1:3,:),2,1);
    % Averaging two unit vectors does not give a unit vector, and averaging two
    % OPPOSED ones gives zero -- which would leave the midpoint direction
    % violating its own norm constraint at the start. Fall back to the left
    % node's direction in that case rather than dividing by eps.
    bad = nrm < 1e-8;
    Um0(1:3,bad) = U0(1:3,find(bad));
    nrm(bad) = 1;
    Um0(1:3,:) = Um0(1:3,:) ./ max(nrm, eps);
    opti.set_initial(Um, Um0);
end

opti.solver('ipopt', struct('print_time',false), ...
    struct('print_level',prnt,'max_iter',maxIter,'tol',1e-10, ...
           'acceptable_tol',1e-8,'linear_solver','mumps'));

out = struct();
try
    sol = opti.solve();
    ok = true;
catch
    sol = opti.debug;  ok = false;
end
out.X  = full(sol.value(X));
out.U  = full(sol.value(U));
out.tf = full(sol.value(TF(1)));
out.tfSpread = max(abs(full(sol.value(TF)) - out.tf));
out.s  = s;
out.ipoptStatus = opti.stats().return_status;
out.success = ok && any(strcmp(out.ipoptStatus, ...
    {'Solve_Succeeded','Solved_To_Acceptable_Level'}));

% diagnostics. The defect is recomputed from the RETURNED numbers using the
% same rule the NLP imposed -- a trapezoid check on a Hermite-Simpson solution
% would report the difference between the two schemes, not a solver error.
out.scheme = scheme;
out.Um = [];  out.Xm = [];
if ~isempty(Um), out.Um = full(sol.value(Um)); end
if ~isempty(Xm), out.Xm = full(sol.value(Xm)); end
Fn = zeros(7,nN);
for k = 1:nN, Fn(:,k) = full(fdyn(out.X(:,k), out.U(:,k))); end
hn = out.tf*repmat(ds,7,1);
if isempty(out.Um)
    Dn = out.X(:,2:end) - out.X(:,1:end-1) - (hn/2).*(Fn(:,1:end-1)+Fn(:,2:end));
else
    % Use the SOLVER'S Xm, not a rebuilt interpolant. Rebuilding it algebraically
    % would satisfy the interpolation constraint by construction and so would
    % hide any violation of it -- the check would be checking itself.
    Xmn = out.Xm;
    Fmap = fdyn.map(N);
    Fmn = full(Fmap(Xmn, out.Um));
    Dn = out.X(:,2:end) - out.X(:,1:end-1) - (hn/6).*(Fn(:,1:end-1) + 4*Fmn + Fn(:,2:end));
    % and report the interpolation residual separately
    XmRef = 0.5*(out.X(:,1:end-1)+out.X(:,2:end)) + (hn/8).*(Fn(:,1:end-1)-Fn(:,2:end));
    out.maxInterp = max(abs(XmRef(:) - out.Xm(:)));
end
out.maxDefect = max(abs(Dn(:)));
if ~isfield(out,'maxInterp'), out.maxInterp = 0; end
out.maxUnit   = max(abs(vecnorm(out.U(1:3,:),2,1) - 1));
if ~isempty(out.Um)
    out.maxUnit = max(out.maxUnit, max(abs(vecnorm(out.Um(1:3,:),2,1) - 1)));
end
out.termErr   = norm(out.X(1:6,end) - rvf(:));
out.thrMin    = min(out.U(4,:));
if ~isempty(out.Um), out.thrMin = min(out.thrMin, min(out.Um(4,:))); end
r2n = vecnorm(out.X(1:3,:) - [1-muStar;0;0], 2, 1);
out.altMinKm  = min(r2n)*lStarKm - rMoonKm;
out.minAltKm  = NaN;  out.altActive = false;
if ~isempty(minAltKm)
    out.minAltKm  = minAltKm;
    out.altActive = (out.altMinKm - minAltKm) < 1;   % within 1 km of binding
end
out.mf        = out.X(7,end);
try
    lamAll = full(sol.value(opti.lam_g));
    kd = find(strcmp({creg.label},'defect'), 1);   % NOT creg(1): Hermite-Simpson
    out.lamDef = reshape(lamAll(creg(kd).rows), 7, N);   % registers interp first
    kf = find(strcmp({creg.label},'bcf'), 1);
    if ~isempty(kf), out.lamBCf = lamAll(creg(kf).rows); else, out.lamBCf = []; end
    k0 = find(strcmp({creg.label},'bc0'), 1);
    if ~isempty(k0), out.lamBC0 = lamAll(creg(k0).rows); else, out.lamBC0 = []; end
catch
    out.lamDef = [];  out.lamBCf = [];  out.lamBC0 = [];
end
if retModel, out.model = struct('opti',opti,'creg',creg,'X',X,'U',U,'TF',TF); end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
