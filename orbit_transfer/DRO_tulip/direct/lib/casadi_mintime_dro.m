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
% NO SUNDMAN REGULARIZATION, and that is deliberate. The GTO->tulip campaign
% needs it because a ~40-revolution Earth spiral passes deep perigee. Measured
% on this transfer: 1.18 revolutions about the Moon, Earth distance never below
% 0.85, closest lunar approach 0.0164 (~4650 km altitude). Plain time-domain
% collocation is appropriate; adding a regularization would be machinery
% imported for its own sake.
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
%         .minAltKm (the constraint applied, NaN if none)
%         .altMinKm (the ACHIEVED minimum altitude at the nodes)
%         .altActive (true if the constraint is within 1 km of binding)
%         .ipoptStatus .maxDefect .maxUnit .termErr .thrMin (lowest throttle --
%         the saturation check) .mf .lamDef [7xN] defect multipliers
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
r0 = size(opti.g,1)+1;
D = X(:,2:end) - X(:,1:end-1) ...
    - repmat(TF(1:end-1),7,1).*(repmat(ds,7,1)/2).*(Fv(:,1:end-1) + Fv(:,2:end));
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

opti.subject_to(X(1:6,1) == rv0(:));
opti.subject_to(X(7,1) == 1);
opti.subject_to(X(1:6,end) == rvf(:));
opti.subject_to(TF(:) > 0);
opti.subject_to(X(7,:).' > 0.05);          % mass stays physical

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
    if retModel, creg(end+1) = struct('label','minAlt','rows',r0:size(opti.g,1)); end
end

opti.minimize(TF(1));

% --- warm start -------------------------------------------------------------
if isempty(X0)
    % Crude: linear in state, mass ramped, thrust along the chord. Deliberately
    % NOT the indirect answer -- if this converges to the same t_f, the two
    % methods agree independently rather than by construction.
    X0 = zeros(7,nN);
    for k = 1:6, X0(k,:) = linspace(rv0(k), rvf(k), nN); end
    X0(7,:) = linspace(1, 0.92, nN);
    ch = rvf(1:3) - rv0(1:3);  ch = ch/max(norm(ch),eps);
    U0 = [repmat(ch(:),1,nN); ones(1,nN)];
    if isempty(tf0), tf0 = 4.0; end
end
opti.set_initial(X, X0);
opti.set_initial(U, U0);
opti.set_initial(TF, tf0*ones(1,nN));

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

% diagnostics
Fn = zeros(7,nN);
for k = 1:nN, Fn(:,k) = full(fdyn(out.X(:,k), out.U(:,k))); end
Dn = out.X(:,2:end) - out.X(:,1:end-1) - out.tf*(repmat(ds,7,1)/2).*(Fn(:,1:end-1)+Fn(:,2:end));
out.maxDefect = max(abs(Dn(:)));
out.maxUnit   = max(abs(vecnorm(out.U(1:3,:),2,1) - 1));
out.termErr   = norm(out.X(1:6,end) - rvf(:));
out.thrMin    = min(out.U(4,:));
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
    out.lamDef = reshape(lamAll(creg(1).rows), 7, N);
catch
    out.lamDef = [];
end
if retModel, out.model = struct('opti',opti,'creg',creg,'X',X,'U',U,'TF',TF); end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
