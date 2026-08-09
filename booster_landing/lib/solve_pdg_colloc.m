function sol = solve_pdg_colloc(P, opts)
% SOLVE_PDG_COLLOC  Min-fuel PDG as a free-tf Hermite-Simpson NLP (IPOPT).
%
% max m(tf) s.t. 3-DOF dynamics, thrust annulus Tmin<=|T|<=Tmax (nonconvex,
% handled directly), glideslope cone, optional pointing cone, r(tf)=0,
% v(tf)=P.vf (ADJUDICATED 2026-08-08: was v(tf)=0, see booster_params.m).
% Time normalized: t = tf*tau, tau in [0,1], tf a decision variable.
%
% ADAPTATION FROM SPEC (documented per task-3 brief): the NLP is built and
% solved in NONDIMENSIONAL variables, not raw SI. A verbatim-SI cold start
% (spec formulation exactly as written) diverges under IPOPT even at
% maxIter=6000 (mass runs past m0, tf leaves its bracket) -- raw SI mixes
% O(1e3) positions, O(1e5) N thrust and O(1e4) kg mass in one Newton step,
% which is textbook-bad Hessian conditioning for a barrier method (see this
% repo's own optimal_control/CLAUDE.md "Scaling Considerations", and Kelly
% 2017 Sec. 2 / Betts 2010). Nondimensionalizing by characteristic scales
% derived from P (Lc, Tc, Mc below; Vc, Fc follow from F=ma consistency)
% converges the identical physics cold in ~30-60 iterations at N=15. Only
% the internal Opti variables are scaled; every field of `sol` is unscaled
% back to SI before it is returned, so the interface below is unchanged.
%
% INPUTS:
%   P    - booster_params struct
%   opts - (optional) .N segments [def P.N], .init prior sol (warm start),
%          .maxIter [def 3000]
% OUTPUTS:
%   sol  - .t .tf .mf .X(7xN+1) .U(3xN+1) .Um(3xN) .lam_defect(7xN)
%          .stats .P   (see plan Task 3 Interfaces; sol.X/.U/.Um/.tf/.t are
%          SI. sol.lam_defect are duals of the NONDIMENSIONAL defect block,
%          extracted via opti.lam_g at rows gDefStart:gDefStart+nDefRows-1
%          -- house sign-bug lesson: never opti.dual.)
%          sol.lam_defect SI conversion (all scales recoverable from sol.P
%          via Lc=norm(P.r0), Tc=0.5*(P.tf_lo+P.tf_hi), Mc=P.m0):
%            rows 1-3 (position defects): lambda_SI = lambda_stored * Mc/Lc
%            rows 4-6 (velocity defects): lambda_SI = lambda_stored * Mc*Tc/Lc  (= Mc/Vc)
%            row  7   (mass defect):      lambda_SI = lambda_stored * 1
%          Rows 4-6 share ONE isotropic scale, so primer-DIRECTION checks
%          against SI thrust are valid on the stored duals unrescaled;
%          cross-block identities (e.g. comparing lambda_r to lambda_v, or
%          lambdadot_v = -lambda_r) are off by a factor of Tc without
%          rescaling.
%
% REFERENCES:
%   [1] Kelly, "An Introduction to Trajectory Optimization," SIAM Rev 2017.
if nargin < 2, opts = struct(); end
if ~isfield(opts,'N'),       opts.N = P.N;        end
if ~isfield(opts,'maxIter'), opts.maxIter = 3000; end
import casadi.*
N    = opts.N;
opti = casadi.Opti();

%% Characteristic scales (nondimensionalization; see ADAPTATION note above):
Lc = norm(P.r0);                  % position scale [m]
Tc = 0.5*(P.tf_lo + P.tf_hi);     % time scale [s]  (bracket midpoint)
Mc = P.m0;                        % mass scale [kg]
Vc = Lc/Tc;                       % velocity scale, F=ma-consistent
Fc = Mc*Lc/Tc^2;                  % force scale,    F=ma-consistent
S  = struct('Lc',Lc,'Tc',Tc,'Mc',Mc,'Vc',Vc,'Fc',Fc, ...
            'ghat', (Tc^2/Lc)*P.gvec, ...          % scaled gravity
            've',   P.Isp*P.g0, ...                % exhaust velocity [m/s]
            'kMdot', Fc*Tc/(Mc*P.Isp*P.g0));        % scaled mass-flow coeff

%% Decision variables (nondimensional):
X  = opti.variable(7, N+1);       % [rhat(3); vhat(3); mhat] at nodes
U  = opti.variable(3, N+1);       % thrust at nodes, /Fc
Um = opti.variable(3, N);         % thrust at HS midpoints, /Fc
tf = opti.variable();             % final time /Tc

%% Dynamics as a CasADi expression (vacuum or drag per P.drag.on):
f = @(x,T) pdg_rhs_casadi(x, T, P, S);

%% Hermite-Simpson defects (register row range for dual extraction):
gDefStart = 1;                    % opti.g rows are appended in order
h = tf / N;
for k = 1:N
    xk = X(:,k);  xk1 = X(:,k+1);
    fk = f(xk, U(:,k));  fk1 = f(xk1, U(:,k+1));
    xm = 0.5*(xk + xk1) + (h/8)*(fk - fk1);
    fm = f(xm, Um(:,k));
    opti.subject_to(xk1 - xk - (h/6)*(fk + 4*fm + fk1) == 0);
end
nDefRows = 7 * N;

%% Path constraints (nodes and midpoints), bounds scaled to match X, U:
Tmin_h = P.Tmin / Fc;  Tmax_h = P.Tmax / Fc;
Uall = [U, Um];
for k = 1:size(Uall,2)
    T2 = sum(Uall(:,k).^2);
    opti.subject_to(Tmin_h^2 <= T2 <= Tmax_h^2);                  % annulus
    if isfinite(P.theta_max_deg)
        opti.subject_to(Uall(3,k)^2 >= cosd(P.theta_max_deg)^2 * T2);
        opti.subject_to(Uall(3,k) >= 0);
    end
end
cotg   = 1 / tand(P.gs_deg);
mdry_h = P.mdry / Mc;
for k = 1:N+1
    opti.subject_to(X(1,k)^2 + X(2,k)^2 <= (cotg * X(3,k))^2);    % glideslope
    opti.subject_to(X(3,k) >= 0);
    opti.subject_to(X(7,k) >= mdry_h);
end

%% Boundary conditions + tf bounds (scaled):
% Terminal velocity ADJUDICATED 2026-08-08 (task-7 fix report round 3):
% X(4:6,end) == vf_h, was zeros(3,1) -- see P.vf's comment in
% booster_params.m for why (Tmin>weight makes v(tf)=0 singular).
r0_h = P.r0 / Lc;  v0_h = P.v0 / Vc;  vf_h = P.vf / Vc;  m0_h = P.m0 / Mc;
opti.subject_to(X(1:3,1) == r0_h);
opti.subject_to(X(4:6,1) == v0_h);
opti.subject_to(X(7,1)   == m0_h);
opti.subject_to(X(1:3,end) == zeros(3,1));
opti.subject_to(X(4:6,end) == vf_h);
opti.subject_to(P.tf_lo/Tc <= tf <= P.tf_hi/Tc);

%% Objective: min fuel == max final mass:
opti.minimize(-X(7,end));

%% Initial guess: straight line, gravity-cancelling thrust, or warm start
%% (built in SI per spec, then divided down into the scaled variables):
if isfield(opts, 'init') && ~isempty(opts.init)
    s0 = opts.init;
    tauN = linspace(0,1,N+1);  tauM = (tauN(1:end-1)+tauN(2:end))/2;
    Xg_si = interp1(s0.t/s0.tf, s0.X.',  tauN, 'pchip').';
    Ug_si = interp1(s0.t/s0.tf, s0.U.',  tauN, 'pchip').';
    Um_si = interp1(s0.t/s0.tf, s0.U.',  tauM, 'pchip').';
    opti.set_initial(X,  [Xg_si(1:3,:)/Lc; Xg_si(4:6,:)/Vc; Xg_si(7,:)/Mc]);
    opti.set_initial(U,  Ug_si/Fc);
    opti.set_initial(Um, Um_si/Fc);
    opti.set_initial(tf, s0.tf/Tc);
else
    tauN  = linspace(0,1,N+1);
    Xg    = [P.r0*(1-tauN); P.v0*(1-tauN); ...
             P.m0 + (P.mdry + 500 - P.m0)*tauN];
    Tg    = repmat(-0.9*P.m0*P.gvec, 1, N+1);      % ~hover thrust, straight up
    opti.set_initial(X, [Xg(1:3,:)/Lc; Xg(4:6,:)/Vc; Xg(7,:)/Mc]);
    opti.set_initial(U, Tg/Fc);
    opti.set_initial(Um, Tg(:,1:N)/Fc);
    opti.set_initial(tf, 30/Tc);
end

%% Solve:
sopts = struct('ipopt', struct('max_iter', opts.maxIter, ...
               'print_level', 3, 'tol', 1e-9));
opti.solver('ipopt', struct('print_time', false), sopts.ipopt);
try
    osol = opti.solve();
    ok   = true;
catch
    osol = opti.debug;                 % return best iterate for diagnosis
    ok   = false;
end

%% Package (unscale back to SI):
tf_h     = full(osol.value(tf));
Xh       = full(osol.value(X));
Uh       = full(osol.value(U));
Umh      = full(osol.value(Um));
sol.tf   = tf_h * Tc;
sol.t    = linspace(0, sol.tf, N+1);
sol.X    = [Xh(1:3,:)*Lc; Xh(4:6,:)*Vc; Xh(7,:)*Mc];
sol.U    = Uh * Fc;
sol.Um   = Umh * Fc;
sol.mf   = sol.X(7,end);
lam      = full(osol.value(opti.lam_g));
sol.lam_defect = reshape(lam(gDefStart : gDefStart+nDefRows-1), 7, N);
sol.stats = struct('success', ok, 'status', opti.stats.return_status, ...
                   'iter', opti.stats.iter_count);
sol.P    = P;
end

function xdot = pdg_rhs_casadi(x, T, P, S)
% PDG_RHS_CASADI  Same RHS as pdg_dynamics, CasADi-symbolic-safe, evaluated
% in the NONDIMENSIONAL variables of the enclosing solve (scale struct S).
% Gate G2 (certify_pdg) proves this matches pdg_dynamics (SI) after the
% solution is unscaled and re-integrated with ode45.
v    = x(4:6);  m = x(7);                    % nondimensional v, m
Tmag = sqrt(sum(T.^2) + 1e-12);              % smooth at |T|=0 (never active)
aD   = 0*v;                                  % vacuum default (symbolic-safe)
if P.drag.on
    Hhat  = P.drag.H / S.Lc;
    kDhat = 0.5*P.drag.rho0*P.drag.Cd*P.drag.A*S.Vc*S.Tc/S.Mc;
    vmag  = sqrt(sum(v.^2) + 1e-12);
    aD    = -kDhat * exp(-x(3)/Hhat) * vmag * v / m;
end
xdot = [v; S.ghat + T/m + aD; -S.kMdot*Tmag];
end
