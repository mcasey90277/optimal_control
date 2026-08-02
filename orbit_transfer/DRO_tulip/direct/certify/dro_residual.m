function R = dro_residual(o, muStar, Tmax, c, opts)
% DRO_RESIDUAL  True continuous-time local error of a direct DRO->tulip solution.
%
% THE POINT OF THIS FUNCTION. A collocation NLP drives its DEFECTS to machine
% precision, but a defect only says the returned numbers satisfy the QUADRATURE
% RULE. It says nothing about whether that rule approximates the dynamics. On
% this problem the two differ by up to seven orders of magnitude, and the gap is
% what let a trapezoidal solve return a trajectory that beat the indirect
% reference minimum time by flying a lunar pass it could not actually fly.
%
% So: for each interval, start at the LEFT node, reconstruct the control the
% transcription implies over that interval, integrate the TRUE equations of
% motion at tight tolerance, and compare against the right node.
%
%   R_k = || Phi(t_{k+1}; t_k, X_k, U(.)) - X_{k+1} ||
%
% This is the local error the scheme is committing. It is the number that decides
% whether a direct solution is physical.
%
% THE CONTROL BETWEEN NODES IS A CHOSEN CONVENTION, NOT A CONSEQUENCE. This
% matters and the earlier version of this comment got it wrong. The transcription
% constrains the dynamics only AT the collocation samples; it says nothing about
% the control between them, so there is no such thing as "the" intersample
% control. What we reconstruct here is the convention CONSISTENT WITH each
% scheme's derivation: linear for trapezoid, and for Hermite-Simpson the
% Lagrange quadratic through the two node controls and the midpoint control,
% since Simpson's rule is derived by assuming the integrand is quadratic through
% those three samples. A different convention would give a different -- also
% defensible -- residual. Report the convention alongside the number.
%
% Two artifacts of the convention are measured and returned rather than hidden:
% the reconstructed direction can shrink toward zero between unit samples
% (.minDirNorm), and the reconstructed throttle can overshoot [0,1]
% (.thrOvershoot). Both are silently repaired inside the RHS; large values mean
% the reconstruction, not the solution, is what is being measured.
%
% INPUTS:
%   o      - solution struct from casadi_mintime_dro. Uses .X [7x(N+1)],
%            .U [4x(N+1)], .s [1x(N+1)], .tf, and, when
%            o.scheme = 'hermite-simpson', .Um [4xN] midpoint controls.
%   muStar - CR3BP mass ratio [scalar]
%   Tmax   - ND thrust acceleration at unit mass fraction [scalar]
%   c      - ND exhaust velocity [scalar]
%   opts   - struct (optional):
%            .relTol [1e-11]  .absTol [1e-13]   ode113 tolerances. These must
%                     stay well below the residuals being measured, or the
%                     measurement reports its own integrator error.
%
% OUTPUTS:
%   R - struct: .Rr/.RrMax/.RrMed  POSITION error [ND length; x lStar for km]
%       .Rv/.RvMax/.RvMed  VELOCITY error [ND; x lStar/tStar for km/s]
%       .Rm/.RmMax         mass-fraction error [dimensionless]
%       .Rx [1xN] the combined 7-component norm, retained for continuity but
%       DIMENSIONALLY MEANINGLESS -- do not convert it to km. .RxMax .RxMed,
%       .kWorst (index of the worst interval), .tMid [1xN] interval midtimes,
%       .scheme (the reconstruction convention used)
%       .minDirNorm  smallest reconstructed direction norm (1 = healthy; near 0
%                    means the quadratic passed close to the origin and the
%                    reconstructed direction is meaningless there)
%       .thrOvershoot largest excursion of the reconstructed throttle outside
%                    [0,1] (0 = healthy)
%       .relTol .absTol
%
% REFERENCES:
%   [1] Betts, "Practical Methods for Optimal Control", Ch. 4 -- discretization
%       error estimation and mesh refinement.
%   [2] orbit_transfer/DRO_tulip/FINDINGS.md -- the measurement this implements.

if nargin < 5, opts = struct(); end
relTol = local_default(opts, 'relTol', 1e-11);
absTol = local_default(opts, 'absTol', 1e-13);
scheme = 'trapezoid';
if isfield(o,'scheme') && ~isempty(o.scheme), scheme = lower(o.scheme); end
hasMid = strncmp(scheme,'hermite',7) && isfield(o,'Um') && ~isempty(o.Um);

% Node times: uniform s*tf normally; the solver's own PHYSICAL times under
% Sundman, where the mesh is deliberately non-uniform in t.
if isfield(o,'tNodes') && ~isempty(o.tNodes)
    t = o.tNodes(:).';
else
    t = o.s(:).' * o.tf;
end
N  = numel(t) - 1;
Rx = nan(1,N);  Rr = nan(1,N);  Rv = nan(1,N);  Rm = nan(1,N);
minDir = inf;  thrOver = 0;
odeo = odeset('RelTol',relTol,'AbsTol',absTol);

for k = 1:N
    a = t(k);  b = t(k+1);  h = max(b-a, eps);
    if hasMid
        % quadratic through u(a), u(mid), u(b) -- the Hermite-Simpson control
        ua = o.U(:,k);  um = o.Um(:,k);  ub = o.U(:,k+1);
        uOf = @(tt) local_quad(ua, um, ub, (tt-a)/h);
    else
        ua = o.U(:,k);  ub = o.U(:,k+1);
        uOf = @(tt) ua + ((tt-a)/h)*(ub-ua);
    end
    % reconstruction health, sampled across the interval
    for w = [0.25 0.5 0.75]
        uw = uOf(a + w*h);
        minDir = min(minDir, norm(uw(1:3)));
        thrOver = max(thrOver, max(0, max(uw(4)-1, -uw(4))));
    end
    [~, Z] = ode113(@(tt,z) local_rhs(z, uOf(tt), muStar, Tmax, c), [a b], o.X(:,k), odeo);
    e = Z(end,:).' - o.X(:,k+1);
    Rx(k) = norm(e);            % mixed norm -- dimensionally meaningless, see below
    Rr(k) = norm(e(1:3));       % POSITION only, convertible to km
    Rv(k) = norm(e(4:6));       % VELOCITY only, convertible to km/s
    Rm(k) = abs(e(7));          % mass fraction
end

% SEPARATE THE COMPONENTS, and this is not pedantry. The original version of
% this function returned only the 7-component norm and the caller multiplied it
% by lStar to get "kilometres". That is dimensionally invalid: it converts a
% quantity built from position (ND length), velocity (ND length/time) and mass
% fraction (dimensionless) as though all three were lengths. External review
% caught it. On the certified N = 1600 solution the mixed norm was dominated by
% the VELOCITY component, so the headline "0.32 km" was not a position error at
% all. Callers should use Rr/Rv/Rm and convert each with its own scale.
[RxMax, kWorst] = max(Rx);
[RrMax, kWorstR] = max(Rr);
R = struct('Rx',Rx, 'RxMax',RxMax, 'RxMed',median(Rx), 'kWorst',kWorst, ...
           'Rr',Rr, 'RrMax',RrMax, 'RrMed',median(Rr), 'kWorstR',kWorstR, ...
           'Rv',Rv, 'RvMax',max(Rv), 'RvMed',median(Rv), ...
           'Rm',Rm, 'RmMax',max(Rm), ...
           'tMid',0.5*(t(1:end-1)+t(2:end)), 'scheme',scheme, ...
           'minDirNorm',minDir, 'thrOvershoot',thrOver, ...
           'relTol',relTol, 'absTol',absTol);
end

% ---------------------------------------------------------------------------
function u = local_quad(ua, um, ub, w)
% LOCAL_QUAD  Lagrange quadratic through ua at w=0, um at w=1/2, ub at w=1.
% INPUTS: ua, um, ub [4x1]; w (normalized position in [0,1])
% OUTPUTS: u [4x1]
u = (2*w-1).*(w-1).*ua + 4*w.*(1-w).*um + w.*(2*w-1).*ub;
end

% ---------------------------------------------------------------------------
function dz = local_rhs(z, u, mu, Tmax, c)
% LOCAL_RHS  CR3BP with thrust. Mirrors pumpkyn.cr3bp.tfMinEoM and the dynamics
% built inside casadi_mintime_dro -- if these two ever diverge the residual
% becomes meaningless, so they are deliberately written the same way.
% INPUTS: z [7x1]; u [4x1]; mu; Tmax; c   OUTPUTS: dz [7x1]
r = z(1:3);  v = z(4:6);  m = z(7);
al = u(1:3);  al = al/max(norm(al),eps);  th = min(max(u(4),0),1);
dd = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
gr = [r(1) - (1-mu)*(r(1)+mu)/dd^3 - mu*(r(1)-1+mu)/rr^3;
      r(2) - (1-mu)*r(2)/dd^3      - mu*r(2)/rr^3;
           - (1-mu)*r(3)/dd^3      - mu*r(3)/rr^3];
dz = [v; gr + [2*v(2); -2*v(1); 0] + (th*Tmax/m)*al; -(Tmax/c)*th];
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
