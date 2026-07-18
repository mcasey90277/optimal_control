function [sigma, X0, U0, dL0, info] = mee_seed(par, opts)
% MEE_SEED  Constant-throttle warm-start seed for the L-domain MEE collocation
% solver (Task 3), built by propagating the L-domain Gauss dynamics
% (lt_mee_rhs) with ode113 from the paper's initial GTO state at L0=pi
% (apogee) and sampling the dense ode113 solution at N+1 uniform-sigma nodes.
% Defect-free by construction: no per-node re-integration or downsampling
% interpolation (the campaign's no-resample lesson) -- L is already the
% collocation solver's independent variable, so dense-output sampling on L
% needs no Sundman-style remap.
%
% Two steering laws for beta (opts.betaMode), both evaluated fresh at every
% ode113 step and every sampled node (never frozen at the initial state):
%   'transverse' - beta = [0;1;0] in the local RTN frame -- a literally
%                  constant control.
%   'tangential' - beta = unit RTN projection of the local velocity
%                  direction (radial and transverse components of v-hat;
%                  normal component 0, since v lies in the instantaneous
%                  orbit plane). Mirrors seed_2body.m's vhat steering law;
%                  the RTN components of "point along velocity" rotate with
%                  the orbit, so this is recomputed at each L, not a single
%                  frozen vector replicated.
% Throttle opts.thr is a scalar held constant at every node in both modes.
%
% Two mutually-exclusive ways to set the integration span:
%   opts.nRev  - fixed span dL0 = 2*pi*nRev.
%   opts.stopP - integrate until P(L) reaches opts.stopP via an ode113
%                terminal event (P - stopP = 0, increasing direction);
%                dL0 = the achieved span.
%
% INPUTS:
%   par  - kepler_lt_params struct (.mu .Tmax .c .LU_km ...), with NO .L
%          field required (set internally per ode113 step / per node)
%   opts - struct: .thr [scalar, 0-1], .betaMode ['transverse'|'tangential'],
%          .N [scalar segments, nodes = N+1], and EITHER .nRev [scalar] OR
%          .stopP [scalar]. Optional .initElems [7x1] [P;ex;ey;hx;hy;m;t] --
%          overrides the initial-node MEE state at L0=pi; absent or empty
%          ([]) falls back to the paper's legacy literal state (Haberkorn
%          et al. GTO at apogee, i0~7deg parameterized as hx=0.0612).
%
% OUTPUTS:
%   sigma - uniform node parameter, 0->1 [(N+1)x1]
%   X0    - sampled MEE states [P;ex;ey;hx;hy;m;t] at each node [7x(N+1)]
%   U0    - sampled controls [beta(3);thr] at each node [4x(N+1)]
%   dL0   - total true-longitude span L(end)-L0 [scalar]
%   info  - struct: .nRev (=dL0/(2*pi)), .tEnd (=X0(7,end)), .mEnd (=X0(6,end))
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004 (initial GTO state at
%       apogee; Gauss dynamics per lt_mee_rhs.m).
%   [2] earth_elliptic_to_geo/seed_2body.m (Cartesian analog: dense-output
%       sampling, vhat steering law, "no-resample" lesson).

L0 = pi;                                    % apogee-start convention (all endpoints)
% Default = the paper's legacy literal, byte-for-byte (NOT recomputed from
% i0=7deg: tan(3.5deg)=0.061163 != the certified 0.0612). A caller wanting a
% different start orbit passes opts.initElems explicitly (run_gergaud builds
% it from P0/e0/i0 only when the user overrides the paper defaults).
if isfield(opts,'initElems') && ~isempty(opts.initElems)
    X_init = opts.initElems(:);
    assert(numel(X_init)==7, 'mee_seed: opts.initElems must be 7x1 [P;ex;ey;hx;hy;m;t]');
else
    X_init = [11625/par.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
end

thr      = opts.thr;
betaMode = opts.betaMode;
N        = opts.N;

odef = @(L, X) local_rhs(X, L, par, thr, betaMode);
oo   = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);

useStopP = isfield(opts, 'stopP') && ~isempty(opts.stopP);
if useStopP
    oo  = odeset(oo, 'Events', @(L, X) stopP_event(X, opts.stopP));
    sol = ode113(odef, [L0, L0 + 2*pi*1000], X_init, oo);
    assert(~isempty(sol.xe), 'mee_seed:noEvent', ...
        'stopP = %.6g not reached within 1000 revolutions', opts.stopP);
else
    sol = ode113(odef, [L0, L0 + 2*pi*opts.nRev], X_init, oo);
end

dL0     = sol.x(end) - L0;
sigma   = linspace(0, 1, N+1).';
Lq      = L0 + sigma*dL0;
Lq(end) = sol.x(end);                  % exact endpoint (avoid deval range error)
X0      = deval(sol, Lq.');

U0 = zeros(4, N+1);
for k = 1:N+1
    U0(:,k) = [local_beta(X0(:,k), Lq(k), par, betaMode); thr];
end

info.nRev = dL0/(2*pi);
info.tEnd = X0(7,end);
info.mEnd = X0(6,end);
end

% -----------------------------------------------------------------------
function dXdL = local_rhs(X, L, par, thr, betaMode)
% ODE right-hand side for ode113: builds U from the steering law, sets
% par.L, and delegates to lt_mee_rhs.
par.L = L;
U     = [local_beta(X, L, par, betaMode); thr];
dXdL  = lt_mee_rhs(X, U, par);
end

function beta = local_beta(X, L, par, betaMode)
% Steering law: unit RTN thrust direction at state X, longitude L.
switch betaMode
    case 'transverse'
        beta = [0; 1; 0];
    case 'tangential'
        [r, v] = elements_to_cart(X(1), X(2), X(3), X(4), X(5), L, par.mu);
        rhat = r / norm(r);
        hvec = cross(r, v);
        nhat = hvec / norm(hvec);
        that = cross(nhat, rhat);
        vhat = v / norm(v);
        beta = [dot(vhat, rhat); dot(vhat, that); 0];
        beta = beta / norm(beta);          % defensive renormalization
    otherwise
        error('mee_seed:betaMode', 'unknown opts.betaMode ''%s''', betaMode);
end
end

function [val, isterm, dirn] = stopP_event(X, stopP)
% Terminal ode113 event: fire when P crosses stopP from below.
val    = X(1) - stopP;
isterm = 1;
dirn   = 1;
end
