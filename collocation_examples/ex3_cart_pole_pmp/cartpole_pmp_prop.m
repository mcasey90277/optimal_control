function [yEnd, PHI] = cartpole_pmp_prop(dt, y0, needSTM, p)
%% Purpose:
%
%   Propagate the cart-pole PMP flow, with the 8x8 state-transition matrix
%   on request -- the propagator contract costate_common/ms_bvp expects as
%   prob.prop.
%
%   The variational equations are integrated ALONGSIDE the state (72
%   equations in all), with the Jacobian taken by COMPLEX STEP through
%   cartpole_pmp_rhs. That is safe here only because the state Jacobian
%   inside the field is generated symbolically: a complex step inside a
%   complex step shares the imaginary channel and silently corrupts the
%   inner derivative.
%
%  ASSUMPTIONS / NOTES:
%
% • Tolerances 1e-12 / 1e-14: the STM feeds a Newton step, and an integrator
%   sloppier than the step it informs turns quadratic convergence into a
%   crawl.
% • THROWS on a non-finite result, per the engine's contract -- the solver
%   converts a throw into a rejected iterate and backtracks. An ode113
%   failure on this problem does NOT raise a MATLAB error and does NOT
%   produce non-finite values: it emits a warning and returns whatever it
%   integrated up to the point the step size collapsed below the minimum,
%   silently short of the requested dt. A finite-value check alone misses
%   that case, so the real gate is whether the returned time grid reached
%   dt.
% • The try/catch around ode113 is a narrow backstop, not a second
%   detector: it relabels only identifiers ode113 itself raises
%   ('MATLAB:ode*', e.g. bad sizes or options) as
%   cartpole_pmp_prop:collapse. A genuine programming error inside
%   cartpole_pmp_rhs/rhs_with_stm (undefined variable, missing p field,
%   dimension mismatch) is rethrown UNCHANGED, with its own identifier --
%   ms_bvp's residual() does a bare catch on prob.prop, so relabelling
%   every error as "collapse" would make a real bug indistinguishable from
%   a rejected iterate and surface only as mysterious non-convergence.
%
%% Inputs:
%
%  dt                       double                  Propagation time (s);
%                                                   0 returns y0 and I
%
%  y0                       [8 x 1]                 [x; lam] at the start
%
%  needSTM                  logical                 Return PHI as well
%
%  p                        struct                  .m1 .m2 .L .g
%
%% Outputs:
%
%  yEnd                     [8 x 1]                 The propagated state
%
%  PHI                      [8 x 8] or []           d(yEnd)/d(y0)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: half a second of the flow, with its STM conditioning:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [yE, PH] = cartpole_pmp_prop(0.5, [0;0;0;0; 0.8; -1.5; 0.4; 0.2], true, pd);
     fprintf('y(0.5) = [%s]\n  cond(PHI) = %.3e\n', sprintf('%.4f ', yE), cond(PH));
     return
end

if dt == 0
    yEnd = y0(:);
    PHI = [];  if needSTM, PHI = eye(8); end
    return
end

oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
try
    if needSTM
        z0 = [y0(:); reshape(eye(8), 64, 1)];
        [T, Z] = ode113(@(t, z) rhs_with_stm(z, p), [0 dt], z0, oo);
        zEnd = Z(end,:).';
        yEnd = zEnd(1:8);
        PHI  = reshape(zEnd(9:72), 8, 8);
    else
        [T, Y] = ode113(@(t, y) cartpole_pmp_rhs(y, p), [0 dt], y0(:), oo);
        yEnd = Y(end,:).';
        PHI  = [];
    end
catch err
    if startsWith(err.identifier, 'MATLAB:ode')
        error('cartpole_pmp_prop:collapse', ...
              'the PMP flow failed to integrate over dt = %g: %s', dt, err.message);
    end
    rethrow(err);
end

reachedEnd = abs(T(end) - dt) <= 1e-9*max(1, abs(dt));
if ~reachedEnd || ~all(isfinite(yEnd)) || (needSTM && ~all(isfinite(PHI(:))))
    error('cartpole_pmp_prop:collapse', ...
          'the PMP flow left the finite range over dt = %g (integrator reached t = %g)', ...
          dt, T(end));
end
end

% ------------------------------------------------------------------------
function dz = rhs_with_stm(z, p)
%% Purpose:
%
%   The state and its variational equations: zdot = [f(y); A(y) PHI], with
%   A taken by complex step through the PMP field.
%
y = z(1:8);
PHI = reshape(z(9:72), 8, 8);
dy = cartpole_pmp_rhs(y, p);
A = zeros(8);
h = 1e-20;
for col = 1:8
    yp = y;  yp(col) = yp(col) + 1i*h;
    A(:,col) = imag(cartpole_pmp_rhs(yp, p))/h;
end
dz = [dy; reshape(A*PHI, 64, 1)];
end
