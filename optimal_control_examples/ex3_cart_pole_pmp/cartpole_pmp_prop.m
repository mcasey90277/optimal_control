function [yEnd, PHI] = cartpole_pmp_prop(dt, y0, needSTM, p)
%% Purpose:
%
%   Propagate the cart-pole PMP flow, with the 8x8 state-transition matrix
%   on request -- the propagator contract oc.ms_bvp expects as
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
% • THROWS cartpole_pmp_prop:collapse when the integration is truncated (the
%   returned grid does not end exactly at dt -- ode113 warns and stops early
%   on a blow-up rather than producing non-finite values) or the result is
%   non-finite. Per the engine's contract a throw is a rejected iterate.
%   Every OTHER error propagates unchanged: nothing is relabelled.
% • REAL base state only when needSTM: the STM's complex step through the
%   field assumes y is real, so this function is not itself complex-steppable
%   on that path. Inputs are validated (cartpole_pmp_prop:input).
% • WHY the reached-dt gate and not a finiteness check alone: an ode113
%   failure on this problem does not raise an error and does not produce
%   non-finite values. It warns, and returns whatever it integrated before
%   the step size collapsed -- finite, plausible, and silently short of dt.
% • No try/catch around ode113: it has nothing to catch on a collapse (see
%   above), and relabelling whatever else is thrown would turn a coding bug
%   in the field into a "rejected iterate" -- oc.ms_bvp's residual() catches
%   everything from prob.prop -- that surfaces only as non-convergence.
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
     here = fileparts(mfilename('fullpath'));
     addpath(here, fullfile(fileparts(here), 'cartpole_common'));
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [yE, PH] = cartpole_pmp_prop(0.5, [0;0;0;0; 0.8; -1.5; 0.4; 0.2], true, pd);
     fprintf('y(0.5) = [%s]\n  cond(PHI) = %.3e\n', sprintf('%.4f ', yE), cond(PH));
     return
end

% A malformed request is refused by name BEFORE any shortcut: the dt == 0
% branch would otherwise hand a non-finite or mis-sized y0 straight back.
if ~(isnumeric(y0) && isreal(y0) && numel(y0) == 8 && all(isfinite(y0(:))))
    error('cartpole_pmp_prop:input', 'y0 must be 8 real finite values');
end
if ~(isnumeric(dt) && isreal(dt) && isscalar(dt) && isfinite(dt))
    error('cartpole_pmp_prop:input', 'dt must be a real finite scalar');
end

if dt == 0
    yEnd = y0(:);
    PHI = [];  if needSTM, PHI = eye(8); end
    return
end

% No try/catch here, deliberately. ode113 does not THROW when this flow blows
% up -- it warns and returns a truncated grid -- so the detectors below are
% the whole collapse contract, and anything that does throw (a coding error
% in the field, a bad option) propagates with its own identifier instead of
% being relabelled as a collapse the engine would quietly backtrack over.
oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
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

% EXACT, not a tolerance band: on success ode113 returns the requested
% endpoint bitwise (measured), and a band in seconds would accept a
% truncation on a short segment, where a tiny time deficit near a collapse
% need not mean a tiny state error.
reachedEnd = (T(end) == dt);
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
