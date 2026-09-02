function [yh, PHI, T, Y] = cr3bp_minfuel_prop(dt, y0, needSTM, Tmax, c, muStar, smooth)
%% Purpose:
%
%   Propagates the smoothed energy->fuel PMP state (and optionally its
%   14x14 STM) for a time dt -- the homotopy sibling of
%   cr3bp_minenergy_prop, in the [yh, PHI] = prop(dt, y0, needSTM) shape
%   ms_bvp expects. Variational equations ride along with A from
%   cr3bp_minfuel_pmp (exact AD). Integrator ode113 at tfMinProp's
%   tolerances (RelTol 1e-10, AbsTol 1e-12).
%
%  ASSUMPTIONS / NOTES:
%
% • Throws (via ode113) on integrator collapse -- ms_bvp converts the
%   throw into a rejected iterate, as required by its contract.
% • The 'huber' family's throttle law is DISCONTINUOUS at Q = 1 (see
%   cr3bp_minfuel_pmp); ode113 handles isolated crossings by step
%   rejection, at a cost the race is designed to measure.
%
%% Inputs:
%
%  dt                       double                  Propagation time, ND
%
%  y0                       [14 x 1]                Initial PMP state
%
%  needSTM                  logical                 Integrate the STM too
%
%  Tmax, c, muStar          double                  As cr3bp_minfuel_pmp
%
%  smooth                   struct                  {family, p} smoothing
%                                                   spec (cr3bp_minfuel_pmp)
%
%% Outputs:
%
%  yh                       [14 x 1]                State at dt
%
%  PHI                      [14 x 14] or []         STM d yh / d y0
%
%  T                        [nT x 1]                Integrator times
%
%  Y                        [nT x 14]               States at T
%
%% Revision History:
%  M. Casey                                                   (c) 09/02/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: one arc at eps = 0.3:
     mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
     y_ = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
           15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
     [yh_, PHI_] = cr3bp_minfuel_prop(0.4, y_, true, T_, c_, mu_, ...
                                      struct('family','eps','p',0.3));
     fprintf('demo: |r(0.4)| = %.6f, m = %.6f, |PHI| = %.3e\n', ...
             sqrt(sum(yh_(1:3).^2)), yh_(7), max(abs(PHI_(:))));
     return
end

y0 = y0(:);
if dt == 0
    yh = y0;  PHI = eye(14);  T = 0;  Y = y0';
    if ~needSTM, PHI = []; end
    return
end
odeOpts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
if needSTM
    z0 = [y0; reshape(eye(14), [], 1)];
    [T, Z] = ode113(@(t, z) rhsSTM(z, Tmax, c, muStar, smooth), [0 dt], z0, odeOpts);
    yh  = Z(end, 1:14)';
    PHI = reshape(Z(end, 15:210), 14, 14);
    Y   = Z(:, 1:14);
else
    [T, Y] = ode113(@(t, y) cr3bp_minfuel_pmp(y, Tmax, c, muStar, smooth), ...
                    [0 dt], y0, odeOpts);
    yh  = Y(end, :)';
    PHI = [];
end
end

% ------------------------------------------------------------------------
function dz = rhsSTM(z, Tmax, c, muStar, smooth)
% RHSSTM  PMP state + STM variational RHS.  INPUTS: z [210x1]; field
% params + smooth.  OUTPUTS: dz [210x1].
[F, A] = cr3bp_minfuel_pmp(z(1:14), Tmax, c, muStar, smooth);
PHI = reshape(z(15:210), 14, 14);
dz = [F; reshape(A*PHI, [], 1)];
end
