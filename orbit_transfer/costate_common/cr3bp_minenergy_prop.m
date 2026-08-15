function [yh, PHI, T, Y] = cr3bp_minenergy_prop(dt, y0, needSTM, Tmax, c, muStar)
%% Purpose:
%
%   Propagates the min-energy PMP state (and optionally its 14x14 STM) for
%   a time dt -- the fixed-tf sibling of pumpkyn.cr3bp.tfMinProp, in the
%   [yh, PHI] = prop(dt, y0, needSTM) shape ms_bvp expects. Variational
%   equations dPHI/dt = A(y) PHI ride along as 196 extra states, with A the
%   exact AD Jacobian from cr3bp_minenergy_pmp. Integrator ode113 at
%   tfMinProp's tolerances (RelTol 1e-10, AbsTol 1e-12).
%
%  ASSUMPTIONS / NOTES:
%
% • Throws (via ode113) on integrator collapse -- ms_bvp converts the
%   throw into a rejected iterate, as required by its contract.
% • dt = 0 returns y0 and the identity.
%
%% Inputs:
%
%  dt                       double                  Propagation time, ND
%
%  y0                       [14 x 1]                Initial PMP state
%
%  needSTM                  logical                 Integrate the STM too
%
%  Tmax, c, muStar          double                  As cr3bp_minenergy_pmp
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
%  M. Casey                                                   (c) 08/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: propagate the cr3bp_minenergy_pmp demo state 0.4 ND with STM.
     mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
     y_ = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
           15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
     tic; [yh_, PHI_] = cr3bp_minenergy_prop(0.4, y_, true, T_, c_, mu_); tw = toc;
     fprintf('demo: |r(0.4)| = %.6f, m = %.6f, |PHI| = %.3e, %.2f s\n', ...
             sqrt(sum(yh_(1:3).^2)), yh_(7), max(abs(PHI_(:))), tw);
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
    [T, Z] = ode113(@(t, z) rhsSTM(z, Tmax, c, muStar), [0 dt], z0, odeOpts);
    yh  = Z(end, 1:14)';
    PHI = reshape(Z(end, 15:210), 14, 14);
    Y   = Z(:, 1:14);
else
    [T, Y] = ode113(@(t, y) cr3bp_minenergy_pmp(y, Tmax, c, muStar), ...
                    [0 dt], y0, odeOpts);
    yh  = Y(end, :)';
    PHI = [];
end
end

% ------------------------------------------------------------------------
function dz = rhsSTM(z, Tmax, c, muStar)
% RHSSTM  State + variational right-hand side.
% INPUTS: z [210x1]; Tmax; c; muStar.  OUTPUTS: dz [210x1].
[F, A] = cr3bp_minenergy_pmp(z(1:14), Tmax, c, muStar);
PHI = reshape(z(15:210), 14, 14);
dz = [F; reshape(A*PHI, [], 1)];
end
