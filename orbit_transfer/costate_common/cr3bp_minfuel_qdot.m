function Qd = cr3bp_minfuel_qdot(y, Tmax, c)
%% Purpose:
%
%   Time derivative of the switch quantity Q = Tmax(|lam_v|/m + lam_m/c)
%   along the CR3BP min-fuel PMP flow, in CLOSED FORM:
%
%       Qdot = - Tmax * (lam_v' lam_r) / (m |lam_v|)
%
%   The throttle drops out exactly: the mass-rate term of d(|lam_v|/m)/dt
%   and lam_m-dot cancel (both are s*Tmax*|lam_v|/(c m^2)), and the Coriolis
%   contribution to lam_v-dot is skew so lam_v' Omega lam_v = 0. Hence Qdot
%   is CONTINUOUS across a throttle switch, n'(F+ - F-) = 0 for the switch
%   normal n = grad Q (measured 9e-16), and this is the exact transversality
%   (saltation denominator) at any crossing -- use it instead of a finite
%   difference of Q along the trajectory. (GPT-6 Astra review #2,
%   2026-09-06; the field itself is cr3bp_minfuel_pmp.)
%
%% Inputs:
%
%  y                        [14 x N]                [r; v; m; lam_r; lam_v;
%                                                   lam_m] (columns)
%
%  Tmax                     double                  ND thrust acceleration
%                                                   at unit mass fraction
%
%  c                        double                  ND exhaust velocity
%                                                   (enters Q, not Qdot;
%                                                   kept for signature
%                                                   symmetry with Q)
%
%% Outputs:
%
%  Qd                       [1 x N]                 dQ/dt
%
%% Revision History:
%  M. Casey                                                   (c) 09/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo:
     y_ = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
           15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
     fprintf('demo: Qdot = %.6f\n', cr3bp_minfuel_qdot(y_, 0.1756418, 8.673746));
     return
end
%#ok<*INUSD>
rho = sqrt(sum(y(11:13, :).^2, 1) + 1e-300);
Qd  = -Tmax * sum(y(11:13, :) .* y(8:10, :), 1) ./ (y(7, :) .* rho);
end
