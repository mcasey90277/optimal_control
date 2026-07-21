function [R, J] = shoot_residual_energy(lambda0, tf, rv0, m0, rvf, Tmax, c, muStar)
% SHOOT_RESIDUAL_ENERGY  Shooting residual (+ Jacobian) for min-energy.
%
% Fixed-final-time minimum-energy TPBVP: integrate the smooth min-energy
% PMP dynamics (LT_PMP_EOM_ENERGY) from the fixed initial state with
% candidate initial costates lambda0, and return the 7 terminal
% conditions: rendezvous in position and velocity plus free-final-mass
% transversality lambda_m(tf) = 0. tf is fixed (no H(tf) row). Jacobian by
% complex step, as in SHOOT_RESIDUAL_TF / SHOOT_RESIDUAL_MINFUEL. No
% smoothing parameter: the min-energy control is already continuous.
%
% INPUTS:
%   lambda0 - initial costates [7x1] (ND)
%   tf      - FIXED transfer time (ND) [scalar]
%   rv0     - initial position/velocity (ND) [1x6 or 6x1]
%   m0      - initial mass FRACTION (1 for a fresh start) [scalar]
%   rvf     - target position/velocity at tf (ND) [1x6 or 6x1]
%   Tmax    - max thrust acceleration at m = 1 (ND) [scalar]
%   c       - exhaust velocity (ND) [scalar]
%   muStar  - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   R - residual [7x1]: [r(tf)-rf; v(tf)-vf; lambda_m(tf)]
%   J - (optional) dR/dlambda0 [7x7]
%
% REFERENCES:
%   [1] Martins, Sturdza, Alonso, ACM TOMS 29(3), 2003 (complex step).
%   [2] Caillau, Gergaud, Noailles, JOTA 2003 (min-energy transfer).

lambda0 = lambda0(:);
opts    = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

R = terminalResidual(lambda0, tf, rv0, m0, rvf, Tmax, c, muStar, opts);

if nargout > 1
    J   = zeros(7, 7);
    hCS = 1e-20;
    for kZ = 1:7
        scale          = max(1, abs(lambda0(kZ)));
        lambdaPert     = lambda0;
        lambdaPert(kZ) = lambdaPert(kZ) + 1i*hCS*scale;
        Rp             = terminalResidual(lambdaPert, tf, rv0, m0, rvf, ...
                                          Tmax, c, muStar, opts);
        J(:, kZ)       = imag(Rp)./(hCS*scale);
    end
end
end

% -------------------------------------------------------------------------
function R = terminalResidual(lambda0, tf, rv0, m0, rvf, Tmax, c, muStar, opts)
y0     = [rv0(:); m0; lambda0(:)];
[~, Y] = ode113(@lt_pmp_eom_energy, [0 tf], y0, opts, Tmax, c, muStar);
yf     = Y(end, :).';
R = [yf(1:6) - rvf(:);   % rendezvous
     yf(14)];            % m(tf) free -> lambda_m(tf) = 0
end
