function [R, J] = shoot_residual_tf(z, rv0, rvf, Tmax, c, muStar)
% SHOOT_RESIDUAL_TF  Shooting residual (+ Jacobian) for the min-time TPBVP.
%
% Integrates the augmented PMP dynamics from the fixed initial state with
% candidate initial costates z(1:7) over the candidate transfer time z(8),
% and returns the 8 terminal conditions of the free-final-time minimum-time
% problem: rendezvous in position and velocity, free-final-mass
% transversality lambda_m(tf) = 0, and free-final-time transversality
% H(tf) = 0.
%
% The Jacobian is computed by COMPLEX-STEP differentiation through the
% integrator for the seven costate columns (machine-accurate, no hand-coded
% 14x14 variational system needed -- requires the EoM to be complex-safe:
% sqrt(sum(x.^2)) not norm, .' not '), and analytically for the tf column:
% dR/dtf = [f_rv(yf); f_lambda_m(yf); 0], since the terminal conditions
% shift with the flow itself and H is conserved along the (autonomous)
% canonical flow.
%
% This problem is EXTREMELY sensitive: ~40 perigee passes amplify initial
% perturbations by ~1e6, so ordinary finite differencing across the
% adaptive integrator produces unusable derivatives (fsolve stalls). The
% complex step sidesteps subtractive cancellation entirely.
%
% INPUTS:
%   z      - decision vector [8x1]: [lambda_r0(3); lambda_v0(3);
%            lambda_m0; tf] (all ND)
%   rv0    - initial position/velocity (ND, rotating frame) [1x6 or 6x1]
%   rvf    - target position/velocity at tf (ND) [1x6 or 6x1]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   R      - residual [8x1]: [r(tf)-rf; v(tf)-vf; lambda_m(tf); H(tf)]
%   J      - (optional) dR/dz [8x8]
%
% REFERENCES:
%   [1] Martins, Sturdza, Alonso, "The Complex-Step Derivative
%       Approximation," ACM TOMS 29(3), 2003.
%   [2] pumpkyn.cr3bp.tfMin (reference: same residual with analytic
%       STM-based Jacobian).

z       = z(:);
lambda0 = z(1:7);
tf      = z(8);

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

R = terminalResidual(lambda0, tf, rv0, rvf, Tmax, c, muStar, opts);

if nargout > 1
    J    = zeros(8, 8);
    hCS  = 1e-20;                       % complex-step size (no cancellation)
    for kZ = 1:7
        scale         = max(1, abs(lambda0(kZ)));
        lambdaPert    = lambda0;
        lambdaPert(kZ) = lambdaPert(kZ) + 1i*hCS*scale;
        Rp            = terminalResidual(lambdaPert, tf, rv0, rvf, ...
                                         Tmax, c, muStar, opts);
        J(:, kZ)      = imag(Rp)./(hCS*scale);
    end
    % tf column: terminal conditions ride the flow; H is conserved
    % (autonomous canonical system), so its tf-derivative is exactly 0.
    yf      = integrateArc(lambda0, tf, rv0, Tmax, c, muStar, opts);
    yfDot   = lt_pmp_eom(tf, yf, Tmax, c, muStar);
    J(:, 8) = [yfDot(1:6); yfDot(14); 0];
end
end

% -------------------------------------------------------------------------
function R = terminalResidual(lambda0, tf, rv0, rvf, Tmax, c, muStar, opts)
yf       = integrateArc(lambda0, tf, rv0, Tmax, c, muStar, opts);
[~, Htf] = lt_pmp_eom(tf, yf, Tmax, c, muStar);
R = [yf(1:6) - rvf(:);   % rendezvous: position + velocity match
     yf(14);             % transversality: final mass free -> lambda_m(tf) = 0
     Htf];               % transversality: final time free -> H(tf) = 0
end

% -------------------------------------------------------------------------
function yf = integrateArc(lambda0, tf, rv0, Tmax, c, muStar, opts)
y0     = [rv0(:); 1; lambda0(:)];
[~, Y] = ode113(@lt_pmp_eom, [0 tf], y0, opts, Tmax, c, muStar);
yf     = Y(end, :).';
end
