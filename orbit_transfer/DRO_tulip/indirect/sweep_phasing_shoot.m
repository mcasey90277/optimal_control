function [r, J, B] = sweep_phasing_shoot(z, rv0, rvf, Tmax, c, muStar)
% SWEEP_PHASING_SHOOT  Shooting kernel for the phasing sweep: residual,
% analytic STM Jacobian, and endpoint sensitivity.
%
% Residual (identical to pumpkyn tfMin's shootingResidual):
%   r = [ y_f(1:6) - rvf ;  lambda_m(t_f) ;  H(t_f) ]
% J = dr/d[lambda0; tf]      -- the Newton matrix (8x8)
% B = dr/d(rv0)              -- endpoint sensitivity (8x6), used by the
%                               continuation TANGENT PREDICTOR: without it a
%                               departure-phase step presents the corrector
%                               with an O(1) residual (initial-state errors
%                               amplify ~1e3x through the flow) and Newton
%                               trials wander into near-singular lunar passes
%                               that stall the integrator.
% All dynamics via pumpkyn's public calls (tfMinProp takes the 14+196
% augmented state; tfMinEoM returns dH/dy). Pumpkyn is not modified.
%
% INPUTS:  z [8x1] = [lambda0(7); tf]; rv0, rvf [1x6]; Tmax; c; muStar
% OUTPUTS: r [8x1]; J [8x8]; B [8x6]
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMin shootingResidual (the assembly mirrored here).
z = z(:);  tf = z(8);
if ~isfinite(tf) || tf <= 0.5 || tf > 20
    % A trial this far afield is never useful and can stall the integrator
    % near the singularity -- fail fast instead.
    r = 1e3*ones(8,1);  J = eye(8)*1e3;  B = zeros(8,6);  return
end
% THE RESIDUAL ALWAYS COMES FROM THE PLAIN 14-STATE PROPAGATION. Propagating
% the 210-dim augmented vector changes the adaptive integrator's error control
% (196 STM components of magnitude up to ~1e3 enter the norm), and the two
% trajectories disagree at the ~1e-3 level -- measured: the corrector
% "converged" its augmented residual while the plain verification residual
% stalled at 1.7e-3 on every edge. Newton tolerates an approximate Jacobian;
% it cannot tolerate an approximate residual. So: r from the plain
% propagation, J (and B) from a second, augmented one.
y0 = [rv0(:); 1; z(1:7)];
[~, Yp] = pumpkyn.cr3bp.tfMinProp(tf, y0, Tmax, c, muStar);
yf = Yp(end,1:14).';
[Ff, Ht, dHdy] = pumpkyn.cr3bp.tfMinEoM(tf, yf, Tmax, c, muStar);
r = [yf(1:6) - rvf(:); yf(14); Ht];
if nargout > 1
    [~, Y] = pumpkyn.cr3bp.tfMinProp(tf, [y0; reshape(eye(14),[],1)], Tmax, c, muStar);
    PHI = reshape(Y(end,15:210), 14, 14);
    J = [PHI(1:6, 8:14),    Ff(1:6);
         PHI(14,  8:14),    Ff(14);
         dHdy*PHI(:, 8:14), dHdy*Ff];
    if nargout > 2
        B = [PHI(1:6, 1:6);
             PHI(14,  1:6);
             dHdy*PHI(:, 1:6)];
    end
end
end
