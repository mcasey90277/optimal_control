function tJ = arc_boundaries_tau(t, rXYZ, M, muStar, pSund)
% ARC_BOUNDARIES_TAU  Arc joint times uniform in the Sundman variable tau.
%
% dt/dtau = kappa = r1^pSund (r1 = Earth distance), so uniform-tau joints
% cluster in physical time near perigee, where trajectory sensitivity
% accrues. Computed along a given reference trajectory.
%
% INPUTS:
%   t      - trajectory times, strictly increasing [1xK]
%   rXYZ   - rotating-frame positions along the trajectory [3xK]
%   M      - number of arcs [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%   pSund  - (optional) Sundman exponent [scalar, default 1.5]
%
% OUTPUTS:
%   tJ - joint times [1x(M+1)], tJ(1)=t(1), tJ(end)=t(end)

if nargin < 5 || isempty(pSund), pSund = 1.5; end
r1    = sqrt(sum((rXYZ - [-muStar; 0; 0]).^2, 1));
kappa = r1.^pSund;
tau   = cumtrapz(t(:), 1./kappa(:));          % strictly increasing (kappa > 0)
tauJ  = linspace(0, tau(end), M+1);
tJ    = interp1(tau, t(:), tauJ);
tJ(1)   = t(1);
tJ(end) = t(end);
end
