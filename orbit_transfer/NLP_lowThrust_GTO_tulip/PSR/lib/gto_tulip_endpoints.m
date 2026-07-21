function [rv0, rvf, tulipTrace] = gto_tulip_endpoints(p)
% GTO_TULIP_ENDPOINTS  Boundary states for the GTO -> south-pole tulip transfer.
%
% Builds the standard problem endpoints in the rotating Earth-Moon CR3BP frame:
% departure from a GTO (350 x 35786 km, argument of perigee -25 deg) and
% arrival at the maximum-ydot point of a (5/6)-resonant south-pole tulip orbit.
% Requires the pumpkyn toolbox on the path (call setup_paths first).
%
% INPUTS:
%   p - parameter struct from CR3BP_LT_PARAMS (uses .muStar,.lStar,.tStar)
%
% OUTPUTS:
%   rv0        - departure state [1x6] (ND, rotating frame) [r v]
%   rvf        - arrival tulip state [1x6] (ND, rotating frame) [r v]
%   tulipTrace - full tulip orbit trace [Kx6] for plotting (optional)
%
% REFERENCES:
%   [1] pumpkyn.cr3bp (Koblick) - orb2eci, fromPCI, getTulip, prop.

muEarth = 6.67384e-20*(1-p.muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;   ecc = (35786-350)/(2*sma);
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0, v0], p.muStar, p.tStar, p.lStar, 1);

[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, yTul]    = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, p.muStar);
[~, idxF]    = max(yTul(:,5));           % maximum-ydot point
rvf = yTul(idxF, 1:6);
tulipTrace = yTul(:, 1:6);
end
