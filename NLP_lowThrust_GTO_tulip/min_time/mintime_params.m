function [rv0, rvf_tul, P] = mintime_params()
% MINTIME_PARAMS  CR3BP low-thrust constants + GTO departure + tulip min-time
% target, self-contained for the min_time module (mirrors gto_tulip_endpoints /
% ztl_endpoints so this directory does not depend on the other modules).
%
% OUTPUTS:
%   rv0     - GTO departure state (350 x 35786 km, argp -25 deg), ND rot [1x6]
%   rvf_tul - tulip max-ydot rendezvous point, ND rot [1x6]
%   P       - constants: .muStar .lStar .tStar .m0kg .g0 .c .Tmax25
%             (Tmax25 = ND thrust accel at m=1 for the nominal 25 mN / 15 kg)
%
% REFERENCES:
%   pumpkyn.cr3bp (Koblick) - orb2eci, fromPCI, getTulip, prop.

P.muStar = 0.012150585609624;
P.lStar  = 389703.264829278;
P.tStar  = 382981.289129055;
P.m0kg   = 15;
P.g0     = 9.80665*P.tStar^2/(1000*P.lStar);
P.c      = (2100/P.tStar)*P.g0;
P.Tmax25 = (0.025/P.m0kg)*P.tStar^2/(P.lStar*1000);

muEarth = 6.67384e-20*(1 - P.muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;   ecc = (35786-350)/(2*sma);
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0, v0], P.muStar, P.tStar, P.lStar, 1);

[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, yTul]    = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, P.muStar);
[~, idxF]    = max(yTul(:,5));
rvf_tul = yTul(idxF, 1:6);
end
