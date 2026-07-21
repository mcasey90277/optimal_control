function [rv0, rvf, P] = ztl_endpoints()
% ZTL_ENDPOINTS  Problem endpoints + constants for the GTO->tulip transfer.
%
% Same construction as thrust_continuation_minfuel_indirect.m and
% run_gto_tulip_indirect.m: GTO departure state (350 x 35786 km, -25 deg
% RAAN, perigee) and the tulip max-ydot rendezvous point, in the ND
% Earth-Moon rotating frame. Factored out so every ZTL stage uses
% identical endpoints.
%
% INPUTS:
%   (none)
% OUTPUTS:
%   rv0 - initial position/velocity, ND rotating frame [1x6]
%   rvf - target position/velocity (tulip max-ydot point), ND [1x6]
%   P   - constants struct:
%         .muStar .lStar .tStar  - CR3BP normalization
%         .m0kg                  - initial mass [kg]
%         .g0 .c                 - ND gravity, ND exhaust velocity
%         .Tmax25                - ND thrust accel at m=1 for 25 mN
%
% REFERENCES:
%   [1] thrust_continuation_minfuel_indirect.m (source of this block).

P.muStar = 0.012150585609624;
P.lStar  = 389703.264829278;
P.tStar  = 382981.289129055;
P.m0kg   = 15;
P.g0     = 9.80665*P.tStar^2/(1000*P.lStar);
P.c      = (2100/P.tStar)*P.g0;
P.Tmax25 = (0.025/P.m0kg)*P.tStar^2/(P.lStar*1000);

muEarth = 6.67384e-20*(1 - P.muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;
ecc = (35786-350)/(2*sma);
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0, v0], P.muStar, P.tStar, P.lStar, 1);

[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, rvTgt]   = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, P.muStar);
[~, idxF]    = max(rvTgt(:,5));
rvf = rvTgt(idxF, :);
end
