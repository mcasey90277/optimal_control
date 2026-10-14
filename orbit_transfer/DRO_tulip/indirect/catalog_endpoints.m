function [B, orbits] = catalog_endpoints(cat_, sheetIdx)
%% Purpose:
%
%   The CERTIFIED endpoint rule rebuilt from a catalog's OWN keys: the
%   departure and arrival orbits (family + parameters stored on the sheet),
%   the phase-fraction closures stateD / stateA / dstateA, and the ND
%   propulsion constants. On the library of record these closures equal
%   arclength_arrival('setup', catalog_setup_request(...)) bitwise, so the
%   shipped tools fly from exactly the states the library was certified on.
%
%  ASSUMPTIONS / NOTES:
%
% • Rule (ladder_endpoints + arclength_arrival 'setup'): get_family_orbit
%   with the sheet's params plus constants.muStar, then costate_common/
%   phase_state with onMissing 'notaknot'.
% • phase_state returns the full row of the orbit table; the closures here
%   return the 6 x 1 state.
% • Tnd = (rungs_N(1)/m0_kg)*tStar^2/(lStar_km*1000) is the catalog's
%   derive.thrust_nd formula at the first rung; cnd = thruster.c_nd.
% • Phases are FRACTIONS of a period and wrap.
% • The wrap seam (s -> 1 vs s = 0, about 1.2e-9 ND) equals the orbit table's own
%   closure error, shared bitwise with the reference rule.
%
%% Inputs:
%
%  cat_                     struct                  catalog (.sheets .constants
%                                                   .thruster .rungs_N)
%  sheetIdx                 double [1]              sheet to read the orbit
%                                                   keys from [1]
%
%% Outputs:
%
%  B                        struct                  .stateD(sD) .stateA(sA)
%                                                   [6x1] .dstateA(sA) .Tnd
%                                                   .cnd .mu .problem (the
%                                                   fields plot_transfer_3d
%                                                   and the study tools read)
%
%  orbits                   struct                  .tD .rvD (departure) .tT
%                                                   .rvT (arrival), one period
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2 || isempty(sheetIdx), sheetIdx = 1; end

sh = cat_.sheets(sheetIdx);
k  = cat_.constants;
muStar = k.muStar;

pD = sh.dep_params;  pD.muStar = muStar;
pA = sh.arr_params;  pA.muStar = muStar;
[tD, rvD] = get_family_orbit(sh.dep_family, pD);
[tT, rvT] = get_family_orbit(sh.arr_family, pA);

ppo = struct('onMissing', 'notaknot');
[xD, ~]        = phase_state(tD, rvD, ppo);
[xA, ~, dxA]   = phase_state(tT, rvT, ppo);
col6 = @(v) reshape(v(1:6), 6, 1);

lStar = k.lStar_km;  tStar = k.tStar_s;
thrustN = cat_.rungs_N(1);
Tnd = (thrustN/cat_.thruster.m0_kg)*tStar^2/(lStar*1000);
cnd = cat_.thruster.c_nd;

B = struct();
B.stateD  = @(s) col6(xD(s));
B.stateA  = @(s) col6(xA(s));
B.dstateA = @(s) col6(dxA(s));
B.Tnd = Tnd;  B.cnd = cnd;  B.mu = muStar;
B.problem = struct('thrustN', thrustN, 'ispS', cat_.thruster.isp_s, 'm0kg', cat_.thruster.m0_kg, ...
    'tauDRO', sh.tauDRO, 'NpTulip', sh.Np, 'pmTulip', sh.pm, 'muStar', muStar, ...
    'lStar', lStar, 'tStar', tStar, 'Tnd', Tnd, 'cnd', cnd, ...
    'periodTulip', sh.period_tulip_nd, 'sD', NaN);
orbits = struct('tD', tD, 'rvD', rvD, 'tT', tT, 'rvT', rvT);
end
