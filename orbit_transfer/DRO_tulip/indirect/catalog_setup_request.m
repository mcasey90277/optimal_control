function so = catalog_setup_request(cat_, sD1)
%% Purpose:
%
%   The arclength_arrival('setup', so) request for this catalog's problem:
%   engine, orbits and the spine's departure phase, CLOSURES ONLY.
%
%  ASSUMPTIONS / NOTES:
%
% • The filler needs the endpoint closures and the propulsion constants,
%   never an anchor; asking for one made it re-polish the shipped 70 mN
%   anchor, which fails for any other engine, orbit pair or departure phase
%   (FINDINGS 78).
% • Reads sheet 1 of the catalog for the orbit keys.
%
%% Inputs:
%
%  cat_                     struct                  catalog
%  sD1                      double                  spine departure phase
%
%% Outputs:
%
%  so                       struct                  request for arclength_arrival('setup', so)
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

sh = cat_.sheets(1);
so = struct('thrustN', cat_.rungs_N(1), 'ispS', cat_.thruster.isp_s, 'm0kg', cat_.thruster.m0_kg, ...
            'tauDRO', sh.tauDRO, 'NpTulip', sh.Np, 'pmTulip', sh.pm, 'sD', sD1, 'physicsOnly', true);
end
