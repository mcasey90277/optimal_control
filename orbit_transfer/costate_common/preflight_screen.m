function [ok, why, minAltKm] = preflight_screen(o, muStar, lStar, floorKm, seedTf)
%% Purpose:
%
%   CHEAP SANITY PRE-CHECK on a direct solution BEFORE any integrator
%   touches it, made a single-home helper (migration #4). A solve can meet
%   the discrete defect test to 1e-9 and still be physically wild;
%   integrating such a trajectory crawls near the lunar singularity
%   WITHOUT BOUND (measured: one cell pinned a catalog run 16 min at 100%
%   CPU). Screens on the discrete nodes only -- altitude floor and a
%   plausible time of flight -- so it costs microseconds.
%
%% Inputs:
%
%  o                        struct                  Direct solution with
%                                                   .X [>=3 x N+1] and .tf
%
%  muStar                   double                  CR3BP mass ratio
%
%  lStar                    double                  Length unit (km)
%
%  floorKm                  double                  Campaign altitude floor
%                                                   (km); the screen trips
%                                                   below HALF of it
%
%  seedTf                   double or []            Previous rung's tf for
%                                                   the plausibility band
%                                                   [0.3, 3] x seedTf; []
%                                                   skips the band
%
%% Outputs:
%
%  ok                       logical                 Safe to integrate
%
%  why                      char                    '' | 'altitude' | 'tf'
%
%  minAltKm                 double                  Node-sampled minimum
%                                                   lunar altitude (km)
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

rMoonKm = 1737.4;
dMoon = vecnorm(o.X(1:3,:) - [1-muStar;0;0], 2, 1)*lStar - rMoonKm;
minAltKm = min(dMoon);
tfBad = ~isempty(seedTf) && (o.tf > 3*seedTf || o.tf < 0.3*seedTf);
ok = true;  why = '';
if minAltKm < 0.5*floorKm
    ok = false;  why = 'altitude';
elseif tfBad
    ok = false;  why = 'tf';
end
end
