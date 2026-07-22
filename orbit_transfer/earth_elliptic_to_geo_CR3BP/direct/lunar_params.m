function pert = lunar_params(par, phi0, gain)
% LUNAR_PARAMS  Lunar third-body constants in the 2-body campaign's units.
%
% Circular Moon in the reference plane (spec D2): geocentric distance D_EM,
% sidereal rate n_M = sqrt((mu_E+mu_M)/D_EM^3), phase phi0 at t=0. Expressed
% in kepler_lt_params canonical units (mu_E = 1, LU = 42165 km GEO radius).
% The continuation scale `gain` multiplies mu_M (spec D5: gain IS the
% mu-continuation knob; gain=1 is the physical Moon).
%
% INPUTS:
%   par  - kepler_lt_params struct (.mu .muKm3s2 .LU_km .TU_s) [struct]
%   phi0 - lunar phase at t=0 [rad, scalar] (spec D6; baseline 0)
%   gain - mu_M continuation scale in [0,1] [scalar]
%
% OUTPUTS:
%   pert - struct .muM (canonical, PHYSICAL value -- gain applied in the
%          RHS, not here) .DM .nM (canonical) .phi0 .gain
%
% NOTE (amendment G): the governing Earth GM is the SOLVER's
% par.muKm3s2 = 398600.47 km^3/s^2 (kepler_lt_params), not the spec's
% 398600.4418 -- a 7e-8 relative difference. This function reads par.muKm3s2
% (never a hardcoded spec constant), so 398600.47 is what actually
% nondimensionalizes .muM/.nM above; test_lunar_params.m:7 bakes this value
% into its mass-ratio assertion.
%
% NOTE (phi0 epoch convention): phi0 is "Moon at angle phi0 in the reference
% plane when the seed's time-state t = X(7) is exactly 0" (spec D6). This is
% only meaningful if the propagated state's initial t is 0 -- true for every
% seed built from mee_seed (t0 = 0 by construction) but NOT independently
% asserted anywhere in the bridge/solve pipeline.
%
% REFERENCES:
%   [1] spec 2026-07-22-elliptic-geo-cr3bp-phase0-design.md sec 3 (constants:
%       mu_M = 4902.800 km^3/s^2, D_EM = 384400 km; ratio 0.0123 consistent
%       with CR3BP mu* = 0.0121506 via mu*/(1-mu*)).
if nargin < 2 || isempty(phi0), phi0 = 0; end
if nargin < 3 || isempty(gain), gain = 1; end
muM_km3s2 = 4902.800;
DM_km     = 384400;
pert.muM  = (muM_km3s2 / par.muKm3s2) * par.mu;          % canonical GM_moon
pert.DM   = DM_km / par.LU_km;                           % canonical distance
pert.nM   = sqrt((par.mu + pert.muM) / pert.DM^3);       % canonical rate
pert.phi0 = phi0;
pert.gain = gain;
end
