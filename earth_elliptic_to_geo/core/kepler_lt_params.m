function p = kepler_lt_params(thrustN, m0kg, ispS)
% KEPLER_LT_PARAMS  Constants + canonical units, Earth 2-body low-thrust problem.
%
% Nondimensionalization: LU = GEO radius 42165 km, TU = sqrt(LU^3/mu_earth), so
% mu = 1, GEO circular speed = 1, GEO period = 2*pi. Mass unit = m0kg.
%
% INPUTS:
%   thrustN - max thrust [N]        (paper cases: 10, 5, 2.5, 1)
%   m0kg    - initial mass [kg]     (paper: 1500)
%   ispS    - specific impulse [s]  (default 2000; the benchmark's exact value
%             is 1994.8 s per Caillau & Noailles 2001 -- see the note at ispS below)
%
% OUTPUTS:
%   p - struct: dimensional anchors .g0 .muKm3s2 .LU_km .TU_s .VU_kms .AU_ms2;
%       nondim .mu(=1) .Tmax (thrust @ m=1) .c (exhaust velocity) .pSund;
%       echo .thrustN .m0kg .ispS
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004 (problem constants).
%   [2] earth_elliptic_to_geo/process/DESIGN.md sec 2 (units decision).
% ISP SOURCE (process/DESIGN.md open item 1 -- CLOSED 2026-07-19): the paper
% (Gergaud-Haberkorn-Martinon, JGCD 2004) never states Isp numerically; its
% benchmark constants come from Caillau & Noailles, "Coplanar control of a
% satellite around the Earth," ESAIM COCV 6 (2001), p.255
% (min_fuel_papers/COCV_2001__6__239_0.pdf), which gives the mass-flow
% coefficient delta = 0.05112 km^-1 s in mdot = -delta*|thrust|. Since
% delta = 1/c = 1/(Isp*g0):  c = 1/delta = 19.562 km/s  =>  Isp = c/g0 = 1994.8 s.
% So Caillau & Noailles use Isp = 1994.8 s. We keep 2000 s as the default (only
% 0.27% above the exact 1994.8 s -> absolute masses run ~0.3 kg high; trajectory
% STRUCTURE is Isp-independent). Pass ispS = 1994.8 for the exact source value.
if nargin < 3, ispS = 2000; end
p.g0      = 9.80665;                          % [m/s^2]
p.muKm3s2 = 398600.47;                        % [km^3/s^2]
p.LU_km   = 42165;                            % GEO radius = terminal P [km]
p.TU_s    = sqrt(p.LU_km^3 / p.muKm3s2);      % => mu = 1
p.VU_kms  = p.LU_km / p.TU_s;
p.AU_ms2  = 1000 * p.VU_kms / p.TU_s;         % acceleration unit [m/s^2]
p.mu      = 1;
p.thrustN = thrustN;  p.m0kg = m0kg;  p.ispS = ispS;
p.Tmax    = (thrustN/m0kg) / p.AU_ms2;        % nondim thrust at m = 1
p.c       = (ispS*p.g0/1000) / p.VU_kms;      % nondim exhaust velocity
p.pSund   = 1.5;                              % Sundman power, dt/dtau = r^pSund
end
