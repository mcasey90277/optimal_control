function P = booster_params()
% BOOSTER_PARAMS  Falcon-9-class 3-DOF landing-burn parameter set.
%
% Single source of truth for physical constants, boundary conditions,
% solver grid and Monte-Carlo settings. Public F9 estimates; every number
% adjustable here and only here.
%
% INPUTS:  none
% OUTPUTS: P - parameter struct (fields documented inline below)
%
% REFERENCES:
%   [1] Blackmore, Acikmese, Scharf, "Minimum-Landing-Error Powered-Descent
%       Guidance for Mars Landing Using Convex Optimization," JGCD 2010.
%   [2] docs/superpowers/specs/2026-08-08-booster-landing-design.md

%% Vehicle (public Falcon 9 block-5 estimates):
P.mdry   = 25600;              % dry mass [kg]
P.m0     = 30000;              % mass at landing-burn start [kg]
P.Tmax   = 845e3;              % one Merlin 1D, sea level [N]
P.Tmin   = 0.40 * P.Tmax;      % ~40 percent min throttle [N]
P.Isp    = 282;                % sea-level Isp [s]
P.g0     = 9.80665;            % standard gravity [m/s^2]
P.gvec   = [0; 0; -P.g0];      % flat-Earth gravity, z up [m/s^2]

%% Boundary conditions (pad at origin, z up):
P.r0     = [500; 100; 2000];   % post-entry-burn position [m]
P.v0     = [-30; 0; -180];     % descending ~180 m/s [m/s]

%% Path constraints:
P.gs_deg        = 30;          % glideslope min elevation angle [deg]
P.theta_max_deg = Inf;         % thrust-pointing cone half-angle, Inf = off

%% Discretization / solver:
P.N      = 60;                 % Hermite-Simpson segments (collocation)
P.Nconv  = 120;                % trapezoid nodes for the convex solver
P.tf_lo  = 10;                 % free-final-time bracket [s]
P.tf_hi  = 50;

%% Atmosphere (Phase 2, OFF by default -- vacuum keeps convexification exact):
P.drag.on   = false;
P.drag.rho0 = 1.225;           % sea-level density [kg/m^3]
P.drag.H    = 8500;            % scale height [m]
P.drag.Cd   = 1.0;             % landing-leg config drag coefficient [-]
P.drag.A    = 10.75;           % reference area, 3.7 m diameter [m^2]

%% Success criteria / Monte Carlo:
P.pad_radius = 15;             % landing accuracy requirement [m]
P.vtd_max    = 2.0;            % max touchdown speed [m/s]
P.seed       = 42;             % rng seed, all random draws
end
