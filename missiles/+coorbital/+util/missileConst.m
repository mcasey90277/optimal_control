function c = missileConst()
%% Purpose:
%
%  Return the Earth and standard-air constants used across the missile
%  trajectory library. Single source of truth -- no routine may hard-code
%  these values.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  c                Struct                      Constants, SI units:
%                                               muE    (m^3/s^2)
%                                               rE     (m) equatorial radius
%                                               omegaE (rad/s) Earth rotation
%                                               rho0   (kg/m^3) sea-level density
%                                               Hscale (m) density scale height
%                                               T0     (K) isothermal temperature
%                                               Rair   (J/kg/K) specific gas constant
%                                               gamAir (-) ratio of specific heats
%                                               g0     (m/s^2) standard gravity
%
%% References:
%   [1] WGS-84, NIMA TR8350.2, 3rd ed., 1997. (muE, rE, omegaE)
%   [2] U.S. Standard Atmosphere, 1976. (rho0, Rair, gamAir, g0)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Earth:
             c.muE = 3.986004418e14;    %m^3/s^2, WGS-84
              c.rE = 6378137;           %m, WGS-84 equatorial radius
          c.omegaE = 7.292115e-5;       %rad/s, WGS-84

%% Standard air:
            c.rho0 = 1.225;             %kg/m^3, sea level, US76
          c.Hscale = 7200;              %m, exponential density scale height
              c.T0 = 250;               %K, isothermal reference (see expAtmos)
            c.Rair = 287.0528;          %J/kg/K, specific gas constant for air
          c.gamAir = 1.4;               %ratio of specific heats
              c.g0 = 9.80665;           %m/s^2, standard gravity
end
