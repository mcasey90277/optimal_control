function [rho,P,T,a] = expAtmos(h)
%% Purpose:
%
%  Isothermal exponential atmosphere. Density decays with a single scale
%  height and temperature is held constant, which makes the model crude above
%  roughly 30 km but gives closed-form entry solutions to validate against.
%  Interface-compatible with richer models (us76Atmos, nrlmsise) so the
%  environment struct can swap them without touching the equations of motion.
%
%% Inputs:
%
%  h                [N x 1]                     Geometric altitude above the
%                                               reference sphere (m)
%
%% Outputs:
%
%  rho              [N x 1]                     Density (kg/m^3)
%
%  P                [N x 1]                     Pressure (Pa)
%
%  T                [N x 1]                     Temperature (K)
%
%  a                [N x 1]                     Speed of sound (m/s)
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Ch. 2.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
              hVec = linspace(0,120e3,400)';
    [rhoV,~,~,~] = coorbital.atmos.expAtmos(hVec);
    figure('color',[1 1 1]);
    semilogx(rhoV,hVec./1000,'linewidth',1.5); grid on;
    xlabel('density (kg/m^3)'); ylabel('altitude (km)');
    title('Isothermal exponential atmosphere');
    [rho,P,T,a] = deal([]);
    return;
end

                 c = coorbital.util.missileConst();

%% Density decays one e-fold per scale height:
               rho = c.rho0.*exp(-h./c.Hscale);

%% Isothermal by construction:
                 T = c.T0.*ones(size(h));

%% Pressure from the ideal gas law, sound speed from the same temperature:
                 P = rho.*c.Rair.*T;
                 a = sqrt(c.gamAir.*c.Rair.*T);
end
