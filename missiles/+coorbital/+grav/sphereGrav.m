function [gr,gLat] = sphereGrav(r,lat)
%% Purpose:
%
%  Point-mass gravitational acceleration. Returns two components so that
%  oblate models (j2Grav, j24Grav) satisfy the same interface -- J2 produces a
%  latitudinal acceleration that a single "local up" scalar cannot represent.
%
%% Inputs:
%
%  r                [N x 1]                     Geocentric radius (m)
%
%  lat              [N x 1]                     Geocentric latitude (rad).
%                                               Unused here; present so the
%                                               signature matches j2Grav.
%
%% Outputs:
%
%  gr               [N x 1]                     Acceleration toward the centre,
%                                               positive downward (m/s^2)
%
%  gLat             [N x 1]                     Northward acceleration (m/s^2);
%                                               identically zero for a point mass
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
              rVec = c.rE + linspace(0,200e3,200)';
       [grV,~] = coorbital.grav.sphereGrav(rVec,zeros(size(rVec)));
    figure('color',[1 1 1]);
    plot((rVec-c.rE)./1000,grV,'linewidth',1.5); grid on;
    xlabel('altitude (km)'); ylabel('g (m/s^2)');
    title('Point-mass gravity');
    [gr,gLat] = deal([]);
    return;
end

                 c = coorbital.util.missileConst();

%% Inverse-square magnitude, positive toward the centre:
                gr = c.muE./r.^2;

%% A point mass has no latitudinal component:
              gLat = zeros(size(r));
end
