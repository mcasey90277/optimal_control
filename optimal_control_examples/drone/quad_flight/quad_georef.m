function [lat, lon, alt] = quad_georef(pLib, lla0)
%% Purpose:
%
%   Georeference a flight: the library's world positions (local-level ENU
%   at a user-defined origin) to WGS-84 latitude, longitude and altitude,
%   through ECEF. The dynamics stay flat-earth, which is right at drone
%   scale; this is the labelling layer only.
%
%   Self-contained on purpose: on this machine the Aerospace Toolbox's
%   lla2ecef/ecef2lla (degrees, M-by-3) shadow mlabTools' (radians,
%   separate arguments) and mlabTools' ned2ecef shadows MATLAB's, so a call
%   by name resolves to a different signature depending on path order.
%   The WGS-84 forward map and the Bowring inverse are a dozen lines each.
%
%% Inputs:
%
%  pLib                     [3 x N]                 ENU positions (m) from
%                                                   the origin
%
%  lla0                     [1 x 3]                 Origin [lat deg, lon
%                                                   deg, alt m (WGS-84)]
%
%% Outputs:
%
%  lat                      [1 x N]                 Geodetic latitude (deg)
%
%  lon                      [1 x N]                 Longitude (deg)
%
%  alt                      [1 x N]                 Altitude above the
%                                                   ellipsoid (m)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: 100 m north and 20 m up from an origin:
     [lat, lon, alt] = quad_georef([0; 100; 20], [33.77, -118.39, 50]);
     fprintf('100 m N, 20 m up of (33.77, -118.39, 50): %.6f deg, %.6f deg, %.2f m\n', lat, lon, alt);
     return
end

if size(pLib, 1) ~= 3
    error('quad_georef:shape', 'pLib must be [3 x N]');
end

a  = 6378137.0;                                      % WGS-84 semi-major axis (m)
f  = 1/298.257223563;                                % flattening
e2 = f*(2 - f);                                      % first eccentricity squared

phi0 = lla0(1)*pi/180;  lam0 = lla0(2)*pi/180;  h0 = lla0(3);
sp = sin(phi0);  cp = cos(phi0);  sl = sin(lam0);  cl = cos(lam0);
N0 = a/sqrt(1 - e2*sp^2);
r0 = [(N0 + h0)*cp*cl; (N0 + h0)*cp*sl; (N0*(1 - e2) + h0)*sp];

% ENU basis at the origin, as ECEF columns [east north up]
Ren = [-sl, -sp*cl, cp*cl;
        cl, -sp*sl, cp*sl;
         0,     cp,    sp];
r = r0 + Ren*pLib;

% Bowring's inverse (one pass is good to millimetres; two is overkill)
x = r(1, :);  y = r(2, :);  z = r(3, :);
b  = a*(1 - f);
ep2 = (a^2 - b^2)/b^2;
rho = sqrt(x.^2 + y.^2);
th  = atan2(z*a, rho*b);
phi = atan2(z + ep2*b*sin(th).^3, rho - e2*a*cos(th).^3);
for pass = 1:2
    N   = a./sqrt(1 - e2*sin(phi).^2);
    h   = rho./cos(phi) - N;
    phi = atan2(z, rho.*(1 - e2*N./(N + h)));
end
N   = a./sqrt(1 - e2*sin(phi).^2);
alt = rho./cos(phi) - N;
lat = phi*180/pi;
lon = atan2(y, x)*180/pi;
end
