function d = greatCircle(lat1,lon1,lat2,lon2)
%% Purpose:
%
%  Central angle between two points on a sphere, via the haversine form,
%  which stays accurate for small separations where the spherical law of
%  cosines loses precision. Multiply the result by a radius to obtain surface
%  arc length. Elementwise, so any broadcast-compatible mix of scalars and
%  [N x 1] columns is accepted.
%
%% Inputs:
%
%  lat1             [N x 1]                     Latitude of point 1 (rad)
%
%  lon1             [N x 1]                     Longitude of point 1 (rad)
%
%  lat2             [N x 1]                     Latitude of point 2 (rad)
%
%  lon2             [N x 1]                     Longitude of point 2 (rad)
%
%% Outputs:
%
%  d                [N x 1]                     Central angle (rad); multiply
%                                               by a radius for arc length
%
%% References:
%   [1] Sinnott, R.W., "Virtues of the Haversine," Sky and Telescope, Vol. 68,
%       No. 2, 1984, p. 159. The haversine form is used here in preference to
%       the spherical law of cosines,
%           cos(d) = sin(lat1)sin(lat2) + cos(lat1)cos(lat2)cos(lon2-lon1),
%       because for small d that expression evaluates to 1 - d^2/2 + O(d^4).
%       Forming it in double precision and then taking acos throws away about
%       half the significant digits, leaving an absolute error of order
%       sqrt(eps) ~ 1e-8 rad; below d ~ 1e-8 rad it returns exactly zero.
%       Haversine works with half-angle sines throughout and keeps full
%       relative accuracy all the way down to d = 0, which matters here
%       because terminal-phase ranges are short. Its own weak point is the
%       near-antipodal case, where it degrades to about sqrt(eps) -- harmless
%       for a ballistic or glide range calculation.
%   [2] Vincenty, T., "Direct and Inverse Solutions of Geodesics on the
%       Ellipsoid with Application of Nested Equations," Survey Review,
%       Vol. 23, No. 176, 1975, pp. 88-93. The atan2 form given there is
%       accurate at both ends; it is not used because it costs more and the
%       antipodal end is not exercised by this library.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
             dDemo = coorbital.util.greatCircle(0,0,0,deg2rad(1));
    fprintf('one degree of equator: %.15f rad, %.4f m of arc\n', ...
            dDemo,c.rE.*dDemo);
             dTiny = coorbital.util.greatCircle(0,0,1e-10,0);
    fprintf('1e-10 rad separation:  haversine %.6e rad, cosine form %.6e rad\n', ...
            dTiny,acos(cos(1e-10)));
                 d = [];
    return;
end

%% Half-angle differences. Every quantity below stays small when the points
%% are close, so nothing is ever subtracted from a number near one:
              dLat = lat2 - lat1;
              dLon = lon2 - lon1;

%% a is the haversine of the central angle, hav(d) = sin^2(d/2). Clamped at
%% one because roundoff can push it a hair above for near-antipodal points,
%% where asin would otherwise return a complex value:
                 a = sin(dLat./2).^2 + cos(lat1).*cos(lat2).*sin(dLon./2).^2;
                 d = 2.*asin(min(1,sqrt(a)));
end
