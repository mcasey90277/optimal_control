function psi0 = greatCircleBearing(lat1,lon1,lat2,lon2)
%% Purpose:
%
%  Initial bearing of the great-circle arc running from point 1 to point 2 on
%  a sphere, measured CLOCKWISE FROM NORTH and wrapped into [0,2*pi). That is
%  exactly the convention of this library's heading state psi, so the result
%  drops straight into a state vector without a sign flip or an offset.
%  Companion to coorbital.util.greatCircle, which gives the central angle
%  along the same arc; together they turn a launch point and a destination
%  into a launch azimuth and a required range. Elementwise, so any
%  broadcast-compatible mix of scalars and [N x 1] columns is accepted.
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
%  psi0             [N x 1]                     Initial bearing at point 1
%                                               (rad), clockwise from north,
%                                               in [0,2*pi)
%
%% Notes:
%
%  INITIAL bearing, not a constant course. A great-circle arc is not a rhumb
%  line: its bearing changes continuously along the track, and the value
%  returned here is the one at point 1 only. Flying a constant heading of
%  psi0 does NOT arrive at point 2.
%
%  Bearing is NOT symmetric under exchange of the two points, and this is the
%  single most common way to get it wrong. Distance is symmetric --
%  greatCircle(A,B) equals greatCircle(B,A) -- but the arc meets the local
%  meridian at a different angle at each end, so bearing(B,A) is generally
%  NOT bearing(A,B) turned through pi. It is only so for a pure meridian leg
%  or a leg along the equator. Worked example, from Los Angeles at 33.94 N
%  118.41 W to New York at 40.64 N 73.78 W: the outbound course is 65.867 deg
%  and the return course is 273.841 deg, whereas reversing the outbound
%  course would give 245.867 deg. The 27.974 deg discrepancy is the whole
%  point. An implementation with lat1 and lat2 transposed reproduces all four
%  cardinal directions correctly and only shows itself on a case like this.
%
%  Coincident points have no arc and no course. atan2(0,0) is defined as
%  zero, so a point to itself returns zero rather than NaN. That keeps a
%  degenerate user input from poisoning a state vector, but zero is a
%  convention here, not an answer.
%
%  The poles are the other degenerate case. At exactly the north pole every
%  direction is south and the returned value is whatever the limit of the
%  atan2 expression happens to be along the approach; it is not meaningful.
%
%  This is a SPHERICAL bearing. On the WGS-84 ellipsoid the geodesic azimuth
%  differs by up to a few tenths of a degree at continental range. That is
%  consistent with the spherical Earth used throughout this library.
%
%% References:
%   [1] Bowditch, N., "The American Practical Navigator," Pub. No. 9,
%       National Geospatial-Intelligence Agency, Chapter on Great-Circle
%       Sailing. The initial course of the great-circle track is obtained
%       from the spherical triangle formed by the pole and the two points;
%       applying the four-parts formula to that triangle and clearing the
%       cotangents gives the atan2 form used below.
%   [2] Sinnott, R.W., "Virtues of the Haversine," Sky and Telescope, Vol. 68,
%       No. 2, 1984, p. 159. Cited for the companion distance routine
%       coorbital.util.greatCircle, whose argument order this function
%       matches exactly.
%
%  The atan2 form is preferred over the arccosine or arcsine forms of the
%  same course angle for two reasons. It resolves the full circle without a
%  separate east/west quadrant test, which is where hand-rolled bearing code
%  usually breaks, and it keeps its relative accuracy for short legs, where
%  both the numerator and the denominator go to zero together and their ratio
%  stays well conditioned.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
              latL = deg2rad( 33.94);   lonL = deg2rad(-118.41);
              latJ = deg2rad( 40.64);   lonJ = deg2rad( -73.78);
             psiLJ = coorbital.util.greatCircleBearing(latL,lonL,latJ,lonJ);
             psiJL = coorbital.util.greatCircleBearing(latJ,lonJ,latL,lonL);
              rngM = c.rE.*coorbital.util.greatCircle(latL,lonL,latJ,lonJ);
    fprintf('LAX -> JFK: initial bearing %.4f deg, range %.1f km\n', ...
            rad2deg(psiLJ),rngM./1000);
    fprintf('JFK -> LAX: initial bearing %.4f deg\n',rad2deg(psiJL));
    fprintf('reversing the outbound course would give %.4f deg, which is\n', ...
            rad2deg(mod(psiLJ + pi,2*pi)));
    fprintf('wrong by %.4f deg: bearing is not symmetric.\n', ...
            rad2deg(abs(mod(psiJL - psiLJ,2*pi) - pi)));
              psi0 = [];
    return;
end

%% Longitude difference. Only the difference enters, so a common longitude
%% offset on both points cannot shift the answer:
              dLon = lon2 - lon1;

%% Numerator and denominator of the spherical four-parts course expression.
%% Both vanish together as point 2 approaches point 1, which is what lets
%% atan2 stay well conditioned for short legs. yB is the eastward component
%% of the initial track direction and xB the northward one, both scaled by
%% the same positive factor, so their ratio is the tangent of the course:
                yB = sin(dLon).*cos(lat2);
                xB = cos(lat1).*sin(lat2) - sin(lat1).*cos(lat2).*cos(dLon);

%% atan2 returns (-pi,pi] measured from north, positive toward the east,
%% which is already the clockwise-from-north sense. Only the negative half
%% needs moving:
              psi0 = atan2(yB,xB);

%% Wrap into [0,2*pi). mod alone is not enough: for a bearing a hair west of
%% north, say -1e-18 rad, the sum 2*pi - 1e-18 rounds to exactly 2*pi in
%% double precision, and mod hands that back. The half-open contract matters
%% because a caller may index or bin on it, so the upper end is folded to
%% zero explicitly:
              psi0 = mod(psi0,2.*pi);
           atTwoPi = psi0 >= 2.*pi;
     psi0(atTwoPi) = 0;
end
