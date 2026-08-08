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
%                                               in [0,2*pi). Raises
%                                               coorbital:greatCircleBearing:degenerateArc
%                                               for coincident or antipodal
%                                               points and
%                                               coorbital:greatCircleBearing:polarOrigin
%                                               at a pole, see Notes
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
%  point.
%
%  That discrepancy is NOT, however, a way to detect a lat1/lat2
%  transposition, and it is worth being precise about why, because the
%  natural assumption is the opposite. Exchanging lat1 with lat2 inside the
%  formula moves BOTH courses -- LAX-JFK becomes 86.159 deg and JFK-LAX
%  becomes 294.133 deg -- and leaves their difference at exactly the same
%  27.974 deg. What a transposition does break is the meridional cardinals
%  and the absolute courses: due north comes back as pi and due south as
%  zero, the two east-west cardinals survive, and every general-position
%  course is wrong. So the things to assert against are absolute bearings on
%  a near-meridional leg, not the reverse-course offset.
%
%% THE DEGENERATE CASES ARE REFUSED, NOT ANSWERED:
%
%  This function is a targeting azimuth. Where the azimuth does not exist, a
%  plausible-looking number is worse than no number at all: it flies. So the
%  three degenerate geometries raise rather than return.
%
%  Write the two components as yB, the eastward part of the initial track
%  direction at point 1, and xB, the northward part. Both are scaled by the
%  same positive factor, and that factor is sin(Delta) for a central angle
%  Delta, so hypot(xB,yB) IS the sine of the angular separation. One test on
%  it therefore catches both ends at once:
%
%      COINCIDENT points have no arc and no course, and hypot is exactly
%      zero. atan2(0,0) is defined as zero, so the old behaviour was to
%      return due north for a point to itself -- a convention, not an answer,
%      and one that puts a launch heading of 0 deg into a state vector.
%
%      ANTIPODAL points are degenerate the opposite way: EVERY great circle
%      joins them, so there is no unique initial course, and hypot is again
%      zero to rounding. From (0,0) to (0,pi) the old form returned exactly
%      90 deg, because sin(pi) evaluates to 1.22e-16 rather than zero against
%      a true-zero denominator. Nudge either latitude by 1e-9 rad and the
%      answer swings the whole way across: +1e-9 gives 7.02e-6 deg and -1e-9
%      gives 179.99999 deg. Near-antipodal geometry stays catastrophically
%      ill-conditioned well outside the refusal band, and no tolerance fixes
%      that -- the refusal only removes the cases with no answer at all.
%
%      AT A POLE the local north direction that the bearing is measured from
%      does not exist. hypot does NOT vanish there, so this one is tested
%      separately, on cos(lat1).
%
%  The tolerance is 1e-12 on both tests. On hypot that is sin(Delta) <= 1e-12,
%  i.e. points within about 6.4 micrometres of coincidence or of the antipode
%  at Earth radius; on cos(lat1) it is a latitude within 1e-12 rad -- 6.4
%  micrometres -- of a pole. Anything a targeting problem can pose is far
%  outside it, and anything inside it has no azimuth to report.
%
%  The identifiers are coorbital:greatCircleBearing:degenerateArc and
%  coorbital:greatCircleBearing:polarOrigin. Elementwise inputs are refused as
%  a batch: if ANY element is degenerate the whole call raises, because a
%  column with one silently wrong azimuth in it is the failure this guard
%  exists to prevent.
%
%  The companion coorbital.util.greatCircle degrades gracefully at all three
%  and returns a meaningful distance; this function cannot, and does not try.
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
%  Michael Casey  Refuse coincident, antipodal and polar inputs  08/07/2026
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

%% Refuse the geometries that have no azimuth, rather than return the
%% rounding artefact atan2 makes of them. See THE DEGENERATE CASES in the
%% header: hypot(xB,yB) is the sine of the angular separation, so one test on
%% it catches coincident AND antipodal points, and the pole -- where the
%% north direction the bearing is measured FROM does not exist -- needs its
%% own test because hypot does not vanish there:
            tolDeg = 1e-12;
              magB = hypot(xB,yB);
if any(magB(:) <= tolDeg)
    error('coorbital:greatCircleBearing:degenerateArc', ...
        ['There is no initial bearing between these points: the track ' ...
         'components vanish to %.3g, at or below the %.3g tolerance, which ' ...
         'means the two points are coincident or antipodal to within about ' ...
         '6.4 micrometres of arc. Every great circle joins an antipodal ' ...
         'pair and none joins a point to itself, so any value returned ' ...
         'here would be an artefact of rounding, and it would fly.'], ...
        min(magB(:)),tolDeg);
end
if any(abs(cos(lat1(:))) <= tolDeg)
    error('coorbital:greatCircleBearing:polarOrigin', ...
        ['A bearing is measured clockwise from local north, and at a pole ' ...
         'there is no local north: |cos(lat1)| came to %.3g, at or below ' ...
         'the %.3g tolerance. Move the origin off the pole.'], ...
        min(abs(cos(lat1(:)))),tolDeg);
end

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
