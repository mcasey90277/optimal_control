function test_greatCircle()
%% Purpose:
%
%  Verify the haversine central angle against cases with exact or independent
%  closed-form answers, confirm it stays accurate at small separations where
%  the cosine form degrades, confirm the argument order really is
%  (lat1,lon1,lat2,lon2) and not (lon,lat), and confirm the surface arc length
%  it produces when scaled by the Earth radius. Covers the spec's
%  great-circle range check.
%
%  Reference values are either exact by symmetry (quarter turns, antipodes),
%  hand-computed literals (arc lengths from the WGS-84 radius), or the
%  SPHERICAL LAW OF COSINES -- a genuinely different formula from the one
%  under test. The cosine form is used as the reference only for separations
%  of order a degree or more, where it is accurate to near machine precision.
%  It is deliberately NOT used as the reference for the small-separation
%  cases, since that is precisely the regime in which it is the wrong one.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% A quarter turn along the equator is pi/2:
                 d = coorbital.util.greatCircle(0,0,0,pi/2);
    assert(abs(d - pi/2) < 1e-12,'equatorial quarter turn wrong: got %.15f',d);

%% Equator to pole is pi/2:
                 d = coorbital.util.greatCircle(0,0,pi/2,0);
    assert(abs(d - pi/2) < 1e-12,'equator-to-pole wrong: got %.15f',d);

%% Coincident points give exactly zero -- not NaN, and not a cancellation
%% residue of order 1e-8 as the cosine form would leave:
                 d = coorbital.util.greatCircle(0.3,1.1,0.3,1.1);
    assert(~isnan(d),'coincident points returned NaN');
    assert(d == 0,'coincident points must give exactly zero, got %g',d);

%% Antipodal points give pi, both across the equator and pole to pole:
                 d = coorbital.util.greatCircle(0,0,0,pi);
    assert(abs(d - pi) < 1e-9,'antipodal separation wrong: got %.15f',d);
                 d = coorbital.util.greatCircle(-pi/2,0,pi/2,0);
    assert(abs(d - pi) < 1e-9,'pole-to-pole separation wrong: got %.15f',d);

%% One degree of equator is rE*pi/180 = 111319.4908 m of surface arc, worked
%% out by hand from the WGS-84 radius. Hard-coded rather than recomputed from
%% c.rE so that a wrong radius in the constants file cannot cancel itself out
%% of both sides of the comparison:
              arcM = c.rE.*coorbital.util.greatCircle(0,0,0,deg2rad(1));
    assert(abs(arcM - 111319.4908) < 1e-3, ...
        'one degree of equator wrong: got %.4f m, expected 111319.4908 m',arcM);

%% Equator to pole is a quarter of the circumference, rE*pi/2 = 10018754.1714 m:
             quadM = c.rE.*coorbital.util.greatCircle(0,0,pi/2,0);
    assert(abs(quadM - 10018754.1714) < 1e-3, ...
        'quarter circumference wrong: got %.4f m, expected 10018754.1714 m',quadM);

%% Argument order is (lat1,lon1,lat2,lon2). Two points 90 deg of LONGITUDE
%% apart at 60 deg N are acos(sin^2(60) + cos^2(60)cos(90)) = acos(0.75)
%% = 0.7227342478 rad apart, NOT pi/2 -- and pi/2 is exactly what the same
%% call returns if latitude and longitude are read in the opposite order:
             dLonQ = coorbital.util.greatCircle(deg2rad(60),0,deg2rad(60),pi/2);
              dRef = acos(sind(60).^2 + cosd(60).^2.*cos(pi/2));
    assert(abs(dLonQ - dRef) < 1e-12, ...
        '90 deg of longitude at 60 N wrong: got %.15f, expected %.15f',dLonQ,dRef);
    assert(abs(dLonQ - pi/2) > 0.5, ...
        'result equals the pure-longitude answer: lat and lon look swapped');

%% Same point made a second way: one degree of LATITUDE is one degree of arc
%% at any longitude, while one degree of LONGITUDE at 60 N is only about half
%% a degree of arc. A lat/lon transposition collapses that distinction:
             dLat1 = coorbital.util.greatCircle(deg2rad(60),0,deg2rad(61),0);
             dLon1 = coorbital.util.greatCircle(deg2rad(60),0,deg2rad(60),deg2rad(1));
    assert(abs(dLat1 - deg2rad(1)) < 1e-14, ...
        'one degree of latitude must be one degree of arc: got %.15f',dLat1);
    assert(abs(dLon1 - 0.5.*deg2rad(1)) < 1e-6, ...
        'one degree of longitude at 60 N should be ~half a degree of arc: got %.15f',dLon1);

%% A general asymmetric pair against the cosine form, which at this
%% separation (order 1 rad) is an accurate independent reference:
              latA = 0.4;    lonA = -1.2;
              latB = -0.7;   lonB =  2.3;
               dAB = coorbital.util.greatCircle(latA,lonA,latB,lonB);
             dABRf = acos(sin(latA).*sin(latB) + ...
                          cos(latA).*cos(latB).*cos(lonB - lonA));
    assert(abs(dAB - dABRf) < 1e-12, ...
        'general pair disagrees with the cosine form: got %.15f, expected %.15f', ...
        dAB,dABRf);

%% Symmetric in its arguments:
               dBA = coorbital.util.greatCircle(latB,lonB,latA,lonA);
    assert(abs(dAB - dBA) < 1e-15, ...
        'distance must be symmetric: %.17f vs %.17f',dAB,dBA);
    assert(dAB > 0,'distinct points must have a positive separation');

%% Small separations stay accurate: one arcsecond of latitude. The cosine
%% form misses this by order 1e-11 rad:
            dSmall = coorbital.util.greatCircle(0,0,arcsecToRad(1),0);
    assert(abs(dSmall - arcsecToRad(1)) < 1e-15, ...
        'small-angle accuracy lost: got %.17e',dSmall);

%% Smaller still, in the regime where the cosine form underflows to exactly
%% zero. The second assert documents that failure mode rather than assuming it:
             dTiny = coorbital.util.greatCircle(0,0,1e-10,0);
    assert(abs(dTiny - 1e-10) < 1e-22, ...
        'sub-nanoradian accuracy lost: got %.17e',dTiny);
    assert(acos(cos(1e-10)) == 0, ...
        'reference check: the cosine form was expected to underflow to zero here');

%% Elementwise over columns, and consistent with the scalar calls:
              latV = [0; 0.1; -0.4; 1.0];
              lonV = [0; 0.2;  0.9; 3.0];
                dV = coorbital.util.greatCircle(latV,lonV,latV,lonV);
    assert(isequal(size(dV),size(latV)),'output must match the input shape');
    assert(all(dV == 0),'self-distance must be exactly zero elementwise');
               dV2 = coorbital.util.greatCircle(latV,lonV,0,0);
    assert(isequal(size(dV2),size(latV)),'scalar broadcast must preserve shape');
    for k = 1:numel(latV)
              dScl = coorbital.util.greatCircle(latV(k),lonV(k),0,0);
        assert(abs(dV2(k) - dScl) < 1e-15, ...
            'vectorized result disagrees with the scalar call at element %d',k);
    end
end

function r = arcsecToRad(a)
%% Purpose:
%
%  Convert arcseconds to radians. Written out rather than folded into the
%  caller so that the small-separation cases read as the angles they are.
%
%% Inputs:
%
%  a                [N x 1]                     Angle (arcsec)
%
%% Outputs:
%
%  r                [N x 1]                     The same angle (rad)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------
                 r = a.*pi./(180*3600);
end
