function test_greatCircleBearing()
%% Purpose:
%
%  Verify the great-circle initial bearing against cases whose answers are
%  fixed by symmetry alone, against one hand-checked real-world pair, and
%  against the two properties that a wrong implementation most easily
%  violates: the wrap into [0,2*pi) and the ASYMMETRY of bearing under
%  argument exchange.
%
%  Every reference value below is either exact by symmetry -- a due-east leg
%  along the equator is pi/2 whatever formula produced it -- or a literal
%  worked out independently of the function under test. Nothing here
%  re-evaluates the atan2 expression that greatCircleBearing itself
%  evaluates, so a sign slip or a transposition inside that expression
%  cannot cancel out of both sides of a comparison.
%
%  The asymmetry block is the load-bearing one. Distance is symmetric:
%  greatCircle(A,B) equals greatCircle(B,A), and test_greatCircle asserts
%  exactly that. Bearing is NOT. On a sphere the two ends of a great-circle
%  arc generally meet the local meridian at different angles, so
%  bearing(B,A) is not bearing(A,B) turned through pi except in the special
%  cases of a pure meridian leg or an equatorial leg. A reader who assumes
%  otherwise, or a caller that reuses one bearing for the return leg, is
%  wrong by 27.974 deg on the LAX-JFK pair -- about 1900 km of cross-range
%  over a 3979 km leg.
%
%  One caution about that offset, which was measured rather than assumed and
%  is recorded here because it is easy to get backwards: the 27.974 deg
%  offset is NOT a transposition detector. Exchanging lat1 with lat2 inside
%  the formula moves both courses and leaves their difference unchanged at
%  27.974 deg. What catches a transposition is the pair of ABSOLUTE course
%  values, and the near-meridional case further down, where a northbound leg
%  transposes into a southbound one.
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
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Due east along the equator. Both points sit on the equator, so the arc IS
%% the equator and it crosses every meridian at a right angle: the bearing is
%% pi/2 at every point of it, exactly, with no rounding to allow for:
               psi = coorbital.util.greatCircleBearing(0,0,0,deg2rad(30));
    assert(psi == pi/2, ...
        'due east on the equator must be exactly pi/2: got %.17f',psi);

%% Due north along the prime meridian. The arc is a meridian and north is by
%% definition zero on a clockwise-from-north convention:
               psi = coorbital.util.greatCircleBearing(0,0,deg2rad(10),0);
    assert(psi == 0, ...
        'due north must be exactly 0: got %.17f',psi);

%% Due south along the same meridian, which is where the clockwise-from-north
%% convention gets pinned down: south is pi, not -pi and not 3*pi/2. A
%% counter-clockwise or from-east convention would not land here:
               psi = coorbital.util.greatCircleBearing(deg2rad(10),0,0,0);
    assert(psi == pi, ...
        'due south must be exactly pi: got %.17f',psi);

%% Due west along the equator. The raw atan2 gives -pi/2; the wrap must carry
%% it to 3*pi/2 = 4.712388980384690 rad. A function that skipped the wrap
%% returns a negative number here and fails outright:
               psi = coorbital.util.greatCircleBearing(0,0,0,deg2rad(-30));
    assert(psi == 3*pi/2, ...
        'due west must wrap to exactly 3*pi/2: got %.17f',psi);
    assert(psi >= 0, ...
        'due west returned a negative bearing: the wrap is missing (%.17f)',psi);

%% LAX to JFK. LAX is 33.94 N, 118.41 W and JFK is 40.64 N, 73.78 W.
%%
%% The reference value 65.867 deg is quoted in the implementation plan
%% docs/plan_2026-08-07_targeting_and_viz.md, where it was obtained by
%% evaluating the spherical initial-course formula at these coordinates and
%% cross-checking it against the published great-circle course for the
%% LAX-JFK leg, which is about 66 deg over a distance of about 3970 km.
%% That leg is the standard worked example in Ed Williams' Aviation
%% Formulary. It was re-derived independently here, outside this library,
%% before this test was written and agreed to 65.8669298403 deg. The 0.01
%% deg window is far tighter than the roughly 28 deg error that an argument
%% transposition produces and far tighter than the 0.4 deg spread between
%% published spherical and ellipsoidal courses for the leg, so it pins the
%% formula without pretending to geodetic accuracy:
              latL = deg2rad( 33.94);   lonL = deg2rad(-118.41);
              latJ = deg2rad( 40.64);   lonJ = deg2rad( -73.78);
             psiLJ = coorbital.util.greatCircleBearing(latL,lonL,latJ,lonJ);
    assert(abs(rad2deg(psiLJ) - 65.867) < 0.01, ...
        'LAX to JFK bearing wrong: got %.6f deg, expected 65.867 deg', ...
        rad2deg(psiLJ));

%% ASYMMETRY. Measured on this pair:
%%
%%      LAX -> JFK   =  65.8669298403 deg   (east-northeast out of Los Angeles)
%%      JFK -> LAX   = 273.8410316064 deg   (just north of due west out of New York)
%%
%% Turning the outbound course through 180 deg gives 245.867 deg, which is
%% 27.974 deg away from the true return course of 273.841 deg. The two are
%% not reverses of one another and must not be treated as such. A reader who
%% assumes bearing(B,A) = bearing(A,B) + pi, or an implementation that has
%% swapped lat1 with lat2, lands on 245.867 deg and is wrong by that same 28
%% deg -- which over a 3979 km leg is roughly a 1900 km cross-range miss.
%%
%% The physical reason is that the LAX-JFK arc is neither a meridian nor the
%% equator, so it meets the local meridian at a different angle at each end.
%% This is the one property of bearing that distance does not share:
             psiJL = coorbital.util.greatCircleBearing(latJ,lonJ,latL,lonL);
    assert(abs(rad2deg(psiJL) - 273.841) < 0.01, ...
        'JFK to LAX bearing wrong: got %.6f deg, expected 273.841 deg', ...
        rad2deg(psiJL));

%% The signed difference between the return course and the reversed outbound
%% course, folded into (-pi,pi] so that a wrap boundary cannot fake a small
%% number. It must be a LARGE angle -- 27.974 deg -- and the assertion is
%% written as a lower bound so that it fails loudly for any implementation
%% that has quietly made bearing symmetric, whether by returning
%% bearing(A,B) + pi for the return leg or by dropping the dependence on
%% which point is which:
           dRevRad = mod(psiJL - (psiLJ + pi) + pi,2*pi) - pi;
    assert(abs(rad2deg(dRevRad) - 27.974) < 0.01, ...
        ['the LAX/JFK reverse-course offset must be 27.974 deg: got %.6f ', ...
         'deg. A value near zero means bearing has been made symmetric.'], ...
        rad2deg(dRevRad));
    assert(abs(rad2deg(dRevRad)) > 20, ...
        'bearing(B,A) came back as the reverse of bearing(A,B): offset %.6f deg', ...
        rad2deg(dRevRad));

%% The reverse-course OFFSET is not by itself a transposition detector, and
%% it would be a mistake to read it as one. Exchanging lat1 with lat2 inside
%% the formula reflects the geometry in a way that moves BOTH courses --
%% LAX-JFK becomes 86.159 deg and JFK-LAX becomes 294.133 deg -- while
%% leaving their difference at exactly the same 27.974 deg. It is the two
%% absolute-course assertions above that catch that mutation, not this one.
%%
%% So here is a case chosen to catch it on its own terms, with no reference
%% value needed beyond common sense. From the equator at the prime meridian
%% to 60 N one degree of longitude east, the leg is 60 deg of latitude long
%% and one degree of longitude wide: it is very nearly DUE NORTH, 0.577 deg.
%% The reverse leg, from 60 N down to the equator, is very nearly DUE SOUTH,
%% 178.845 deg. A transposed implementation returns exactly those two
%% numbers interchanged, so it reports the northbound leg as heading south:
             psiUp = coorbital.util.greatCircleBearing(0,0,deg2rad(60),deg2rad(1));
             psiDn = coorbital.util.greatCircleBearing(deg2rad(60),0,0,deg2rad(1));
    assert(abs(rad2deg(psiUp) - 0.5773) < 0.001, ...
        ['equator to 60 N one degree east is nearly due north and must be ', ...
         '0.5773 deg: got %.6f deg. A value near 179 deg means lat1 and ', ...
         'lat2 have been exchanged.'],rad2deg(psiUp));
    assert(abs(rad2deg(psiDn) - 178.8453) < 0.001, ...
        ['60 N to the equator one degree east is nearly due south and must ', ...
         'be 178.8453 deg: got %.6f deg.'],rad2deg(psiDn));
    assert(psiUp < deg2rad(90) && psiDn > deg2rad(90), ...
        ['a northbound leg must have a northerly bearing and a southbound ', ...
         'leg a southerly one: got %.6f and %.6f deg'], ...
        rad2deg(psiUp),rad2deg(psiDn));

%% That pair is also a second, independent instance of the asymmetry:
%% reversing 0.5773 deg gives 180.5773 deg, which is 1.732 deg away from the
%% true return course of 178.8453 deg:
           dRevUpD = rad2deg(mod(psiDn - (psiUp + pi) + pi,2*pi) - pi);
    assert(abs(dRevUpD - (-1.732)) < 0.001, ...
        'the 60 N leg reverse-course offset must be -1.732 deg: got %.6f deg', ...
        dRevUpD);

%% A pure meridian leg IS the special case where the reversal does hold, and
%% saying so here keeps the assertion above from reading as a claim that the
%% offset is always large. 10 N to 20 N on one meridian: out is 0, back is pi:
             psiMN = coorbital.util.greatCircleBearing(deg2rad(10),0,deg2rad(20),0);
             psiMS = coorbital.util.greatCircleBearing(deg2rad(20),0,deg2rad(10),0);
    assert(psiMN == 0 && psiMS == pi, ...
        'meridian leg must reverse exactly: got %.17f and %.17f',psiMN,psiMS);

%% A point to itself is degenerate -- there is no arc and no course -- but it
%% must not poison a caller with a NaN that then propagates into a state
%% vector. atan2(0,0) is defined as zero, and that is the answer taken:
              psiS = coorbital.util.greatCircleBearing(0.3,1.1,0.3,1.1);
    assert(~isnan(psiS),'a point to itself returned NaN');
    assert(isfinite(psiS),'a point to itself returned a non-finite bearing');
    assert(psiS >= 0 && psiS < 2*pi, ...
        'a point to itself returned %.17f, outside [0,2*pi)',psiS);

%% Wrapping over a sweep. Every combination below must land in [0,2*pi), the
%% upper end EXCLUDED. Sub-degree westward legs are included deliberately:
%% they produce a bearing a hair below zero before wrapping, and 2*pi minus a
%% hair rounds to exactly 2*pi in double precision, so a bare mod(psi,2*pi)
%% can hand back 2*pi itself and quietly break a half-open contract:
             latSw = deg2rad(-80:20:80);
             lonSw = deg2rad(-170:40:170);
    for ka = 1:numel(latSw)
        for kb = 1:numel(lonSw)
            for kc = 1:numel(latSw)
                for kd = 1:numel(lonSw)
              psiW = coorbital.util.greatCircleBearing(latSw(ka),lonSw(kb), ...
                                                        latSw(kc),lonSw(kd));
                    assert(psiW >= 0 && psiW < 2*pi, ...
                        'bearing out of [0,2*pi): got %.17f for (%g,%g)->(%g,%g) rad', ...
                        psiW,latSw(ka),lonSw(kb),latSw(kc),lonSw(kd));
                end
            end
        end
    end

%% The specific rounding hazard, isolated. A westward leg of 1e-18 rad of
%% longitude from the equator to 1 rad north gives a raw bearing of about
%% -1.4e-18 rad; 2*pi added to that rounds straight back to 2*pi:
            psiEps = coorbital.util.greatCircleBearing(0,0,1,-1e-18);
    assert(psiEps < 2*pi, ...
        'a hair-west leg returned exactly 2*pi (%.17f); the wrap must exclude the upper end', ...
        psiEps);
    assert(psiEps >= 0, ...
        'a hair-west leg returned a negative bearing: %.17f',psiEps);

%% Elementwise, matching the broadcasting contract of coorbital.util.greatCircle:
%% a column in, a column of the same shape out, and scalars broadcast against
%% columns. Each element must equal the scalar call on that element:
              latA = [0; 0.1; -0.4; 1.0];
              lonA = [0; 0.2;  0.9; 3.0];
              latB = [0.5; -0.3; 0.8; -1.1];
              lonB = [1.0;  2.2; -2.0; 0.4];
              psiV = coorbital.util.greatCircleBearing(latA,lonA,latB,lonB);
    assert(isequal(size(psiV),size(latA)), ...
        'output must match the input shape');
    for k = 1:numel(latA)
              psiK = coorbital.util.greatCircleBearing(latA(k),lonA(k), ...
                                                       latB(k),lonB(k));
        assert(psiV(k) == psiK, ...
            'vectorized result disagrees with the scalar call at element %d',k);
    end

%% Scalar broadcast, same rule:
             psiV2 = coorbital.util.greatCircleBearing(latA,lonA,0,0);
    assert(isequal(size(psiV2),size(latA)), ...
        'scalar broadcast must preserve shape');
    for k = 1:numel(latA)
              psiK = coorbital.util.greatCircleBearing(latA(k),lonA(k),0,0);
        assert(psiV2(k) == psiK, ...
            'broadcast result disagrees with the scalar call at element %d',k);
    end

%% And the shapes greatCircle accepts, greatCircleBearing must accept, since
%% the two are called side by side to build a targeting solution:
                dV = coorbital.util.greatCircle(latA,lonA,latB,lonB);
    assert(isequal(size(dV),size(psiV)), ...
        'bearing and distance must broadcast to the same shape');
end
