function test_pitchProgram()
%% Purpose:
%
%  Verify coorbital.guide.pitchProgram, the prescribed-pitch-attitude boost
%  guidance law, with the same linear-interpolation-inside /
%  hold-constant-outside schedule semantics as coorbital.guide.prescribed,
%  plus an optional symmetric alphaMax clamp.
%
%  Every expected value below is a hand-worked literal (arithmetic shown in
%  the comment immediately above the assertion), never a value obtained by
%  re-running pitchProgram's own theta-minus-gamma or clamp formula on its
%  own output.
%
%% THE BANKED CASES ARE CHECKED BY SUBSTITUTION, NOT BY REPEATING THE INVERSE:
%
%  At zero bank the law is alpha = theta - gamma and a hand literal is the
%  whole story. At nonzero bank the law inverts the attitude relation
%
%      sin(theta) = sin(gamma)*cos(alpha) + cos(gamma)*cos(sigma)*sin(alpha)
%
%  through an arcsine, and re-running that same inversion in the test would
%  be exactly the self-referential check this file's other assertions avoid:
%  a sign slip inside the arcsine form would cancel out of both sides. So the
%  banked assertions substitute the RETURNED alpha back into the relation
%  above -- the forward direction, which the function never evaluates -- and
%  require the recovered attitude to be the scheduled theta. That is an
%  independent check of the inversion and it is simultaneously the check that
%  the schedule interpolated to the right theta.
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

%% Common schedule: three grid points, theta and sigma both change at every
%  point so interpolation tests cannot pass by one series being flat:
              schedA = struct('tGrid',[0 30 60], ...
                              'theta',deg2rad([30 60 90]), ...
                              'sigma',deg2rad([0 -20 20]));

%% =========================================================================
%% Grid point + state dependence: theta = 30 deg at t = 0, gamma = 10 deg.
%  Hand arithmetic: alpha = 30 - 10 = 20 deg:
%% =========================================================================
                 x1 = zeros(7,1);
              x1(5) = deg2rad(10);
                 u1 = coorbital.guide.pitchProgram(0,x1,schedA);
    assert(abs(u1(1) - deg2rad(20)) < 1e-9, ...
        'grid-point alpha = %.9f rad, expected 20 deg = %.9f rad', ...
        u1(1),deg2rad(20));

%% =========================================================================
%% Linear interpolation, isolated with gamma = 0: t = 15 is halfway between
%  tGrid(1) = 0 (theta = 30 deg) and tGrid(2) = 30 (theta = 60 deg).
%  Hand arithmetic: theta = 30 + (60-30)*(15-0)/(30-0) = 30 + 15 = 45 deg;
%  sigma = 0 + (-20-0)*(15-0)/(30-0) = -10 deg.
%
%  THIS POINT IS BANKED, so alpha is NOT 45 deg. With gamma = 0 the attitude
%  relation reduces to sin(theta) = cos(sigma)*sin(alpha), so
%  alpha = asin(sin(45 deg)/cos(10 deg)) = asin(0.70710678/0.98480775)
%        = asin(0.71801581) = 45.890840951349 deg -- 0.89 deg MORE incidence
%  than the plain difference, because only cos(10 deg) of it lifts the nose.
%  The old alpha = theta - gamma law returned 45 deg here and was wrong:
%% =========================================================================
                 x0 = zeros(7,1);
                 u2 = coorbital.guide.pitchProgram(15,x0,schedA);
    assert(abs(u2(2) - deg2rad(-10)) < 1e-9, ...
        'midpoint sigma = %.9f rad, expected -10 deg = %.9f rad', ...
        u2(2),deg2rad(-10));
    assert(abs(u2(1) - deg2rad(45.890840951349)) < 1e-9, ...
        'midpoint alpha = %.9f rad, expected 45.890841 deg = %.9f rad', ...
        u2(1),deg2rad(45.890840951349));
    assert(u2(1) > deg2rad(45), ...
        ['a banked pitch command must cost MORE incidence than theta - ' ...
         'gamma = 45 deg, not the same: got %.9f deg'],rad2deg(u2(1)));

%% ...and the same point checked by SUBSTITUTION into the attitude relation,
%  which is the direction the function never evaluates. Recovering theta from
%  the returned alpha must give back the interpolated 45 deg:
              th2R = attitudeOf(u2(1),u2(2),x0(5));
    assert(abs(th2R - deg2rad(45)) < 1e-12, ...
        ['substituting the returned alpha back into the attitude relation ' ...
         'recovers %.9f deg, not the scheduled 45 deg'],rad2deg(th2R));

%% A second midpoint, with a nonzero gamma, between tGrid(2) = 30
%  (theta = 60 deg) and tGrid(3) = 60 (theta = 90 deg), at t = 45.
%  Hand arithmetic: theta = 60 + (90-60)*(45-30)/(60-30) = 60 + 15 = 75 deg;
%  sigma = -20 + (20-(-20))*(45-30)/(60-30) = -20 + 20 = 0 deg;
%  alpha = 75 - 15 = 60 deg:
                 x3 = zeros(7,1);
              x3(5) = deg2rad(15);
                 u3 = coorbital.guide.pitchProgram(45,x3,schedA);
    assert(abs(u3(1) - deg2rad(60)) < 1e-9, ...
        'second midpoint alpha = %.9f rad, expected 60 deg = %.9f rad', ...
        u3(1),deg2rad(60));
    assert(abs(u3(2) - deg2rad(0)) < 1e-9, ...
        'second midpoint sigma = %.9f rad, expected 0 deg = %.9f rad', ...
        u3(2),deg2rad(0));

%% =========================================================================
%% Outside the grid, before the first point: t = -10 must clamp to the
%  tGrid(1) values (theta = 30 deg, sigma = 0 deg), same as prescribed.
%  gamma = 0, so alpha = 30 - 0 = 30 deg:
%% =========================================================================
                 u4 = coorbital.guide.pitchProgram(-10,x0,schedA);
    assert(abs(u4(1) - deg2rad(30)) < 1e-9, ...
        'pre-grid alpha = %.9f rad, expected 30 deg = %.9f rad', ...
        u4(1),deg2rad(30));
    assert(abs(u4(2) - deg2rad(0)) < 1e-9, ...
        'pre-grid sigma = %.9f rad, expected 0 deg = %.9f rad', ...
        u4(2),deg2rad(0));

%% Outside the grid, after the last point: t = 100 must clamp to the
%  tGrid(end) values, theta = 90 deg and sigma = 20 deg. That pair is
%  UNREACHABLE and must be refused, which is also what proves the hold
%  happened -- nothing but the held tGrid(end) values produces this refusal.
%  A 90 deg attitude at gamma = 0 needs the thrust axis pointed straight up,
%  i.e. a unit radial component; banked 20 deg the axis can only deliver
%  cos(20 deg) = 0.93969 of one, whatever the incidence. The OLD law returned
%  a nominal-looking 90 deg here for an attitude the vehicle never reaches:
             threw = false;
    try
        coorbital.guide.pitchProgram(100,x0,schedA);
    catch errU
             threw = true;
        assert(strcmp(errU.identifier,'coorbital:pitchProgram:unreachableAttitude'), ...
            'post-grid hold raised %s, expected the unreachable-attitude error', ...
            errU.identifier);
    end
    assert(threw, ...
        ['a 90 deg attitude commanded at 20 deg of bank and gamma = 0 is ' ...
         'unreachable -- cos(20 deg) = 0.93969 < 1 -- and must be refused, ' ...
         'not answered with a plausible 90 deg']);

%% The post-grid hold pinned NUMERICALLY as well, on a schedule whose held
%  terminal attitude is reachable. theta = 50 deg, sigma = 20 deg, gamma = 0:
%  sin(50 deg) = 0.76604444, cos(20 deg) = 0.93969262, so
%  alpha = asin(0.76604444/0.93969262) = asin(0.81520747) = 54.607881 deg:
            schedPost = struct('tGrid',[0 60],'theta',deg2rad([30 50]), ...
                               'sigma',deg2rad([0 20]));
                u5b = coorbital.guide.pitchProgram(100,x0,schedPost);
    assert(abs(u5b(2) - deg2rad(20)) < 1e-9, ...
        'post-grid sigma = %.9f rad, expected 20 deg = %.9f rad', ...
        u5b(2),deg2rad(20));
    assert(abs(rad2deg(u5b(1)) - 54.607881128795) < 1e-9, ...
        'post-grid alpha = %.9f deg, expected 54.607881 deg',rad2deg(u5b(1)));
              th5R = attitudeOf(u5b(1),u5b(2),x0(5));
    assert(abs(th5R - deg2rad(50)) < 1e-12, ...
        ['substituting the returned alpha back into the attitude relation ' ...
         'recovers %.9f deg, not the held 50 deg'],rad2deg(th5R));

%% =========================================================================
%% alphaMax clamp, POSITIVE bound. theta = 90 deg commanded at gamma = 0
%  would give alpha = 90 deg; alphaMax = 20 deg must clamp it to exactly
%  20 deg:
%% =========================================================================
            schedClampP = struct('tGrid',[0 1],'theta',deg2rad([90 90]), ...
                                 'sigma',deg2rad([0 0]),'alphaMax',deg2rad(20));
                 u6 = coorbital.guide.pitchProgram(0,x0,schedClampP);
    assert(abs(u6(1) - deg2rad(20)) < 1e-9, ...
        'positive clamp gave alpha = %.9f rad, expected +20 deg = %.9f rad', ...
        u6(1),deg2rad(20));

%% alphaMax clamp, NEGATIVE bound. theta = -90 deg commanded at gamma = 0
%  would give alpha = -90 deg; alphaMax = 20 deg must clamp it to exactly
%  -20 deg. A clamp built from min() alone (dropping the lower bound) would
%  leave this at -90 deg and fail here even though the positive-bound test
%  above passes:
            schedClampN = struct('tGrid',[0 1],'theta',deg2rad([-90 -90]), ...
                                 'sigma',deg2rad([0 0]),'alphaMax',deg2rad(20));
                 u7 = coorbital.guide.pitchProgram(0,x0,schedClampN);
    assert(abs(u7(1) - deg2rad(-20)) < 1e-9, ...
        'negative clamp gave alpha = %.9f rad, expected -20 deg = %.9f rad', ...
        u7(1),deg2rad(-20));

%% =========================================================================
%% alphaMax clamp with NONZERO gamma, POSITIVE bound. This separates
%  clamping alpha (correct) from clamping theta and then subtracting gamma
%  (a plausible bug): theta = 90 deg, gamma = 10 deg, alphaMax = 20 deg.
%  Clamping alpha: unclamped alpha = 90 - 10 = 80 deg, clamped to 20 deg.
%  Clamping theta instead: clamp(90,20) = 20, then 20 - 10 = 10 deg. The two
%  forms disagree (20 vs 10 deg), so gamma = 0 alone cannot catch this bug:
%% =========================================================================
                 xcp = zeros(7,1);
              xcp(5) = deg2rad(10);
                u6b = coorbital.guide.pitchProgram(0,xcp,schedClampP);
    assert(abs(u6b(1) - deg2rad(20)) < 1e-9, ...
        ['nonzero-gamma positive clamp gave alpha = %.9f rad, expected ' ...
         '20 deg = %.9f rad (clamping theta instead of alpha would give ' ...
         '10 deg)'],u6b(1),deg2rad(20));

%% alphaMax clamp with NONZERO gamma, NEGATIVE bound. theta = -90 deg,
%  gamma = 10 deg, alphaMax = 20 deg. Clamping alpha: unclamped
%  alpha = -90 - 10 = -100 deg, clamped to -20 deg. Clamping theta instead:
%  clamp(-90,20) = -20, then -20 - 10 = -30 deg. Again the two forms
%  disagree (-20 vs -30 deg):
                u7b = coorbital.guide.pitchProgram(0,xcp,schedClampN);
    assert(abs(u7b(1) - deg2rad(-20)) < 1e-9, ...
        ['nonzero-gamma negative clamp gave alpha = %.9f rad, expected ' ...
         '-20 deg = %.9f rad (clamping theta instead of alpha would give ' ...
         '-30 deg)'],u7b(1),deg2rad(-20));

%% =========================================================================
%% alphaMax ABSENT: no clamping, even for a physically absurd 100 deg
%  commanded angle of attack. This must equal 100 deg exactly, unclamped:
%% =========================================================================
            schedNoClamp = struct('tGrid',[0 1],'theta',deg2rad([100 100]), ...
                                  'sigma',deg2rad([0 0]));
                 u8 = coorbital.guide.pitchProgram(0,x0,schedNoClamp);
    assert(abs(u8(1) - deg2rad(100)) < 1e-9, ...
        'unclamped alpha = %.9f rad, expected 100 deg = %.9f rad', ...
        u8(1),deg2rad(100));

%% =========================================================================
%% Single-point schedule: held constant regardless of t (interp1 needs two
%  points, so this exercises the dedicated branch). theta = 25 deg,
%  gamma = 8 deg, sigma = 5 deg. BANKED, so the answer is not the 17 deg the
%  plain difference would give: with R = hypot(sin 8, cos 8 cos 5) and
%  phi = atan2(sin 8, cos 8 cos 5), alpha = asin(sin 25 / R) - phi
%  = 17.069948487809 deg, a shade more incidence than 25 - 8:
%% =========================================================================
            schedSingle = struct('tGrid',0,'theta',deg2rad(25),'sigma',deg2rad(5));
                 x8 = zeros(7,1);
              x8(5) = deg2rad(8);
                u9a = coorbital.guide.pitchProgram(0,x8,schedSingle);
                u9b = coorbital.guide.pitchProgram(500,x8,schedSingle);
    assert(abs(rad2deg(u9a(1)) - 17.069948487809) < 1e-9, ...
        'single-point alpha at t=0 = %.9f deg, expected 17.069948 deg', ...
        rad2deg(u9a(1)));
    assert(abs(rad2deg(u9b(1)) - 17.069948487809) < 1e-9, ...
        'single-point alpha at t=500 = %.9f deg, expected 17.069948 deg', ...
        rad2deg(u9b(1)));
              th9R = attitudeOf(u9a(1),u9a(2),x8(5));
    assert(abs(th9R - deg2rad(25)) < 1e-12, ...
        ['substituting the returned alpha back into the attitude relation ' ...
         'recovers %.9f deg, not the scheduled 25 deg'],rad2deg(th9R));
    assert(abs(u9a(2) - deg2rad(5)) < 1e-9, ...
        'single-point sigma at t=0 = %.9f rad, expected 5 deg = %.9f rad', ...
        u9a(2),deg2rad(5));
    assert(abs(u9b(2) - deg2rad(5)) < 1e-9, ...
        'single-point sigma at t=500 = %.9f rad, expected 5 deg = %.9f rad', ...
        u9b(2),deg2rad(5));

%% =========================================================================
%% STATE DEPENDENCE is the point of this file: hold t fixed at the grid
%  point t = 0 (theta = 30 deg from schedA) and vary gamma across several
%  values. Hand arithmetic for each: alpha = 30 - gamma (deg):
%    gamma =   0 deg -> alpha = 30 -   0 = 30 deg
%    gamma =   5 deg -> alpha = 30 -   5 = 25 deg
%    gamma =  10 deg -> alpha = 30 -  10 = 20 deg
%    gamma = -10 deg -> alpha = 30 - (-10) = 40 deg
%    gamma =  20 deg -> alpha = 30 -  20 = 10 deg
%    gamma = -25 deg -> alpha = 30 - (-25) = 55 deg
%% =========================================================================
             gammasDeg = [0 5 10 -10 20 -25];
      expectedAlphaDeg = [30 25 20 40 10 55];
    for kg = 1:numel(gammasDeg)
                 xg = zeros(7,1);
              xg(5) = deg2rad(gammasDeg(kg));
                 ug = coorbital.guide.pitchProgram(0,xg,schedA);
        assert(abs(ug(1) - deg2rad(expectedAlphaDeg(kg))) < 1e-9, ...
            'gamma = %g deg: alpha = %.9f rad, expected %g deg = %.9f rad', ...
            gammasDeg(kg),ug(1),expectedAlphaDeg(kg),deg2rad(expectedAlphaDeg(kg)));
    end

%% =========================================================================
%% THE ZERO-BANK REDUCTION, BIT FOR BIT. Every trajectory this library ships
%  flies its boost phase unbanked, and their headline numbers are pinned to
%  the last printed digit, so the corrected attitude relation must not merely
%  AGREE with alpha = theta - gamma at sigma = 0, it must return the identical
%  double. It does not do so if routed through the general expression:
%  hypot(sin(gamma),cos(gamma)) is not bit-exactly one and
%  atan2(sin(gamma),cos(gamma)) is not bit-exactly gamma, and measured over
%  200000 random (theta,gamma) pairs in +/-89 deg only 62.9 percent of the
%  results matched exactly, with a worst case of 9.3e-15 rad. That is why
%  pitchProgram writes the zero-bank branch out separately, and this is the
%  assertion that keeps it written out. EQUALITY, deliberately, not a
%  tolerance:
%% =========================================================================
             thGrid = deg2rad(-89:7:89);
             gaGrid = deg2rad(-89:11:89);
    for ka = 1:numel(thGrid)
        for kb = 1:numel(gaGrid)
            schedZ = struct('tGrid',0,'theta',thGrid(ka),'sigma',0);
                xz = zeros(7,1);
             xz(5) = gaGrid(kb);
                uz = coorbital.guide.pitchProgram(0,xz,schedZ);
            assert(uz(1) == thGrid(ka) - gaGrid(kb), ...
                ['the zero-bank reduction is not bit-exact at theta = %.6f ' ...
                 'deg, gamma = %.6f deg: got %.17e against %.17e, a %.3e ' ...
                 'rad difference'],rad2deg(thGrid(ka)),rad2deg(gaGrid(kb)), ...
                uz(1),thGrid(ka) - gaGrid(kb),uz(1) - (thGrid(ka) - gaGrid(kb)));
        end
    end

%% ...and the counterpart: at 45 deg of bank the answer must DIFFER from the
%  old law by a margin no rounding could account for. theta = 30 deg,
%  gamma = 10 deg: the old law gave 20 deg flat; the corrected relation gives
%  30.159383967200 deg, 0.159 deg = 2.78e-3 rad more incidence to buy the
%  same nose-up attitude with the thrust axis rolled 45 deg out of plane.
%  Without this, a mutation reverting the banked branch to theta - gamma
%  would pass the whole file:
            schedB45 = struct('tGrid',0,'theta',deg2rad(30),'sigma',deg2rad(45));
                x45 = zeros(7,1);
             x45(5) = deg2rad(10);
                u45 = coorbital.guide.pitchProgram(0,x45,schedB45);
    assert(abs(rad2deg(u45(1)) - 30.159383967200) < 1e-9, ...
        '45 deg bank alpha = %.12f deg, expected 30.159384 deg',rad2deg(u45(1)));
    assert(abs(u45(1) - deg2rad(20)) > 1e-3, ...
        ['at 45 deg of bank alpha must not still be theta - gamma = 20 deg: ' ...
         'got %.12f deg'],rad2deg(u45(1)));
             th45R = attitudeOf(u45(1),u45(2),x45(5));
    assert(abs(th45R - deg2rad(30)) < 1e-12, ...
        ['substituting the returned alpha back into the attitude relation ' ...
         'recovers %.9f deg, not the scheduled 30 deg'],rad2deg(th45R));
end

function theta = attitudeOf(alpha,sigma,gamma)
%% Purpose:
%
%  Pitch attitude actually achieved by an angle of attack alpha flown at bank
%  sigma on a flight path gamma. This is the FORWARD attitude relation,
%
%      sin(theta) = sin(gamma)*cos(alpha) + cos(gamma)*cos(sigma)*sin(alpha),
%
%  obtained by projecting the thrust/body axis
%  cos(alpha)*vHat + sin(alpha)*(cos(sigma)*eGamma + sin(sigma)*ePsi) onto the
%  local radial. coorbital.guide.pitchProgram solves the same relation the
%  other way, through an arcsine and an atan2 phase, and never evaluates this
%  form -- which is the point: substituting its answer in here is an
%  independent check of the inversion rather than a repetition of it.
%
%% Inputs:
%
%  alpha            [1 x 1]                     Angle of attack (rad)
%
%  sigma            [1 x 1]                     Bank angle (rad)
%
%  gamma            [1 x 1]                     Flight path angle (rad)
%
%% Outputs:
%
%  theta            [1 x 1]                     Achieved pitch attitude (rad),
%                                               principal branch in
%                                               [-pi/2,pi/2]
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------
             theta = asin(sin(gamma).*cos(alpha) ...
                          + cos(gamma).*cos(sigma).*sin(alpha));
end
