function test_pitchProgram()
%% Purpose:
%
%  Verify coorbital.guide.pitchProgram, the prescribed-pitch-attitude boost
%  guidance law: alpha = theta(t) - gamma, gamma = x(5), with the same
%  linear-interpolation-inside / hold-constant-outside schedule semantics
%  as coorbital.guide.prescribed, plus an optional symmetric alphaMax
%  clamp.
%
%  Every expected value below is a hand-worked literal (arithmetic shown in
%  the comment immediately above the assertion), never a value obtained by
%  re-running pitchProgram's own theta-minus-gamma or clamp formula on its
%  own output.
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
%  sigma = 0 + (-20-0)*(15-0)/(30-0) = -10 deg:
%% =========================================================================
                 x0 = zeros(7,1);
                 u2 = coorbital.guide.pitchProgram(15,x0,schedA);
    assert(abs(u2(1) - deg2rad(45)) < 1e-9, ...
        'midpoint alpha = %.9f rad, expected 45 deg = %.9f rad', ...
        u2(1),deg2rad(45));
    assert(abs(u2(2) - deg2rad(-10)) < 1e-9, ...
        'midpoint sigma = %.9f rad, expected -10 deg = %.9f rad', ...
        u2(2),deg2rad(-10));

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
%  tGrid(end) values (theta = 90 deg, sigma = 20 deg). gamma = 0, so
%  alpha = 90 - 0 = 90 deg:
                 u5 = coorbital.guide.pitchProgram(100,x0,schedA);
    assert(abs(u5(1) - deg2rad(90)) < 1e-9, ...
        'post-grid alpha = %.9f rad, expected 90 deg = %.9f rad', ...
        u5(1),deg2rad(90));
    assert(abs(u5(2) - deg2rad(20)) < 1e-9, ...
        'post-grid sigma = %.9f rad, expected 20 deg = %.9f rad', ...
        u5(2),deg2rad(20));

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
%  gamma = 8 deg -> alpha = 25 - 8 = 17 deg at both t = 0 and t = 500:
%% =========================================================================
            schedSingle = struct('tGrid',0,'theta',deg2rad(25),'sigma',deg2rad(5));
                 x8 = zeros(7,1);
              x8(5) = deg2rad(8);
                u9a = coorbital.guide.pitchProgram(0,x8,schedSingle);
                u9b = coorbital.guide.pitchProgram(500,x8,schedSingle);
    assert(abs(u9a(1) - deg2rad(17)) < 1e-9, ...
        'single-point alpha at t=0 = %.9f rad, expected 17 deg = %.9f rad', ...
        u9a(1),deg2rad(17));
    assert(abs(u9b(1) - deg2rad(17)) < 1e-9, ...
        'single-point alpha at t=500 = %.9f rad, expected 17 deg = %.9f rad', ...
        u9b(1),deg2rad(17));
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
end
