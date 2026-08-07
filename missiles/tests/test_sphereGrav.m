function test_sphereGrav()
%% Purpose:
%
%  Verify point-mass gravity: correct surface value, inverse-square falloff,
%  no latitude dependence, and a zero latitudinal component.
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

%% Surface gravity is mu/r^2 = 9.7983 m/s^2, worked out by hand from the
%% WGS-84 constants (NOT g0 = 9.80665; the ~0.09% difference is real, from
%% oblateness and rotation) so a formula error cannot hide behind the
%% implementation's own arithmetic:
       [gr,gLat] = coorbital.grav.sphereGrav(c.rE,0);
    assert(abs(gr - 9.7983) < 1e-4,'surface gravity wrong: got %.6f m/s^2',gr);
    assert(gr > 0,'gr must be positive downward');
    assert(gLat == 0,'spherical gravity has no latitudinal component');

%% Inverse square: doubling the radius quarters the acceleration, which a
%% wrong exponent (e.g. r^3) would fail:
      [gr2,~] = coorbital.grav.sphereGrav(2*c.rE,0);
    assert(abs(gr2 - gr/4) < 1e-12,'not inverse-square');

%% No latitude dependence for a point mass:
      [grEq,~] = coorbital.grav.sphereGrav(c.rE,0);
      [grPo,~] = coorbital.grav.sphereGrav(c.rE,pi/2);
    assert(abs(grEq - grPo) < 1e-12,'point-mass gravity must not vary with latitude');

%% Vectorizes:
              rVec = c.rE + (0:10000:100000)';
      [grV,gLatV] = coorbital.grav.sphereGrav(rVec,zeros(size(rVec)));
    assert(isequal(size(grV),size(rVec)),'gr not [N x 1]');
    assert(isequal(size(gLatV),size(rVec)),'gLat not [N x 1]');
    assert(all(diff(grV) < 0),'gravity must decrease with radius');
end
