function test_terminalConstraint()
%% Purpose:
%
%  Verify the five shared terminal constraints: that they vanish on a state
%  built to sit exactly on the target manifold, that each residual measures
%  the quantity it claims to with the scaling it claims, that the roll
%  reference cancels from the constraint set, that a rigid rotation about the
%  target plane normal leaves the constraint content untouched -- the freedom
%  VOA's transversality condition exists to price -- and that every degenerate
%  input raises rather than returning a plausible number.
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
%  Michael Casey                                                08/09/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% A target condition, and a state constructed to satisfy it exactly. With
%% uHat perpendicular to wHat, h = r x v = Rd*Vd*cos(gammaD)*wHat identically,
%% so the plane residuals are zero by construction and not by arithmetic
%% accident:
            des.Rd = c.rE + 200e3;
            des.Vd = 7200;
        des.gammaD = deg2rad(2.5);
              incl = deg2rad(38);
          des.wHat = [0; -sin(incl); cos(incl)];
              uHat = [1; 0; 0];
              sHat = cross(des.wHat,uHat);
              sHat = sHat./sqrt(sHat.'*sHat);
               rOn = des.Rd.*uHat;
               vOn = des.Vd.*(sin(des.gammaD).*uHat + cos(des.gammaD).*sHat);

%% The manifold itself:
          [cOn,ach] = coorbital.guide.terminalConstraint(rOn,vOn,des);
    assert(numel(cOn) == 5 && iscolumn(cOn),'c must be a 5 x 1 column');
    assert(max(abs(cOn)) < 1e-14, ...
        'An on-manifold state must give a zero residual; worst was %.3g', ...
        max(abs(cOn)));
    assert(abs(ach.rMag - des.Rd) < 1e-6,'achieved radius should be Rd');
    assert(abs(ach.vMag - des.Vd) < 1e-9,'achieved speed should be Vd');
    assert(abs(ach.gamma - des.gammaD) < 1e-14,'achieved gamma should be gammaD');
    assert(abs(ach.planeAngle) < 1e-14,'achieved plane should be the desired one');
    assert(abs(ach.hMag - des.Rd.*des.Vd.*cos(des.gammaD)) < 1e-3, ...
        'achieved |h| should be Rd*Vd*cos(gammaD)');

%% Each of the first three residuals measures what it says, with the scaling
%% it says. A one per cent radius error is a residual of exactly 0.01 -- a
%% literal, not the implementation's own formula fed back on itself:
           [cR,~] = coorbital.guide.terminalConstraint(1.01.*rOn,vOn,des);
    assert(abs(cR(1) - 0.01) < 1e-12,'c(1) should be r/Rd - 1 = 0.01');
           [cV,~] = coorbital.guide.terminalConstraint(rOn,1.02.*vOn,des);
    assert(abs(cV(2) - 0.02) < 1e-12,'c(2) should be V/Vd - 1 = 0.02');

%% The flight-path-angle residual: hold r and V at their desired values and
%% swing the velocity to a different gamma. With c(1) and c(2) zero the third
%% residual is sin(gamma) - sin(gammaD):
             gamB = deg2rad(6.5);
              vBad = des.Vd.*(sin(gamB).*uHat + cos(gamB).*sHat);
           [cG,~] = coorbital.guide.terminalConstraint(rOn,vBad,des);
    assert(abs(cG(1)) < 1e-14 && abs(cG(2)) < 1e-14, ...
        'the gamma perturbation must not move the radius or speed residual');
    assert(abs(cG(3) - (sin(gamB) - sin(des.gammaD))) < 1e-14, ...
        'c(3) should be sin(gamma) - sin(gammaD)');

%% The plane residuals: rotate the whole on-manifold state rigidly about q1,
%% which tilts the achieved angular momentum out of the target plane by the
%% rotation angle and changes nothing else. hypot(c(4),c(5)) must be sin(tilt):
              bRef = [0; 0; 1];
             q1Raw = cross(bRef,des.wHat);
                q1 = q1Raw./sqrt(q1Raw.'*q1Raw);
    for tilt = deg2rad([0.1 1 5 20])
              rotM = axisRot(q1,tilt);
           [cT,achT] = coorbital.guide.terminalConstraint(rotM*rOn,rotM*vOn,des);
        assert(max(abs(cT(1:3))) < 1e-13, ...
            'a rigid rotation must not move the radius, speed or gamma residual');
        assert(abs(hypot(cT(4),cT(5)) - sin(tilt)) < 1e-12, ...
            'the plane residual magnitude should be sin(tilt) at %.3f deg', ...
            rad2deg(tilt));
        assert(abs(achT.planeAngle - tilt) < 1e-12, ...
            'the reported plane angle should be the tilt');
    end

%% THE ROTATIONAL FREEDOM. Rotating the terminal state rigidly about the
%% target plane normal slides the burnout point around the target orbit, which
%% the five constraints say nothing about. c(1:3) and hypot(c(4),c(5)) must be
%% untouched -- for an OFF-manifold state as well, so this is not just the
%% zero residual being invariant. This is the variation whose transversality
%% condition is VOA's seventh boundary condition:
%% The out-of-plane kick is along wHat, not along q1: a velocity increment
%% parallel to the position vector would leave the angular momentum direction
%% untouched and this whole block would test nothing:
              rOff = 1.03.*rOn;
              vOff = 0.98.*vOn + 40.*des.wHat;
           [cA,~] = coorbital.guide.terminalConstraint(rOff,vOff,des);
    for spin = deg2rad([7 55 123 250])
              rotW = axisRot(des.wHat,spin);
           [cB,~] = coorbital.guide.terminalConstraint(rotW*rOff,rotW*vOff,des);
        assert(max(abs(cB(1:3) - cA(1:3))) < 1e-13, ...
            'a spin about wHat must leave c(1:3) unchanged');
        assert(abs(hypot(cB(4),cB(5)) - hypot(cA(4),cA(5))) < 1e-13, ...
            'a spin about wHat must leave the plane residual magnitude unchanged');
    end

%% THE ROLL REFERENCE CANCELS FROM THE SET. Different bRef choices move c(4)
%% and c(5) individually and leave the surface c = 0 -- and its magnitude --
%% alone. Both facts are asserted, because only checking the invariant would
%% not notice a bRef that was being ignored:
          desRef = des;
       desRef.bRef = [0; 0; 1];
         [cRef1,~] = coorbital.guide.terminalConstraint(rOff,vOff,desRef);
       desRef.bRef = [1; 0; 0];
         [cRef2,~] = coorbital.guide.terminalConstraint(rOff,vOff,desRef);
    assert(abs(hypot(cRef1(4),cRef1(5)) - hypot(cRef2(4),cRef2(5))) < 1e-13, ...
        'bRef must not change the plane residual magnitude');
    assert(max(abs(cRef1(1:3) - cRef2(1:3))) < 1e-14, ...
        'bRef must not touch the first three residuals');
    assert(max(abs(cRef1(4:5) - cRef2(4:5))) > 1e-6, ...
        'these two bRef choices should move c(4) and c(5) individually');
%% And on the manifold, every bRef gives zero:
       desRef.bRef = [0.3; -0.9; 0.2];
         [cRef3,~] = coorbital.guide.terminalConstraint(rOn,vOn,desRef);
    assert(max(abs(cRef3)) < 1e-14,'the manifold is bRef-independent');

%% wHat need not arrive normalised:
          desScal = des;
       desScal.wHat = 1234.5.*des.wHat;
         [cScal,~] = coorbital.guide.terminalConstraint(rOn,vOn,desScal);
    assert(max(abs(cScal)) < 1e-14,'wHat must be normalised internally');

%% Degenerate and malformed inputs raise, with their own identifiers:
    assertRaises(@() coorbital.guide.terminalConstraint([1;2],vOn,des), ...
        'coorbital:terminalConstraint:badVector');
    assertRaises(@() coorbital.guide.terminalConstraint([NaN;0;0],vOn,des), ...
        'coorbital:terminalConstraint:badVector');
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,[Inf;0;0],des), ...
        'coorbital:terminalConstraint:badVector');
             desBad = des;
           desBad.Rd = 0;
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,vOn,desBad), ...
        'coorbital:terminalConstraint:badTarget');
             desBad = des;
           desBad.Vd = -1;
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,vOn,desBad), ...
        'coorbital:terminalConstraint:badTarget');
             desBad = rmfield(des,'gammaD');
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,vOn,desBad), ...
        'coorbital:terminalConstraint:badTarget');
             desBad = des;
       desBad.gammaD = pi/2;
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,vOn,desBad), ...
        'coorbital:terminalConstraint:radialBurnout');
             desBad = des;
         desBad.wHat = [0; 0; 0];
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,vOn,desBad), ...
        'coorbital:terminalConstraint:badPlane');
             desBad = des;
         desBad.bRef = des.wHat;
    assertRaises(@() coorbital.guide.terminalConstraint(rOn,vOn,desBad), ...
        'coorbital:terminalConstraint:badRef');

%% The self-demo runs. evalc keeps its printout out of the suite's output:
    evalc('coorbital.guide.terminalConstraint()');
end

function rotM = axisRot(axisVec,ang)
%% Rodrigues rotation matrix about a unit axis, used only to build test
%% geometries; the library has no rotation utility and this test must not
%% invent an interface for one:
              aHat = axisVec(:)./sqrt(axisVec(:).'*axisVec(:));
             skewM = [      0, -aHat(3),  aHat(2); ...
                      aHat(3),       0, -aHat(1); ...
                     -aHat(2),  aHat(1),       0];
              rotM = eye(3) + sin(ang).*skewM + (1-cos(ang)).*(skewM*skewM);
end

function assertRaises(fn,idWanted)
%% Assert that fn throws with the expected identifier:
              threw = false;
    try
        fn();
    catch err
              threw = true;
        assert(strcmp(err.identifier,idWanted), ...
            'expected %s, got %s',idWanted,err.identifier);
    end
    assert(threw,'expected %s, but nothing was raised',idWanted);
end
