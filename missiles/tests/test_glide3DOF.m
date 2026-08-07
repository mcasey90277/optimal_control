function test_glide3DOF()
%% Purpose:
%
%  Two exact checks on the glide equations of motion. In vacuum with no lift
%  the specific orbital energy must be conserved to solver tolerance, and with
%  drag present the mechanical energy must decrease at exactly the rate the
%  drag force does work.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();

%% Vacuum environment: zero density, point-mass gravity, no lift:
        envVac.atmos = @(h) deal(zeros(size(h)),zeros(size(h)), ...
                                 c.T0.*ones(size(h)),ones(size(h)));
         envVac.grav = @coorbital.grav.sphereGrav;
         envVac.aero = @(al,ma,vv) deal(0,0);
       envVac.omegaE = 0;

%% A lofted ballistic arc:
                x0 = [c.rE + 100e3; 0; 0; 5000; deg2rad(20); deg2rad(90)];
              odeF = @(t,x) coorbital.eom.glide3DOF(t,x,[0;0],veh,envVac);
              opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
           [~,X] = ode45(odeF,[0 300],x0,opts);

%% Specific energy V^2/2 - mu/r is conserved in vacuum:
              rArc = X(:,1);
              vArc = X(:,4);
              eArc = 0.5.*vArc.^2 - c.muE./rArc;
            eDrift = abs(eArc - eArc(1))./abs(eArc(1));
    assert(max(eDrift) < 1e-8, ...
        'vacuum energy drifted by %.3e relative; expected < 1e-8',max(eDrift));

%% With drag, energy must fall at exactly the rate drag does work:
          env.atmos = @coorbital.atmos.expAtmos;
           env.grav = @coorbital.grav.sphereGrav;
           env.aero = @coorbital.aero.constLD;
         env.omegaE = 0;

                xD = [c.rE + 40e3; 0; 0; 3000; deg2rad(-2); deg2rad(90)];
             xdotD = coorbital.eom.glide3DOF(0,xD,[0;0],veh,env);

%% Drag power per unit mass, computed independently of the EOM:
        [rhoD,~,~,~] = coorbital.atmos.expAtmos(xD(1) - c.rE);
          [~,CDd] = coorbital.aero.constLD(0,0,veh);
             dragA = 0.5*rhoD*xD(4)^2*veh.Sref*CDd/veh.mass;
          [grD,~] = coorbital.grav.sphereGrav(xD(1),xD(3));

%% dV/dt must equal -dragAccel - g*sin(gamma) exactly:
            vdotEx = -dragA - grD*sin(xD(5));
    assert(abs(xdotD(4) - vdotEx) < 1e-12, ...
        'dV/dt = %.9f does not match -D/m - g sin(gamma) = %.9f',xdotD(4),vdotEx);

%% Kinematics: with gamma = 0 the radius must be stationary:
                xF = [c.rE + 40e3; 0; 0; 3000; 0; deg2rad(90)];
             xdotF = coorbital.eom.glide3DOF(0,xF,[0;0],veh,env);
    assert(abs(xdotF(1)) < 1e-12,'dr/dt must vanish at gamma = 0');

%% Heading due east at the equator moves longitude, not latitude:
    assert(xdotF(2) > 0,'eastward flight must increase longitude');
    assert(abs(xdotF(3)) < 1e-12,'eastward flight at the equator must not change latitude');

%% Rotating-frame Jacobi integral: the rotating analogue of vacuum energy
%  conservation. With gravity spherical, no lift, and rotation switched on,
%  V^2/2 - muE/r - (1/2)*omegaE^2*r^2*cos(lat)^2 is the Jacobi integral for a
%  particle in a uniformly rotating frame and must be conserved exactly. An
%  off-equator start with a heading that is neither due north nor due east
%  exercises the centrifugal terms nonzero. NOTE: dJ/dt depends only on Vdot,
%  rdot, and latdot -- never on gammadot or psidot -- so this conservation
%  law is structurally blind to the true (2*om, workless) Coriolis terms,
%  which live only in gammadot and psidot. The instantaneous spot check
%  below closes that gap by recomputing those terms directly:
             envRot = envVac;
         envRot.omegaE = c.omegaE;

                x0R = [c.rE + 100e3; 0; 0.3; 5000; deg2rad(20); 1.0];
              odeFR = @(t,x) coorbital.eom.glide3DOF(t,x,[0;0],veh,envRot);
             [~,XR] = ode45(odeFR,[0 300],x0R,opts);

                rR = XR(:,1);
              latR = XR(:,3);
                vR = XR(:,4);
            jacobi = 0.5.*vR.^2 - c.muE./rR - 0.5.*c.omegaE.^2.*rR.^2.*cos(latR).^2;
            jDrift = abs(jacobi - jacobi(1))./abs(jacobi(1));
    assert(max(jDrift) < 1e-8, ...
        'rotating Jacobi integral drifted by %.3e relative; expected < 1e-8',max(jDrift));

%% Rotating-branch instantaneous spot check: independently recompute the
%  full gammadot and psidot rotating terms (both the 2*om Coriolis pieces
%  and the om^2 centrifugal pieces) at one off-equator, non-cardinal-heading
%  state, using the model functions directly rather than the EOM's own
%  internal expressions. This is what actually certifies the Coriolis terms,
%  since the Jacobi integral above cannot see them:
                omS = c.omegaE;
              latR2 = 0.3;
              psiR2 = 1.0;
            gammaR2 = deg2rad(20);
                VR2 = 5000;
                rR2 = c.rE + 100e3;
             xdotR2 = coorbital.eom.glide3DOF(0,[rR2;0;latR2;VR2;gammaR2;psiR2],[0;0],veh,envRot);

            [gr2,~] = coorbital.grav.sphereGrav(rR2,latR2);
    gammadotNonrotEx = -(gr2/VR2 - VR2/rR2)*cos(gammaR2);
      psidotNonrotEx = VR2*cos(gammaR2)*sin(psiR2)*tan(latR2)/rR2;

       gammadotRotEx = 2*omS*cos(latR2)*sin(psiR2) ...
                        + omS^2*rR2*cos(latR2)* ...
                          (cos(gammaR2)*cos(latR2) + sin(gammaR2)*cos(psiR2)*sin(latR2))/VR2;
         psidotRotEx = 2*omS*(sin(latR2) - cos(latR2)*cos(psiR2)*tan(gammaR2)) ...
                        + omS^2*rR2*sin(latR2)*cos(latR2)*sin(psiR2)/(VR2*cos(gammaR2));

    assert(abs(xdotR2(5) - (gammadotNonrotEx + gammadotRotEx)) < 1e-12, ...
        'rotating gammadot = %.9f does not match independent calc = %.9f', ...
        xdotR2(5),gammadotNonrotEx + gammadotRotEx);
    assert(abs(xdotR2(6) - (psidotNonrotEx + psidotRotEx)) < 1e-12, ...
        'rotating psidot = %.9f does not match independent calc = %.9f', ...
        xdotR2(6),psidotNonrotEx + psidotRotEx);

%% Guard tests: each singular state must fail loudly with the documented
%  error identifier, not silently return NaN or fall through with no error:
              xPole = [c.rE + 40e3; 0; pi/2; 3000; deg2rad(-2); deg2rad(90)];
         threwPole = false;
    try
        coorbital.eom.glide3DOF(0,xPole,[0;0],veh,env);
    catch errPole
        threwPole = true;
        assert(strcmp(errPole.identifier,'coorbital:glide3DOF:polarSingularity'), ...
            'wrong error identifier at the pole: %s',errPole.identifier);
    end
    assert(threwPole,'glide3DOF must throw at the pole (lat = pi/2)');

              xVert = [c.rE + 40e3; 0; 0; 3000; pi/2; deg2rad(90)];
         threwVert = false;
    try
        coorbital.eom.glide3DOF(0,xVert,[0;0],veh,env);
    catch errVert
        threwVert = true;
        assert(strcmp(errVert.identifier,'coorbital:glide3DOF:verticalFlight'), ...
            'wrong error identifier in vertical flight: %s',errVert.identifier);
    end
    assert(threwVert,'glide3DOF must throw at gamma = pi/2');

              xSlow = [c.rE + 40e3; 0; 0; 0.5; deg2rad(-2); deg2rad(90)];
         threwSlow = false;
    try
        coorbital.eom.glide3DOF(0,xSlow,[0;0],veh,env);
    catch errSlow
        threwSlow = true;
        assert(strcmp(errSlow.identifier,'coorbital:glide3DOF:zeroSpeed'), ...
            'wrong error identifier at low speed: %s',errSlow.identifier);
    end
    assert(threwSlow,'glide3DOF must throw below V = 1 m/s');

%% Lift and off-equator spot check: independent recomputation of gammadot,
%  psidot, and the cos(lat) placement in londot, at a state with everything
%  turned on -- nonzero latitude, angle of attack, and bank angle:
               latS = 0.4;
               psiS = 1.0;
             gammaS = -0.05;
                 VS = 5000;
                 hS = 50e3;
                 rS = c.rE + hS;
             alphaS = 0.1;
             sigmaS = 0.4;
                 xS = [rS; 0; latS; VS; gammaS; psiS];
              xdotS = coorbital.eom.glide3DOF(0,xS,[alphaS;sigmaS],veh,env);

%% Independent recomputation using the model functions directly, not the
%  EOM's own internal expressions:
    [rhoS,~,~,aSndS] = coorbital.atmos.expAtmos(hS);
           [CLs,~] = coorbital.aero.constLD(alphaS,VS./aSndS,veh);
           [grS,~] = coorbital.grav.sphereGrav(rS,latS);
              qbarS = 0.5*rhoS*VS^2;
                aLS = qbarS*veh.Sref*CLs/veh.mass;

         gammadotEx = aLS*cos(sigmaS)/VS - (grS/VS - VS/rS)*cos(gammaS);
           psidotEx = aLS*sin(sigmaS)/(VS*cos(gammaS)) ...
                      + VS*cos(gammaS)*sin(psiS)*tan(latS)/rS;

    assert(abs(xdotS(5) - gammadotEx) < 1e-10, ...
        'gammadot = %.9f does not match independent calc = %.9f',xdotS(5),gammadotEx);
    assert(abs(xdotS(6) - psidotEx) < 1e-10, ...
        'psidot = %.9f does not match independent calc = %.9f',xdotS(6),psidotEx);

%% Pin the cos(lat) placement in londot: londot*r*cos(lat) must equal
%  V*cos(gamma)*sin(psi) to machine tolerance:
          londotLHS = xdotS(2)*rS*cos(latS);
          londotRHS = VS*cos(gammaS)*sin(psiS);
    assert(abs(londotLHS - londotRHS)/abs(londotRHS) < 1e-10, ...
        'londot*r*cos(lat) = %.9f does not match V cos(gamma) sin(psi) = %.9f', ...
        londotLHS,londotRHS);
end
