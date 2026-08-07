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
end
