function traj = run_glide()
%% Purpose:
%
%  Propagate a prescribed-control hypersonic glide from entry interface to a
%  terminal altitude and report range, flight time, and peak conditions.
%  Everything a routine run needs to change lives in the USER PARAMETERS
%  block; nothing below it should require editing.
%
%  Angles and distances in the USER PARAMETERS block are in degrees and
%  kilometres because that is how a user thinks about them. They are
%  converted to the library's SI units (m, m/s, rad, s) immediately after the
%  block, in one place, and nothing below that point works in any other unit.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun:
%                                               t        [N x 1] (s)
%                                               x        [N x 6] state
%                                               u        [N x 2] control (rad)
%                                               phaseIdx [N x 1]
%                                               junction [0 x 1] struct
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Resolve paths so the script runs from anywhere:
              here = fileparts(mfilename('fullpath'));
    addpath(here);
    addpath(fullfile(here,'..'));

%% ========================= USER PARAMETERS ==============================
%% Entry interface state -- where and how the glide begins:
            hEntry = 60;               %km,  entry altitude above the sphere [30 .. 120]
            vEntry = 6000;             %m/s, planet-relative speed at entry [1000 .. 7500]
        gammaEntry = -1.0;             %deg, flight path angle; NEGATIVE is descending [-10 .. 0]
          psiEntry = 90;               %deg, heading clockwise from north [0 .. 360]; 90 = due east
          latEntry = 0;                %deg, geocentric latitude [-89 .. 89]; the poles are singular
          lonEntry = 0;                %deg, longitude [-180 .. 180]

%% Termination -- when to stop integrating:
             hStop = 5;                %km, terminal altitude, DESCENDING crossing only [0 .. hEntry)
           tMaxSim = 4000;             %s,  integration horizon; raise it if the glide is cut short

%% Prescribed control schedule -- held constant for the whole flight here:
        alphaAngle = 0;                %deg, angle of attack [0 .. 30]; ignored by the constLD aero model
         bankAngle = 0;                %deg, bank [-90 .. 90]; 0 = all lift up, +/-90 = all lift into the turn

%% Vehicle -- point this at any file in this folder returning a veh struct:
         vehicleFn = @vehicle_hgv;     %handle, see vehicle_hgv.m and coorbital.util.vehicleDefaults

%% Model selection -- swap a handle here to raise fidelity; each replacement
%% must keep the same signature as the one it replaces:
           atmosFn = @coorbital.atmos.expAtmos;   %[rho,P,T,a] = atmosFn(h)
            gravFn = @coorbital.grav.sphereGrav;  %[gr,gLat]   = gravFn(r,lat)
            aeroFn = @coorbital.aero.constLD;     %[CL,CD]     = aeroFn(alpha,mach,veh)
         earthSpin = false;            %true enables the Coriolis and centrifugal terms

%% Output:
         showPlots = true;             %false to skip the figures, e.g. under matlab -batch
%% ======================= END USER PARAMETERS ============================

%% Convert the user block to library SI units. This is the ONLY unit
%% conversion in the file; everything past this point is m, m/s, rad and s:
           hEntryM = hEntry.*1000;
            hStopM = hStop.*1000;
      gammaEntryRd = deg2rad(gammaEntry);
        psiEntryRd = deg2rad(psiEntry);
        latEntryRd = deg2rad(latEntry);
        lonEntryRd = deg2rad(lonEntry);
          alphaRad = deg2rad(alphaAngle);
           bankRad = deg2rad(bankAngle);

%% Sanity-check the user block before spending time in the integrator:
    assert(hStopM < hEntryM, ...
        'hStop (%.1f km) must be below hEntry (%.1f km) or the run ends at t = 0.', ...
        hStop,hEntry);
    assert(tMaxSim > 0,'tMaxSim must be positive.');
    assert(vEntry > 1,'vEntry must exceed 1 m/s; the equations are singular below that.');

                 c = coorbital.util.missileConst();
               veh = vehicleFn();

%% Assemble the environment from the handles chosen above:
         env.atmos = atmosFn;
          env.grav = gravFn;
          env.aero = aeroFn;
        env.omegaE = 0;
    if earthSpin
        env.omegaE = c.omegaE;
    end

%% Build the control schedule and the single glide phase:
             sched = struct('tGrid',[0 tMaxSim], ...
                            'alpha',[alphaRad alphaRad], ...
                            'sigma',[bankRad  bankRad]);
           ph.eom = @coorbital.eom.glide3DOF;
         ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
     ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,hStopM);
         ph.tspan = [0 tMaxSim];

                x0 = [c.rE + hEntryM; ...
                      lonEntryRd; ...
                      latEntryRd; ...
                      vEntry; ...
                      gammaEntryRd; ...
                      psiEntryRd];

%% Propagate:
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);

%% Derived quantities for the summary:
                nS = numel(traj.t);
               hKm = (traj.x(:,1) - c.rE)./1000;
                 V = traj.x(:,4);
    [rho,~,~,aSnd] = atmosFn(traj.x(:,1) - c.rE);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;

%% Aerodynamic coefficients sample by sample, so a table-driven aero model
%% reports correctly here and not just the constant one:
               CLv = zeros(nS,1);
               CDv = zeros(nS,1);
    for k = 1:nS
        [CLv(k),CDv(k)] = aeroFn(traj.u(k,1),mach(k),veh);
    end
             aLift = qbar.*veh.Sref.*CLv./veh.mass;
             aDrag = qbar.*veh.Sref.*CDv./veh.mass;

%% Sensed aerodynamic load, the quantity a structures engineer calls the
%% deceleration; gravity is not sensed and is excluded:
             nLoad = sqrt(aLift.^2 + aDrag.^2)./c.g0;

%% Great-circle ground range from the entry point to the last sample:
            angRad = coorbital.util.greatCircle(latEntryRd,lonEntryRd, ...
                                                traj.x(end,3),traj.x(end,2));
           rangeKm = c.rE.*angRad./1000;

%% Peaks, with the times at which they occur:
      [qMax,kQ]    = max(qbar);
      [machMax,kM] = max(mach);
      [nMax,kN]    = max(nLoad);
      [hMax,kH]    = max(hKm);

%% Label for the rotation switch in the model line below:
           spinTxt = 'OFF';
    if earthSpin
           spinTxt = 'ON';
    end

%% Diagnose WHY the propagation stopped, so a truncated run can never be read
%% as a completed glide:
             hEndM = traj.x(end,1) - c.rE;
              tEnd = traj.t(end);
       stopNominal = false;
    if abs(hEndM - hStopM) < 1
           stopWhy = sprintf('reached the %.1f km terminal altitude (nominal)',hStop);
       stopNominal = true;
    elseif tEnd >= tMaxSim - 1e-6
           stopWhy = sprintf(['hit the %.0f s integration horizon at %.2f km, ' ...
                              'still flying'],tMaxSim,hEndM./1000);
    else
           stopWhy = sprintf(['stopped early at t = %.2f s, h = %.2f km ' ...
                              '(integrator failure or an unmodelled event)'],tEnd,hEndM./1000);
    end

%% Report:
    fprintf('\n');
    fprintf('===== Glide summary =====================================\n');
    fprintf('  termination        %s\n',stopWhy);
    if ~stopNominal
        fprintf('  *** CAUTION ***    the trajectory below is TRUNCATED, not a completed glide\n');
    end
    fprintf('\n');
    fprintf('  Entry conditions\n');
    fprintf('    altitude         %10.2f km\n',hEntry);
    fprintf('    speed            %10.2f %-5s  (Mach %.2f)\n',vEntry,'m/s',mach(1));
    fprintf('    flight path      %10.2f %-5s  (negative is descending)\n',gammaEntry,'deg');
    fprintf('    heading          %10.2f %-5s  (clockwise from north)\n',psiEntry,'deg');
    fprintf('    latitude         %10.4f deg\n',latEntry);
    fprintf('    longitude        %10.4f deg\n',lonEntry);
    fprintf('\n');
    fprintf('  Overall\n');
    fprintf('    flight time      %10.2f s\n',tEnd);
    fprintf('    ground range     %10.2f %-5s  (great circle on the r = %.1f km sphere)\n', ...
            rangeKm,'km',c.rE./1000);
    fprintf('    central angle    %10.4f deg\n',rad2deg(angRad));
    fprintf('    peak altitude    %10.2f %-5s  (at t = %.1f s)\n',hMax,'km',traj.t(kH));
    fprintf('    samples          %10d       (ode45 adaptive steps)\n',nS);
    fprintf('\n');
    fprintf('  Terminal conditions\n');
    fprintf('    altitude         %10.3f km\n',hEndM./1000);
    fprintf('    speed            %10.2f %-5s  (Mach %.2f)\n',V(end),'m/s',mach(end));
    fprintf('    flight path      %10.2f deg\n',rad2deg(traj.x(end,5)));
    fprintf('    heading          %10.2f deg\n',rad2deg(traj.x(end,6)));
    fprintf('    latitude         %10.4f deg\n',rad2deg(traj.x(end,3)));
    fprintf('    longitude        %10.4f deg\n',rad2deg(traj.x(end,2)));
    fprintf('\n');
    fprintf('  Peak loads\n');
    fprintf('    dynamic pressure %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qMax./1000,'kPa',traj.t(kQ),hKm(kQ));
    fprintf('    Mach number      %10.2f       (at t = %.1f s)\n',machMax,traj.t(kM));
    fprintf('    deceleration     %10.2f %-5s  (%.2f m/s^2 sensed aero load, at t = %.1f s)\n', ...
            nMax,'g',nMax.*c.g0,traj.t(kN));
    fprintf('\n');
    fprintf('  Models: atmos %s | grav %s | aero %s | Earth rotation %s\n', ...
            func2str(atmosFn),func2str(gravFn),func2str(aeroFn),spinTxt);
    fprintf('  Vehicle: %s, mass %.1f kg, Sref %.3f m^2, CL %.3f, L/D %.2f (PLACEHOLDER values)\n', ...
            func2str(vehicleFn),veh.mass,veh.Sref,veh.CL,veh.LD);
    fprintf('=========================================================\n\n');

%% Plots:
if showPlots
    figure('color',[1 1 1],'name','Glide time histories');
    subplot(3,1,1); plot(traj.t,hKm,'linewidth',1.5); grid on;
    ylabel('altitude (km)'); title('Prescribed-control glide');
    subplot(3,1,2); plot(traj.t,V,'linewidth',1.5); grid on;
    ylabel('speed (m/s)');
    subplot(3,1,3); plot(traj.t,qbar./1000,'linewidth',1.5); grid on;
    ylabel('q (kPa)'); xlabel('time (s)');

    figure('color',[1 1 1],'name','Ground track');
    plot(rad2deg(traj.x(:,2)),rad2deg(traj.x(:,3)),'linewidth',1.5); grid on; hold on;
    plot(lonEntry,latEntry,'o','markersize',8,'linewidth',1.5);
    plot(rad2deg(traj.x(end,2)),rad2deg(traj.x(end,3)),'s','markersize',8,'linewidth',1.5);
    xlabel('longitude (deg)'); ylabel('latitude (deg)');
    legend('ground track','entry','terminus','location','best');
    title(sprintf('Ground track, %.0f km great-circle range',rangeKm));

    figure('color',[1 1 1],'name','Mach and load factor');
    subplot(2,1,1); plot(traj.t,mach,'linewidth',1.5); grid on;
    ylabel('Mach (-)'); title('Speed regime and sensed aerodynamic load');
    subplot(2,1,2); plot(traj.t,nLoad,'linewidth',1.5); grid on;
    ylabel('load factor (g)'); xlabel('time (s)');
end

%% Hand the trajectory back only when the caller asked for it. Typing
%% "run_glide" at the prompt should leave the summary on screen, not bury it
%% under a dump of the whole struct as ans:
if nargout == 0
    clear traj;
end
end
