function test_phaseRun()
%% Purpose:
%
%  Verify the phase driver terminates on an altitude event, returns a
%  monotonic time vector, records the junction state between phases, and
%  hands the terminal state of one phase to the next without a discontinuity.
%  Also exercises the two pieces the driver leans on: the one-sided altitude
%  event, which must ignore an ascending crossing, and the prescribed control
%  schedule, which must interpolate inside its grid and clamp outside it.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;

%% ---------------------------------------------------------------------
%% eventAltitude: the documented three-output ODE event interface. The
%  direction must be -1 so only descending crossings are located; the
%  behavioural consequence is tested on a lofted arc further down:
%% ---------------------------------------------------------------------
             xProbe = [c.rE + 32e3; 0; 0; 4000; 0; deg2rad(90)];
    [vEv,isTermEv,dirEv] = coorbital.prop.eventAltitude(0,xProbe,20e3);
    assert(abs(vEv - 12e3) < 1e-9, ...
        'event value = %.6f, expected altitude above hStop = 12000 m',vEv);
    assert(isTermEv == 1,'the altitude event must be terminal');
    assert(dirEv == -1,'direction must be -1 so ascending crossings are ignored');

%% The event value must vanish exactly at the stop altitude:
    [vZero,~,~] = coorbital.prop.eventAltitude(0,[c.rE + 20e3;0;0;4000;0;0],20e3);
    assert(abs(vZero) < 1e-9,'event value must vanish at hStop, got %.6e',vZero);

%% ---------------------------------------------------------------------
%% prescribed: linear inside the grid, clamped outside, [2 x 1] always:
%% ---------------------------------------------------------------------
             schTri = struct('tGrid',[0 10 20],'alpha',[0 0.2 0.1], ...
                             'sigma',[0 -0.4 0.4]);
              uMid = coorbital.guide.prescribed(5,[],schTri);
    assert(isequal(size(uMid),[2 1]),'prescribed must return [2 x 1]');
    assert(abs(uMid(1) - 0.1) < 1e-12,'alpha at t = 5 must be 0.1, got %.9f',uMid(1));
    assert(abs(uMid(2) + 0.2) < 1e-12,'sigma at t = 5 must be -0.2, got %.9f',uMid(2));

              uLo = coorbital.guide.prescribed(-500,[],schTri);
    assert(max(abs(uLo - [0;0])) < 1e-12,'schedule must clamp to its first point');
              uHi = coorbital.guide.prescribed(500,[],schTri);
    assert(max(abs(uHi - [0.1;0.4])) < 1e-12,'schedule must clamp to its last point');

%% A one-point schedule is a constant, not an interp1 error:
             schPt = struct('tGrid',0,'alpha',0.3,'sigma',-0.2);
              uPt = coorbital.guide.prescribed(7,[],schPt);
    assert(max(abs(uPt - [0.3;-0.2])) < 1e-12,'a one-point schedule must be constant');

%% ---------------------------------------------------------------------
%% A level-bank glide that ends when it reaches 20 km:
%% ---------------------------------------------------------------------
             sched = struct('tGrid',[0 5000],'alpha',[0 0],'sigma',[0 0]);
        phase1.eom = @coorbital.eom.glide3DOF;
      phase1.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
  phase1.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
      phase1.tspan = [0 4000];

                x0 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
              traj = coorbital.prop.phaseRun(phase1,x0,veh,env);

%% Time is monotonic and the run stopped on the event, not the horizon:
    assert(all(diff(traj.t) > 0),'time must increase monotonically');
    assert(traj.t(end) < 4000,'integration hit the horizon instead of the event');

%% It stopped at 20 km:
              hEnd = traj.x(end,1) - c.rE;
    assert(abs(hEnd - 20e3) < 1,'terminated at %.1f m, expected 20000 m',hEnd);

%% Shapes are as documented, and every field agrees with numel(t):
                nS = numel(traj.t);
    assert(isequal(size(traj.t),[nS 1]),'t must be [N x 1]');
    assert(isequal(size(traj.x),[nS 6]),'state must be [N x 6]');
    assert(isequal(size(traj.u),[nS 2]),'control must be [N x 2]');
    assert(isequal(size(traj.phaseIdx),[nS 1]),'phaseIdx must be [N x 1]');
    assert(all(traj.phaseIdx == 1),'a single-phase run must be labelled phase 1');

%% A single-phase run has no junctions, and the empty must still be a
%  0 x 1 struct array carrying the documented fields:
    assert(isstruct(traj.junction),'junction must be a struct array');
    assert(isequal(size(traj.junction),[0 1]), ...
        'a single-phase run must give a 0 x 1 junction array, got %s', ...
        mat2str(size(traj.junction)));
    assert(all(isfield(traj.junction,{'t','x'})), ...
        'the empty junction array must still declare fields t and x');

%% The initial condition is preserved exactly at the first sample:
    assert(max(abs(traj.x(1,:)' - x0)) == 0,'the first sample must be x0 exactly');
    assert(traj.t(1) == 0,'a run starting at tspan(1) = 0 must begin at t = 0');

%% ---------------------------------------------------------------------
%% Two phases: the second continues from the first without a jump:
%% ---------------------------------------------------------------------
            phase2 = phase1;
  phase2.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,5e3);
      phase2.tspan = [0 4000];
             traj2 = coorbital.prop.phaseRun([phase1 phase2],x0,veh,env);

    assert(numel(traj2.junction) == 1,'expected one junction between two phases');
    assert(isequal(size(traj2.junction),[1 1]), ...
        'junction must be a [(P-1) x 1] array, got %s',mat2str(size(traj2.junction)));

%% Time never runs backwards, and no instant is recorded twice -- the shared
%  sample at the handoff must appear exactly once:
    assert(all(diff(traj2.t) > 0), ...
        'time must strictly increase; a zero step means a duplicated sample');

%% Locate the junction by its recorded TIME rather than by any index
%  convention. It must match exactly one row of the trajectory:
              tJun = traj2.junction(1).t;
              rows = find(abs(traj2.t - tJun) < 1e-9);
    assert(numel(rows) == 1, ...
        'the junction instant appears in %d rows; expected exactly 1',numel(rows));

%% ...and the state recorded there must be the state on that row:
              iJun = rows(1);
              jump = abs(traj2.x(iJun,:) - traj2.junction(1).x(:)');
    assert(max(jump) < 1e-6,'state jumped by %.3e across the junction',max(jump));
    assert(isequal(size(traj2.junction(1).x),[6 1]), ...
        'the junction state must be stored as [6 x 1]');

%% Shapes stay consistent for the multi-phase case:
               nS2 = numel(traj2.t);
    assert(isequal(size(traj2.x),[nS2 6]),'state must be [N x 6]');
    assert(isequal(size(traj2.u),[nS2 2]),'control must be [N x 2]');
    assert(isequal(size(traj2.phaseIdx),[nS2 1]),'phaseIdx must be [N x 1]');
    assert(isequal(unique(traj2.phaseIdx)',[1 2]), ...
        'a two-phase run must label its samples 1 and 2 only');
    assert(all(diff(traj2.phaseIdx) >= 0),'phase labels must not go backwards');

%% Both phases really ran: the first ended at 20 km, the second at 5 km:
    assert(abs(traj2.x(iJun,1) - c.rE - 20e3) < 1, ...
        'phase 1 must hand over at 20 km, handed over at %.1f m', ...
        traj2.x(iJun,1) - c.rE);
    assert(abs(traj2.x(end,1) - c.rE - 5e3) < 1, ...
        'phase 2 must terminate at 5 km, terminated at %.1f m', ...
        traj2.x(end,1) - c.rE);

%% Independent reference: re-integrate the second phase on its own from the
%  RECORDED junction state and confirm the driver's tail matches, both in
%  elapsed time and in terminal state. This is what pins the time offset --
%  the driver's last sample must sit tJun + tRef(end) on the cumulative
%  clock, not at tRef(end):
              odeRef = @(t,x) coorbital.eom.glide3DOF(t,x, ...
                              coorbital.guide.prescribed(t,x,sched),veh,env);
             optsRef = odeset('RelTol',1e-10,'AbsTol',1e-10, ...
                              'Events',@(t,x) coorbital.prop.eventAltitude(t,x,5e3));
       [tRef,xRef] = ode45(odeRef,[0 4000],traj2.junction(1).x,optsRef);

              dtDrv = traj2.t(end) - tJun;
    assert(abs(dtDrv - tRef(end)) < 1e-6, ...
        'phase 2 spanned %.9f s on the cumulative clock but %.9f s standalone', ...
        dtDrv,tRef(end));
    assert(max(abs(traj2.x(end,:) - xRef(end,:))) < 1e-6, ...
        'the driver terminal state differs from the standalone re-integration');

%% The cumulative clock must also start where phase 1 started:
    assert(traj2.t(1) == 0,'the cumulative clock must start at 0');
    assert(traj2.t(end) > tJun,'the second phase must advance the clock');

%% ---------------------------------------------------------------------
%% The event is one-sided. A lofted arc launched BELOW the trigger altitude
%  climbs through it, coasts to apogee far above it, and only then descends
%  back through it. A two-sided or ascending event would cut the run off on
%  the way up, leaving a peak altitude barely above the trigger and a
%  positive terminal flight path angle:
%% ---------------------------------------------------------------------
        phLoft.eom = @coorbital.eom.glide3DOF;
      phLoft.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
  phLoft.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,50e3);
      phLoft.tspan = [0 3000];
             xLoft = [c.rE + 30e3; 0; 0; 4000; deg2rad(15); deg2rad(90)];
           trLoft = coorbital.prop.phaseRun(phLoft,xLoft,veh,env);

              hMax = max(trLoft.x(:,1)) - c.rE;
    assert(hMax > 60e3, ...
        'peak altitude %.1f m: the event fired on the ascending crossing',hMax);
    assert(trLoft.x(end,5) < 0, ...
        'terminal flight path angle %.4f rad: the run must end descending', ...
        trLoft.x(end,5));
    assert(abs(trLoft.x(end,1) - c.rE - 50e3) < 1, ...
        'the lofted arc must stop at 50 km, stopped at %.1f m', ...
        trLoft.x(end,1) - c.rE);
    assert(trLoft.t(end) < 3000,'the lofted arc hit the horizon, not the event');

%% ---------------------------------------------------------------------
%% The control history is the schedule evaluated on the returned grid at
%  time SINCE PHASE START, which is what coorbital.guide.prescribed
%  documents. With a ramp in both channels, evaluating on the cumulative
%  clock instead would be visible in phase 2:
%% ---------------------------------------------------------------------
            schRamp = struct('tGrid',[0 600],'alpha',[0 0.2],'sigma',[0 0.1]);
        phRamp.eom = @coorbital.eom.glide3DOF;
      phRamp.guide = @(t,x) coorbital.guide.prescribed(t,x,schRamp);
  phRamp.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,40e3);
      phRamp.tspan = [0 2000];
            phRamp2 = phRamp;
 phRamp2.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,25e3);
            trRamp = coorbital.prop.phaseRun([phRamp phRamp2],x0,veh,env);

%% Phase-local clock, rebuilt from the recorded junction time:
              tLoc = trRamp.t;
             isPh2 = trRamp.phaseIdx == 2;
       tLoc(isPh2) = trRamp.t(isPh2) - trRamp.junction(1).t;
             tClip = min(max(tLoc,0),600);
             aWant = 0.2.*tClip./600;
             sWant = 0.1.*tClip./600;

    assert(max(abs(trRamp.u(:,1) - aWant)) < 1e-12, ...
        'alpha history off by %.3e from the schedule at phase-local time', ...
        max(abs(trRamp.u(:,1) - aWant)));
    assert(max(abs(trRamp.u(:,2) - sWant)) < 1e-12, ...
        'sigma history off by %.3e from the schedule at phase-local time', ...
        max(abs(trRamp.u(:,2) - sWant)));

%% A nonzero bank must actually reach the dynamics, or the control plumbing
%  is decorative: the banked run must have turned out of the equatorial plane:
    assert(max(abs(trRamp.x(:,3))) > deg2rad(0.5), ...
        'the banked run never left the equator; sigma did not reach the EOM');
end
