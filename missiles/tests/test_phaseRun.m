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

              iJun = rows(1);
    assert(isequal(size(traj2.junction(1).x),[6 1]), ...
        'the junction state must be stored as [6 x 1]');

%% Continuity, checked against a genuinely INDEPENDENT reference. Comparing
%  the junction record against traj2.x would be circular -- both descend from
%  the same ode45 output -- so the same trajectory is integrated ONCE over the
%  whole span with a DIFFERENT solver at a tighter tolerance, and the junction
%  record is measured against that. Note the reference time is taken from the
%  phase LABELS in the output, never from the junction record itself, so a
%  junction recorded at the wrong instant cannot hide:
            iLast1 = find(traj2.phaseIdx == 1,1,'last');
           iFirst2 = find(traj2.phaseIdx == 2,1,'first');
    assert(iFirst2 == iLast1 + 1, ...
        'the phase labels must change exactly once, at the boundary');

            odeOne = @(t,x) coorbital.eom.glide3DOF(t,x, ...
                            coorbital.guide.prescribed(t,x,sched),veh,env);
           optsOne = odeset('RelTol',1e-12,'AbsTol',1e-12, ...
                            'Events',@(t,x) coorbital.prop.eventAltitude(t,x,5e3));
            solOne = ode89(odeOne,[0 8000],x0,optsOne);

%% Scale each component by its own magnitude so metres and radians are
%  compared on equal terms:
            xScale = max(abs(traj2.junction(1).x),1);
             dLast = max(abs(deval(solOne,traj2.t(iLast1)) ...
                             - traj2.junction(1).x)./xScale);
            dFirst = max(abs(deval(solOne,traj2.t(iFirst2)) ...
                             - traj2.junction(1).x)./xScale);
    assert(min(dLast,dFirst) < 1e-8, ...
        ['the junction state is not the one-shot trajectory state at the ' ...
         'phase boundary: %.3e before, %.3e after'],dLast,dFirst);

%% The recorded junction TIME must also be consistent with that same
%  independent solution:
             dAtJn = max(abs(deval(solOne,tJun) ...
                             - traj2.junction(1).x)./xScale);
    assert(dAtJn < 1e-8, ...
        'the junction state disagrees with the one-shot solution by %.3e',dAtJn);

%% The one-shot run must also reach the same terminal state as the split run,
%  or the split has changed the physics rather than merely partitioned it:
            dFinal = max(abs(deval(solOne,solOne.x(end)) - traj2.x(end,:)') ...
                         ./max(abs(traj2.x(end,:)'),1));
    assert(dFinal < 1e-8, ...
        'splitting the run into two phases changed the terminal state by %.3e',dFinal);

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
%% A phase whose tspan does not start at zero must contribute only its own
%  DURATION to the cumulative clock. Shifting phase 2 from [0 4000] to
%  [10 4010] leaves the span length and the dynamics untouched, so the whole
%  cumulative trajectory must be unchanged -- no 10 s gap opened at the
%  junction, no 10 s added to the elapsed time:
%% ---------------------------------------------------------------------
          phase2Sh = phase2;
    phase2Sh.tspan = [10 4010];
             trShf = coorbital.prop.phaseRun([phase1 phase2Sh],x0,veh,env);

    assert(abs(trShf.t(end) - traj2.t(end)) < 1e-6, ...
        ['shifting phase 2 tspan to [10 4010] changed the elapsed time by ' ...
         '%.6f s; the clock must follow the phase duration, not tspan(1)'], ...
        abs(trShf.t(end) - traj2.t(end)));
    assert(abs(trShf.junction(1).t - traj2.junction(1).t) < 1e-9, ...
        'the junction time must not depend on the next phase tspan(1)');
    assert(all(diff(trShf.t) > 0),'time must strictly increase');

%% No gap at the junction: the step across the boundary must be an ordinary
%  integration step, identical to the unshifted run's step, not ~10 s:
           iLastSh = find(trShf.phaseIdx == 1,1,'last');
          iFirstSh = find(trShf.phaseIdx == 2,1,'first');
            gapShf = trShf.t(iFirstSh) - trShf.t(iLastSh);
            gapRef = traj2.t(iFirst2)  - traj2.t(iLast1);
    assert(abs(gapShf - gapRef) < 1e-6, ...
        ['a %.6f s step opened at the junction against %.6f s in the ' ...
         'unshifted run; tspan(1) leaked into the cumulative clock'], ...
        gapShf,gapRef);
    assert(gapShf < 0.5*max(diff(traj2.t)) + 1, ...
        'the junction step %.6f s is not an ordinary integration step',gapShf);

%% The shifted run must still record its junction on exactly one row:
          rowsShf = find(abs(trShf.t - trShf.junction(1).t) < 1e-9);
    assert(numel(rowsShf) == 1, ...
        'the shifted junction instant appears in %d rows; expected 1',numel(rowsShf));

%% ---------------------------------------------------------------------
%% Solver tolerances are overridable through env, defaulting to 1e-10. A
%  loose override must visibly change the step selection while still landing
%  on the event:
%% ---------------------------------------------------------------------
           envLoose = env;
 envLoose.odeRelTol = 1e-6;
 envLoose.odeAbsTol = 1e-6;
            trLoose = coorbital.prop.phaseRun(phase1,x0,veh,envLoose);

    assert(numel(trLoose.t) < 0.5*numel(traj.t), ...
        ['a 1e-6 tolerance produced %d samples against %d at 1e-10; the ' ...
         'override did not reach ode45'],numel(trLoose.t),numel(traj.t));
    assert(abs(trLoose.x(end,1) - c.rE - 20e3) < 1, ...
        'the loose run must still terminate on the 20 km event');
    assert(abs(trLoose.t(end) - traj.t(end)) < 0.1, ...
        'the loose run diverged from the tight run by %.4f s', ...
        abs(trLoose.t(end) - traj.t(end)));

%% Absent those fields the default must still be 1e-10, i.e. the untouched
%  env must reproduce the tight run exactly:
    assert(numel(traj.t) > 2*numel(trLoose.t), ...
        'the default tolerance must remain the tight 1e-10');

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
