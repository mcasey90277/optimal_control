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

%% CLOCK TRAP, closed deliberately. The reference is ONE phase spanning both
%  of the driver's, so its solver time is the cumulative clock -- whereas the
%  driver evaluates the schedule on each phase's own clock. Feeding
%  coorbital.guide.prescribed the reference's time would therefore apply the
%  WRONG control in the second half, and the resulting mismatch would read as
%  a continuity failure. Rather than rely on the schedule happening to be
%  constant, that is asserted and the control is then held fixed, which makes
%  the reference independent of any time convention. Anyone extending this to
%  a time-varying schedule must split the reference at the junction instead:
    assert(all(sched.alpha == sched.alpha(1)) && all(sched.sigma == sched.sigma(1)), ...
        'the one-shot reference is only valid for a constant schedule');
            uConst = coorbital.guide.prescribed(0,[],sched);

            odeOne = @(t,x) coorbital.eom.glide3DOF(t,x,uConst,veh,env);
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
%% THREE phases, with the shifted phase in the MIDDLE. The two-phase check
%  above cannot defend the tspan(1) fix on its own: when the shifted phase is
%  LAST, a corrupted time accumulator is computed but never read, because
%  only a SUBSEQUENT phase consumes it. A middle phase is the position that
%  bites. Descend 60 -> 40 -> 25 -> 10 km, with the middle leg given
%  tspan = [10 4010]:
%% ---------------------------------------------------------------------
              ph3a = phase1;
    ph3a.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,40e3);
        ph3a.tspan = [0 4000];
              ph3b = ph3a;
    ph3b.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,25e3);
        ph3b.tspan = [10 4010];
              ph3c = ph3a;
    ph3c.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,10e3);
        ph3c.tspan = [0 4000];
               tr3 = coorbital.prop.phaseRun([ph3a ph3b ph3c],x0,veh,env);

    assert(numel(tr3.junction) == 2,'three phases must give two junctions');
    assert(isequal(unique(tr3.phaseIdx)',[1 2 3]),'samples must be labelled 1, 2, 3');
    assert(all(diff(tr3.t) > 0),'time must strictly increase across both junctions');

%% All three legs really ran, each stopping on its own trigger:
    assert(abs(tr3.junction(1).x(1) - c.rE - 40e3) < 1,'leg 1 must hand over at 40 km');
    assert(abs(tr3.junction(2).x(1) - c.rE - 25e3) < 1,'leg 2 must hand over at 25 km');
    assert(abs(tr3.x(end,1) - c.rE - 10e3) < 1,'leg 3 must terminate at 10 km');

%% Independent reference: integrate each leg standalone and measure its
%  duration off its OWN solver clock. Nothing below touches the driver's
%  accumulator, so the elapsed-time identity is a genuine external check:
            optsAt = @(hs) odeset('RelTol',1e-10,'AbsTol',1e-10, ...
                            'Events',@(t,x) coorbital.prop.eventAltitude(t,x,hs));
              fRef = @(t,x) coorbital.eom.glide3DOF(t,x,uConst,veh,env);
         [tR1,xR1] = ode45(fRef,[0 4000],  x0,          optsAt(40e3));
         [tR2,xR2] = ode45(fRef,[10 4010], xR1(end,:)', optsAt(25e3));
         [tR3,  ~] = ode45(fRef,[0 4000],  xR2(end,:)', optsAt(10e3));
               dur = [tR1(end)-tR1(1); tR2(end)-tR2(1); tR3(end)-tR3(1)];
            cumDur = cumsum(dur);

%% Total elapsed time is the sum of the three legs' actual durations. A
%  middle phase starting at tspan(1) = 10 must contribute its duration, not
%  its duration plus 10:
    assert(abs(tr3.t(end) - cumDur(3)) < 1e-6, ...
        ['elapsed time %.6f s against %.6f s from the three measured leg ' ...
         'durations; tspan(1) leaked into the accumulator'],tr3.t(end),cumDur(3));

%% Each junction sits at the running sum of the preceding legs, and appears
%  on exactly one row:
           medStep = median(diff(tr3.t));
    for kj = 1:2
        assert(abs(tr3.junction(kj).t - cumDur(kj)) < 1e-6, ...
            'junction %d recorded at %.6f s against %.6f s measured', ...
            kj,tr3.junction(kj).t,cumDur(kj));
             rows3 = find(abs(tr3.t - tr3.junction(kj).t) < 1e-9);
        assert(numel(rows3) == 1, ...
            'junction %d appears in %d rows; expected exactly 1',kj,numel(rows3));

%% ...and no anomalous jump in diff(t) at either boundary. Each gap must be
%  an ordinary integration step, not the ~10 s a leaked tspan(1) would open:
            iLast3 = find(tr3.phaseIdx == kj,1,'last');
           iFirst3 = find(tr3.phaseIdx == kj+1,1,'first');
        assert(iFirst3 == iLast3 + 1,'phase labels %d/%d must be contiguous',kj,kj+1);
              gap3 = tr3.t(iFirst3) - tr3.t(iLast3);
        assert(gap3 < 5*medStep, ...
            ['boundary %d opened a %.6f s gap against a %.6f s median step; ' ...
             'that is not an ordinary integration step'],kj,gap3,medStep);
    end

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

%% ---------------------------------------------------------------------
%% A seven-state phase runs and returns a seven-column trajectory. The
%  driver must take its state width from x0, not from the six the glide
%  EOM happens to have; coorbital.eom.massConstant lifts that six-state
%  EOM into the seven-state chain a powered phase needs:
%% ---------------------------------------------------------------------
         schedZero = struct('tGrid',[0 5000],'alpha',[0 0],'sigma',[0 0]);
              eom7 = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
           ph7.eom = eom7;
         ph7.guide = @(t,x) coorbital.guide.prescribed(t,x,schedZero);
     ph7.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
         ph7.tspan = [0 4000];
               x07 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90); veh.mass];
             traj7 = coorbital.prop.phaseRun(ph7,x07,veh,env);

    assert(size(traj7.x,2) == 7,'expected a 7-column trajectory, got %d',size(traj7.x,2));
    assert(all(abs(traj7.x(:,7) - veh.mass) < 1e-12), ...
        'mass must be exactly constant under massConstant');

%% The other outputs must stay consistent with the widened state, and the
%  six flight states must be untouched by the lift. The six-state run of the
%  identical phase is the independent reference -- x(7) is inert in
%  glide3DOF, which reads vehicle mass from veh, so widening the state must
%  not perturb the trajectory:
               nS7 = numel(traj7.t);
    assert(isequal(size(traj7.x),[nS7 7]),'state must be [N x 7]');
    assert(isequal(size(traj7.u),[nS7 2]),'control must still be [N x 2]');
    assert(isequal(size(traj7.phaseIdx),[nS7 1]),'phaseIdx must be [N x 1]');

%% Compared at the TERMINUS, which both runs reach through the same 20 km
%  event, so the comparison needs no shared time grid. Component-scaled so
%  metres and radians are weighed alike. The two runs do in fact land on
%  identical grids today, because ode45 gives a zero-derivative component
%  zero weight in its error norm -- but that is a property of the step
%  controller, not of this driver, so it is not asserted. What IS asserted is
%  that the step COUNT stays the same order, which still catches a mass state
%  that has started driving the integration:
           xScale7 = max(abs(traj.x(end,1:6)'),1);
             dEnd7 = max(abs(traj7.x(end,1:6)' - traj.x(end,1:6)')./xScale7);
    assert(dEnd7 < 1e-9, ...
        'widening the state changed the terminal flight state by %.3e',dEnd7);
    assert(abs(traj7.t(end) - traj.t(end)) < 1e-6, ...
        'widening the state changed the flight time by %.3e s', ...
        abs(traj7.t(end) - traj.t(end)));
    assert(abs(nS7 - numel(traj.t)) < 0.05*numel(traj.t), ...
        ['the seven-state run took %d samples against %d for the identical ' ...
         'six-state run; the inert mass state is driving the integration'], ...
        nS7,numel(traj.t));

%% massConstant must hand the base EOM SIX states, not the whole widened
%  vector. glide3DOF is structurally blind to this -- it unpacks components 1
%  through 6 and ignores anything after -- so the slice is pinned with a
%  probe EOM that REPORTS what it was given instead: its first derivative is
%  the length of the state it saw, and the rest echo that state back.
%
%  Component 7 of the probe input is veh.mass rather than the 7 that the
%  running index would give. massConstant now REFUSES to run when the carried
%  mass x(7) disagrees with the veh.mass its wrapped equations actually divide
%  by, which is the guard that catches a staging jump applied to the state and
%  not to the vehicle. A 7 kg mass against a 900 kg vehicle is exactly the
%  disagreement that guard exists to reject. Nothing below reads component 7 --
%  the slice assertions cover 1 through 6 -- so pinning it to the consistent
%  value changes what this probe proves not at all:
             probe = @(~,x,~,~,~) [numel(x); x(2:end)];
          eomProbe = coorbital.eom.massConstant(probe);
             xdotP = eomProbe(0,[(1:6)'; veh.mass],[0;0],veh,env);

    assert(numel(xdotP) == 7, ...
        ['massConstant returned %d derivatives for a 7-state input; it must ' ...
         'pass 6 states down and append exactly one zero'],numel(xdotP));
    assert(xdotP(1) == 6, ...
        ['the base EOM was handed %d states; massConstant must slice x(1:6) ' ...
         'and pass the base exactly 6'],xdotP(1));
    assert(isequal(xdotP(2:6),(2:6)'), ...
        'the base EOM must receive components 1 through 6 in order, unaltered');
    assert(xdotP(7) == 0, ...
        'the appended mass derivative must be exactly zero, got %.3e',xdotP(7));

%% ---------------------------------------------------------------------
%% The CONTROL width is measured too, not assumed to be 2. A powered phase
%  carries a throttle channel alongside the two aerodynamic angles; the
%  driver must return all three columns, and the extra channel -- which
%  glide3DOF ignores -- must arrive unaltered:
%% ---------------------------------------------------------------------
            thrCmd = 0.73;
            guide3 = @(t,x) [coorbital.guide.prescribed(t,x,schedZero); thrCmd];
              ph3u = ph7;
        ph3u.guide = guide3;
             traj3 = coorbital.prop.phaseRun(ph3u,x07,veh,env);

    assert(size(traj3.u,2) == 3, ...
        ['expected a 3-column control history from a 3-channel guide, got ' ...
         '%d; the control width is still hard-coded'],size(traj3.u,2));
    assert(isequal(size(traj3.u),[numel(traj3.t) 3]),'control must be [N x 3]');
    assert(max(abs(traj3.u(:,3) - thrCmd)) < 1e-12, ...
        'the third control channel must be carried through unaltered');
    assert(max(max(abs(traj3.u(:,1:2)))) < 1e-12, ...
        'the two angle channels must still come from the zero schedule');

%% ---------------------------------------------------------------------
%% STAGE SEPARATION. A phase may carry an optional link, xNext = link(xEnd),
%  applied to its terminal state to form the next phase's initial state.
%  Dropping 500 kg of spent booster must move component 7 by exactly 500 and
%  leave the six flight states untouched ACROSS THE JUNCTION ITSELF. Every
%  expected number below comes from the test's OWN inputs -- the 900 kg
%  initial mass and the 500 kg drop -- never from the driver's bookkeeping:
%% ---------------------------------------------------------------------
%  PER-PHASE VEHICLE, AND WHY THE CHAIN IS BUILT THIS WAY. phaseRun carries
%  ONE vehicle struct for a whole chain, but separation changes the vehicle.
%  The post-separation phase therefore binds its own REBUILT struct inside its
%  equation-of-motion closure -- the same construction BM/run_ballistic uses --
%  so both phases run the REAL coorbital.eom.massConstant, guard and all, and
%  the carried mass agrees with the flown mass on both sides of the jettison.
%
%  DO NOT "SIMPLIFY" THIS BY DROPPING THE REBUILD. glide3DOF divides by
%  veh.mass and never reads x(7), so handing phase 2 the pre-separation 900 kg
%  struct makes it fly at the wrong weight while its own mass history says
%  otherwise. That configuration used to pass this test, on the stated ground
%  that "mass is inert in glide3DOF". Mass is inert in those EQUATIONS and
%  emphatically not inert in the PHYSICS: rebuilding phase 2 at 400 kg moves
%  the terminus by 329 km and the terminal speed by 329 m/s. Confusing the two
%  is exactly what massConstant's guard now refuses to let through -- see the
%  block after this one, which pins that refusal.
             mDrop = 500;
            vehPre = veh;
           vehPost = veh;
      vehPost.mass = veh.mass - mDrop;
            phLnkA = ph7;
        phLnkA.eom = @(t,x,u,vArg,eArg) eom7(t,x,u,vehPre,eArg);
  phLnkA.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,40e3);
       phLnkA.link = @(x) [x(1:6); x(7) - mDrop];
            phLnkB = phLnkA;
        phLnkB.eom = @(t,x,u,vArg,eArg) eom7(t,x,u,vehPost,eArg);
  phLnkB.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
       phLnkB.link = [];
            trLink = coorbital.prop.phaseRun([phLnkA phLnkB],x07,veh,env);

%% The cumulative clock never runs backwards across the junction:
    assert(all(diff(trLink.t) >= 0), ...
        'the cumulative time vector is not monotonically non-decreasing');

    assert(numel(trLink.junction) == 1,'expected one junction between two phases');
    assert(isequal(size(trLink.junction(1).x),[7 1]), ...
        'the junction state must be [7 x 1] in a seven-state chain, got %s', ...
        mat2str(size(trLink.junction(1).x)));

%% Mass is piecewise constant at the two values the test itself prescribes:
             isP1L = trLink.phaseIdx == 1;
             isP2L = trLink.phaseIdx == 2;
    assert(max(abs(trLink.x(isP1L,7) - veh.mass)) < 1e-9, ...
        'phase 1 must fly at the %.1f kg initial mass throughout',veh.mass);
    assert(max(abs(trLink.x(isP2L,7) - (veh.mass - mDrop))) < 1e-9, ...
        'phase 2 must fly at %.1f kg after separation',veh.mass - mDrop);

%% ...and the step between them, located by the phase LABELS rather than by
%  any index convention, is exactly the 500 kg the link removes:
            iLastA = find(isP1L,1,'last');
           iFirstB = find(isP2L,1,'first');
    assert(iFirstB == iLastA + 1,'the phase labels must change exactly once');
             dMass = trLink.x(iLastA,7) - trLink.x(iFirstB,7);
    assert(abs(dMass - mDrop) < 1e-9, ...
        'component 7 dropped %.9f kg across the junction, expected exactly %.1f', ...
        dMass,mDrop);

%% THE ORDERING THE HEADER PROMISES. junction(k).x is the state AS HANDED TO
%  THE NEXT PHASE, i.e. AFTER the link -- the far side of the staging jump.
%  Recording it before the link instead would leave the pre-separation mass
%  here, and a reader debugging a staging event would read the wrong side:
      mJunRecorded = trLink.junction(1).x(7);
    assert(abs(mJunRecorded - (veh.mass - mDrop)) < 1e-9, ...
        ['junction(1).x(7) = %.6f kg, but the junction must record the state ' ...
         'HANDED TO THE NEXT PHASE, i.e. AFTER the %.1f kg separation: ' ...
         '%.1f kg, not the %.1f kg the ending phase carried'], ...
        mJunRecorded,mDrop,veh.mass - mDrop,veh.mass);

%% The link touched component 7 and nothing else. The junction record and the
%  last phase-1 sample sit at the same instant on the cumulative clock, so the
%  six flight states must agree there to 1e-9:
              d16J = max(abs(trLink.junction(1).x(1:6) - trLink.x(iLastA,1:6)'));
    assert(d16J < 1e-9, ...
        'the link perturbed the six flight states by %.3e; it must be the identity there', ...
        d16J);

%% The UNSTAGED reference: the same two legs with no link at all, flying the
%  whole way at 900 kg. It pins the two documented ways of asking for the
%  identity -- an ABSENT link field, and link = [] -- and it is the control
%  against which the staged run's divergence is measured below:
              ph7a = ph7;
    ph7a.terminate = phLnkA.terminate;
              ph7b = ph7;
    ph7b.terminate = phLnkB.terminate;
            trBare = coorbital.prop.phaseRun([ph7a ph7b],x07,veh,env);

    assert(~isfield(ph7a,'link'),'the bare reference must carry no link field at all');
    assert(abs(trBare.junction(1).x(7) - veh.mass) < 1e-12, ...
        'a phase with no link field must hand its mass on unchanged');
    assert(max(abs(trBare.x(:,7) - veh.mass)) < 1e-12, ...
        'an unlinked chain must fly the whole way at the initial mass');

%% UP TO SEPARATION the two runs are the same flight -- phase 1 is 900 kg in
%  both -- so the junction states must coincide exactly. This localises every
%  difference below to the staging event itself, and it catches a rebuilt
%  vehicle leaking backwards into phase 1:
             dJn16 = max(abs(trLink.junction(1).x(1:6) - trBare.junction(1).x(1:6)));
    assert(dJn16 < 1e-9, ...
        ['the staged and unstaged runs disagree at the junction by %.3e; ' ...
         'phase 1 is 900 kg in both and must be the identical flight'],dJn16);

%% AFTER SEPARATION they must DIVERGE, and that is the whole point. An earlier
%  version of this test asserted the opposite -- that the terminus agreed to
%  1e-9 -- which was only ever true because phase 2 was being handed the
%  pre-separation vehicle and so flew at 900 kg regardless of what its mass
%  state said. With phase 2 genuinely rebuilt at 400 kg the aerodynamic
%  accelerations go as 1/m, and the lighter body decelerates harder:
            dEnd16 = max(abs(trLink.x(end,1:6) - trBare.x(end,1:6)));
    assert(dEnd16 > 1, ...
        ['the staged and unstaged runs agree at the terminus to %.3e. Phase 2 ' ...
         'is not actually flying at the rebuilt %.1f kg -- the per-phase ' ...
         'vehicle binding is not reaching the equations of motion.'], ...
        dEnd16,vehPost.mass);
    assert(trLink.x(end,4) < trBare.x(end,4) - 100, ...
        ['the %.1f kg staged body arrived at %.2f m/s against %.2f m/s for the ' ...
         '%.1f kg unstaged one. Drag acceleration goes as 1/m, so the lighter ' ...
         'body must arrive SLOWER, not faster.'], ...
        vehPost.mass,trLink.x(end,4),trBare.x(end,4),veh.mass);

%% ...and the post-separation leg must be flying at 400 kg specifically, not
%  merely at something other than 900. The independent reference is that same
%  leg flown ALONE from the junction state with the rebuilt vehicle: if phase 2
%  of the chain really used vehPost, the two must land on the same terminus.
%  The 900 kg version of the identical leg is run too, so the check is shown to
%  discriminate rather than assumed to:
              solo = ph7b;
          solo.eom = phLnkB.eom;
            trSolo = coorbital.prop.phaseRun(solo,trLink.junction(1).x,vehPost,env);
             xSc16 = max(abs(trLink.x(end,1:6)),1);
              dSolo = max(abs(trSolo.x(end,1:6) - trLink.x(end,1:6))./xSc16);
    assert(dSolo < 1e-9, ...
        ['the chain''s phase 2 and the same leg flown alone at %.1f kg differ ' ...
         'by %.3e; phase 2 did not fly the rebuilt vehicle'],vehPost.mass,dSolo);

             xSolo9 = trLink.junction(1).x;
          xSolo9(7) = veh.mass;
             solo9  = solo;
         solo9.eom  = phLnkA.eom;
            trSolo9 = coorbital.prop.phaseRun(solo9,xSolo9,vehPre,env);
             dSolo9 = max(abs(trSolo9.x(end,1:6) - trLink.x(end,1:6))./xSc16);
    assert(dSolo9 > 1e-3, ...
        ['the 400 kg and 900 kg versions of the post-separation leg agree to ' ...
         '%.3e, so the reference above cannot tell them apart'],dSolo9);

%% link = [] on the LAST phase is the identity and must not disturb the run,
%  which is what lets a struct array mix linked and unlinked phases:
    assert(abs(trLink.x(end,7) - (veh.mass - mDrop)) < 1e-9, ...
        'an empty link must be the identity, leaving %.1f kg at the terminus', ...
        veh.mass - mDrop);

%% Both legs really ran, each stopping on its own trigger:
    assert(abs(trLink.junction(1).x(1) - c.rE - 40e3) < 1, ...
        'the staged leg 1 must hand over at 40 km');
    assert(abs(trLink.x(end,1) - c.rE - 20e3) < 1, ...
        'the staged leg 2 must terminate at 20 km');

%% ---------------------------------------------------------------------
%% THE massConstant MASS GUARD BITES. Take the staged chain above and REMOVE
%  THE PER-PHASE REBUILD -- both phases fall back to the single chain-wide
%  900 kg vehicle phaseRun supplies, while the link still drops the state to
%  400 kg. That is the configuration this whole block exists to outlaw, and
%  it must now REFUSE rather than fly silently at the pre-separation weight:
%% ---------------------------------------------------------------------
         phGuardA = phLnkA;
     phGuardA.eom = eom7;
         phGuardB = phLnkB;
     phGuardB.eom = eom7;
            errG  = [];
    try
        coorbital.prop.phaseRun([phGuardA phGuardB],x07,veh,env);
    catch errG
    end
    assert(~isempty(errG), ...
        ['a chain that jettisons 500 kg without rebuilding the vehicle struct ' ...
         'must raise; massConstant let it fly at the wrong weight instead']);
    assert(strcmp(errG.identifier,'coorbital:massConstant:massMismatch'), ...
        ['expected coorbital:massConstant:massMismatch, got "%s". The guard ' ...
         'must be the thing that fires, not an incidental downstream failure.'], ...
        errG.identifier);
    assert(contains(errG.message,'400.000000000') && ...
           contains(errG.message,'900.000000000'), ...
        ['the mismatch message must name BOTH masses so the reader can see ' ...
         'which side is stale. Message was:\n%s'],errG.message);

%% ...and it must NOT fire on a consistent chain. The same two phases with no
%  jettison run clean, so the guard is discriminating rather than simply
%  hostile to seven-state chains -- trBare above already flew that path:
    assert(max(abs(trBare.x(:,7) - veh.mass)) < 1e-12, ...
        'the consistent guarded chain must have run to completion at %.1f kg',veh.mass);

%% The tolerance must clear an ODE event solve's residual on the mass state,
%  which is a fraction of a nanogram, while rejecting a real staging jump. A
%  perturbation just inside the documented 1e-9 relative budget must pass and
%  one just outside must fail, so the threshold is pinned rather than assumed:
              tolG = 1e-9.*max(1,veh.mass);
             xNear = x07;
          xNear(7) = veh.mass + 0.4.*tolG;
             xFar  = x07;
          xFar(7)  = veh.mass + 2.0.*tolG;
    eom7(0,xNear,[0;0],veh,env);
             threwFar = false;
    try
        eom7(0,xFar,[0;0],veh,env);
    catch
             threwFar = true;
    end
    assert(threwFar, ...
        'a %.3e kg disagreement must exceed the %.3e kg budget and raise', ...
        2.0.*tolG,tolG);

%% THE ABSOLUTE FLOOR in the budget is documented behaviour, so it gets a test
%  that can actually see it. The rule is 1e-9 RELATIVE but never tighter than
%  1e-9 ABSOLUTE, and the floor is the half that nothing above exercises,
%  because every vehicle in this suite masses hundreds of kilograms. It is
%  kept rather than dropped for a concrete reason: a purely relative budget
%  collapses with the mass, and at 0.25 kg it would be 2.5e-10 kg -- TIGHTER
%  than the 3.4e-10 kg residual an ODE burnout-event solve already leaves on a
%  900 kg stack. The guard would then start rejecting ordinary roundoff
%  instead of real mass errors, which is the failure mode the floor exists to
%  prevent. Pinned with a disagreement that sits ABOVE the unfloored budget
%  and BELOW the floored one, so this call succeeds only because of the floor:
            vehSml = veh;
       vehSml.mass = 0.25;
             xSml  = x07;
          xSml(7)  = vehSml.mass + 6e-10;
         threwNear = false;
    try
        eom7(0,xSml,[0;0],vehSml,env);
    catch
         threwNear = true;
    end
    assert(~threwNear, ...
        ['a %.1e kg disagreement on a %.2f kg vehicle was REJECTED. The ' ...
         'budget must never fall below 1e-9 kg absolute; a purely relative ' ...
         'rule would give %.1e kg here, tighter than the %.2e kg residual a ' ...
         'burnout event solve leaves on a 900 kg stack.'], ...
        6e-10,vehSml.mass,1e-9.*vehSml.mass,3.4e-10);

%% ...and the floor is not a blanket amnesty: past it, a light vehicle is
%  still policed:
             xSmlF  = x07;
          xSmlF(7)  = vehSml.mass + 3e-9;
          threwSml  = false;
    try
        eom7(0,xSmlF,[0;0],vehSml,env);
    catch
          threwSml  = true;
    end
    assert(threwSml, ...
        ['a %.1e kg disagreement on a %.2f kg vehicle must still exceed the ' ...
         '1e-9 kg floor and raise'],3e-9,vehSml.mass);

%% A link that changes the state dimension breaks the chain and must say so:
            phBadL = phLnkA;
       phBadL.link = @(x) x(1:6);
              hitL = false;
    try
        coorbital.prop.phaseRun([phBadL phLnkB],x07,veh,env);
    catch errL
              hitL = strcmp(errL.identifier,'coorbital:phaseRun:linkWidth');
    end
    assert(hitL,'a link that drops a state must raise coorbital:phaseRun:linkWidth');

%% ---------------------------------------------------------------------
%% A call with NO phases returns an empty traj rather than throwing. This is
%  degenerate, but it is what the driver did before the widths were measured,
%  and measuring nu must not cost it: there is no phase 1 to ask for a guide
%  evaluation, and asking anyway raised a bare subscript error:
%% ---------------------------------------------------------------------
            trNone = coorbital.prop.phaseRun(struct([]),x07,veh,env);
    assert(isempty(trNone.t),'a zero-phase call must return an empty time vector');
    assert(isempty(trNone.x),'a zero-phase call must return an empty state history');
    assert(isempty(trNone.u),'a zero-phase call must return an empty control history');
    assert(isempty(trNone.phaseIdx),'a zero-phase call must return no phase labels');
    assert(isequal(size(trNone.junction),[0 1]), ...
        'a zero-phase call must give a 0 x 1 junction array, got %s', ...
        mat2str(size(trNone.junction)));
end
