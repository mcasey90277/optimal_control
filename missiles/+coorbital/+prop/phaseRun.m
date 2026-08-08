function traj = phaseRun(phases,x0,veh,env)
%% Purpose:
%
%  Integrate a sequence of trajectory phases, handing the terminal state of
%  each phase to the next and recording the junction states so a downstream
%  optimizer can enforce them as linkage constraints.
%
%% Inputs:
%
%  phases           [1 x P] struct              Each with fields:
%                                               eom       function handle
%                                               guide     function handle
%                                               terminate function handle
%                                               tspan     [1 x 2] (s)
%                                               link      function handle,
%                                                         OPTIONAL, see Notes
%                                               veh       struct, OPTIONAL,
%                                                         the vehicle THIS
%                                                         phase flies; absent
%                                                         or [] falls back to
%                                                         the chain-wide veh
%                                                         argument, see Notes
%
%  x0               [nx x 1]                    Initial state. Its length nx
%                                               sets the state dimension for
%                                               the whole call
%
%  veh              Struct                      Vehicle parameters, used by
%                                               every phase that does not
%                                               carry its own veh field
%
%  env              Struct                      Environment model handles:
%                                               atmos, grav, aero, omegaE, and
%                                               prop for any POWERED phase.
%                                               prop is read only by
%                                               coorbital.eom.boost3DOF; an
%                                               unpowered chain may omit it and
%                                               glide3DOF ignores it if present.
%                                               Optional odeRelTol and
%                                               odeAbsTol override the ode45
%                                               tolerances, default 1e-10 each.
%                                               odeAbsTol may be a VECTOR of
%                                               nx per-component tolerances,
%                                               see Notes
%
%% Outputs:
%
%  traj             Struct                      Trajectory:
%                                               t        [N x 1] (s), cumulative
%                                               x        [N x nx] state
%                                               u        [N x nu] control (rad
%                                                        for the angle channels)
%                                               phaseIdx [N x 1] 1-based phase
%                                               junction [(P-1) x 1] struct
%                                                        with fields t and x,
%                                                        x being the state
%                                                        AFTER the link, see
%                                                        Notes
%                                               phaseEnd [P x 1] string, per
%                                                        phase, "event" or
%                                                        "tspan", see Notes
%                                               endedOnEvent
%                                                        [P x 1] logical, the
%                                                        same fact as a mask
%
%% Notes:
%
%  TWO CLOCKS, AND WHICH ONE EACH PARTY SEES. There is exactly one convention
%  and it is worth stating in one place, because "time since phase start" is
%  ambiguous English and has been read both ways.
%
%      OUTPUT clock, traj.t   -- CUMULATIVE and zero-based. Each phase is
%                                referenced to its own tspan(1) and offset by
%                                the elapsed time of the phases before it, so
%                                a phase given tspan = [10 50] contributes
%                                40 s and opens no gap at the junction.
%
%      MODEL clock, the t handed to eom, guide and terminate
%                             -- the PHASE'S OWN tspan VALUE, exactly as the
%                                integrator produced it. It is NOT shifted to
%                                zero. A phase given tspan = [10 50] evaluates
%                                its guide at 10 through 50, not 0 through 40.
%
%  So "time since phase start" in the input block of every eom, guide, event
%  and propulsion model in this library means the phase's tspan clock, which
%  coincides with elapsed time within the phase only for the usual
%  tspan = [0 T]. Every schedule that ships is written against that clock and
%  every chain in this library starts each phase at tspan(1) = 0, where the
%  two readings agree. The behaviour is deliberately NOT changed here: the
%  schedules are written against it, and re-basing the model clock would move
%  every time-dependent command in the library.
%
%  A GUIDE MUST BE A PURE FUNCTION OF (t,x), and the reason is arithmetic, not
%  taste. Per phase it is evaluated:
%
%      - ONCE at (tspan(1),x0), for phase 1 only, to measure nu;
%      - as many times as ode45's stages require, which is not knowable in
%        advance and is not the same as the number of output samples;
%      - then ONCE MORE PER RETURNED SAMPLE, numel(tk) calls, on the
%        reconstruction pass that fills traj.u -- all of them AFTER every
%        solver-stage call for that phase has already been made.
%
%  A guide that latches, counts, integrates or otherwise remembers therefore
%  reports a traj.u that is not the control the integrator actually flew, and
%  worse, the reconstruction pass leaves its hidden state advanced before the
%  next phase begins. Carry controller state in the integrated state vector
%  and read it back out of traj.x; do not hide it in the closure.
%
%  The sample at a phase boundary is recorded ONCE and carries the OUTGOING
%  phase's control. A control discontinuity across a boundary therefore does
%  not appear in traj.u -- only the value just before the handoff is stored.
%  A consumer needing both sides must evaluate the incoming phase's guide at
%  the junction time itself.
%
%  Every phase in one call shares ONE state dimension, nx = numel(x0), because
%  the phases are concatenated into a single traj.x. A six-state unpowered
%  glide therefore cannot be chained directly onto a seven-state powered
%  phase: wrap the six-state EOM with coorbital.eom.massConstant, which
%  appends dm/dt = 0, and run the whole chain seven-state. The control width
%  nu is likewise measured, not assumed: phase 1's guide is called ONCE
%  before the phase loop, at (tspan(1), x0), purely to size traj.u. See the
%  purity requirement above for the full evaluation count.
%
%  A phase may carry an optional link handle, xNext = link(xEnd), applied to
%  that phase's terminal state to form the initial state of the NEXT phase --
%  stage separation dropping spent booster mass, for instance. Absent, or set
%  to [], it is the identity; [] is how a struct array mixes linked and
%  unlinked phases, MATLAB requiring every element of a struct array to carry
%  the same fields. The link must preserve nx.
%
%  The link is applied BEFORE the junction is recorded, so junction(k).x is
%  the state AS HANDED TO THE NEXT PHASE: the value on the FAR side of a
%  staging jump, not the terminal state of the phase that just ended. The
%  near side is the last traj.x row labelled phase k, which sits at the same
%  instant on the cumulative clock. A link on the FINAL phase is still
%  evaluated and validated -- it can raise linkWidth, and any side effect it
%  has still happens -- but it does not affect the output, there being no
%  next phase to seed and no junction to record.
%
%  A phase may carry an optional veh struct, ph.veh, which is the vehicle
%  THAT phase flies. Absent or [], the chain-wide veh argument is used and the
%  behaviour is identical to a call that never heard of the field -- so this
%  is purely additive. It exists because a staged chain does not fly one
%  vehicle: after separation the airframe is a different object, mass AND
%  Sref, CL, LD, and a phase link can only change the STATE, never the
%  parameters. Before this field existed every chain entry script worked
%  around it by binding a per-phase vehicle inside its own EOM closure and
%  ignoring the forwarded argument; those closures still work unchanged, and
%  coorbital.eom.massConstant still guards the mass half of the divergence.
%
%  A phase records HOW IT ENDED, because ode45 reaching tspan(end) and ode45
%  stopping on a terminal event are different outcomes and the state history
%  alone does not distinguish them. traj.phaseEnd(kp) is "event" when a
%  located event coincides with the phase's last returned time and "tspan"
%  otherwise; traj.endedOnEvent is the same fact as a logical mask. Neither
%  is an error: a time-limited phase is a legitimate design, and the entry
%  scripts do their own termination diagnosis. The point is that they no
%  longer have to re-derive it from the trajectory.
%
%  odeAbsTol may be a VECTOR of nx per-component tolerances, which odeset
%  accepts directly; the default is the SCALAR 1e-10 for every component.
%  A scalar absolute tolerance across a state whose components are a radius in
%  metres (6.4e6), angles in radians (order 1) and a mass in kilograms (1e3 to
%  1e4) is not a physical error budget -- 1e-10 is 0.1 nm on the radius and
%  0.1 mg on the mass, and only the tightest of the seven does any work. The
%  scalar default is kept because every pinned number in the test suite was
%  measured under it, and a per-component vector is offered instead, e.g.
%
%      env.odeAbsTol = [1e-6; 1e-12; 1e-12; 1e-9; 1e-12; 1e-12; 1e-9];
%
%  which reads as a micrometre of radius, a picoradian of angle, a nanometre
%  per second of speed and a microgram of mass. Supply nx entries: six for an
%  unpowered chain, seven with the mass state.
%
%  A call with no phases returns an empty traj and does not integrate.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Michael Casey  Any state dimension; optional phase link      08/07/2026
%  Michael Casey  Optional per-phase veh; record phase end      08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
             sched = struct('tGrid',[0 5000],'alpha',[0 0],'sigma',[0 0]);
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,0);
          ph.tspan = [0 4000];
                x0 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);
    figure('color',[1 1 1]);
    plot(traj.t,(traj.x(:,1)-c.rE)./1000,'linewidth',1.5); grid on;
    xlabel('time (s)'); ylabel('altitude (km)');
    return;
end

%% One cell per phase, concatenated once at the end. Segment lengths are not
%% known until ode45 returns, so the cells are the preallocation:
               nPh = numel(phases);
              tSeg = cell(nPh,1);
              xSeg = cell(nPh,1);
              uSeg = cell(nPh,1);
            idxSeg = cell(nPh,1);
          junction = repmat(struct('t',[],'x',[]),max(nPh-1,0),1);
          phaseEnd = strings(nPh,1);
      endedOnEvent = false(nPh,1);
              tOff = 0;
             xCurr = x0(:);

%% State and control widths are MEASURED, never assumed, so a seven-state
%% powered chain, or a control carrying a throttle channel, needs no edit
%% here. nu costs one guide evaluation, which a zero-phase call must not
%% make -- there is no phase 1 to ask, and such a call returns empty:
                nx = numel(xCurr);
                nu = 0;
if nPh > 0
                nu = numel(phases(1).guide(phases(1).tspan(1),xCurr));
end

%% Solver tolerances, overridable through env so a later task can trade
%% accuracy for speed without editing this file:
            relTol = 1e-10;
            absTol = 1e-10;
if isfield(env,'odeRelTol')
            relTol = env.odeRelTol;
end
if isfield(env,'odeAbsTol')
            absTol = env.odeAbsTol;
end

%% Integrate each phase in turn:
for kp = 1:nPh
                ph = phases(kp);

%% The vehicle THIS phase flies. Optional and purely additive: absent or [],
%% the chain-wide argument is used and nothing changes. A staged chain needs
%% it because a phase link can move the STATE across separation but can never
%% replace the parameters -- mass, Sref, CL, LD -- that the equations divide by:
          phaseVeh = veh;
    if isfield(ph,'veh') && ~isempty(ph.veh)
          phaseVeh = ph.veh;
    end
              odeF = @(t,x) ph.eom(t,x,ph.guide(t,x),phaseVeh,env);
              opts = odeset('RelTol',relTol,'AbsTol',absTol, ...
                            'Events',ph.terminate);
   [tk,xk,te,~,ie] = ode45(odeF,ph.tspan,xCurr,opts);

%% How the phase ended. Reaching tspan(end) and stopping on a terminal event
%% are different outcomes that the state history alone cannot tell apart, and
%% a caller that has to re-derive it from the trajectory will re-derive it
%% differently in each script. Recorded, not enforced: a time-limited phase is
%% a legitimate design. The test is that a LOCATED event coincides with the
%% last returned sample, which a non-terminal event elsewhere in the arc does
%% not satisfy:
  endedOnEvent(kp) = ~isempty(ie) && ~isempty(te) && te(end) == tk(end);
    if endedOnEvent(kp)
      phaseEnd(kp) = "event";
    else
      phaseEnd(kp) = "tspan";
    end

%% Rebuild the control history on the returned grid. A phase disagreeing with
%% the measured control width could not be concatenated into one traj.u, so it
%% is caught here rather than deep inside the assignment:
                uk = zeros(numel(tk),nu);
    for kt = 1:numel(tk)
               ukt = ph.guide(tk(kt),xk(kt,:)');
        if numel(ukt) ~= nu
            error('coorbital:phaseRun:controlWidth', ...
                'Phase %d returned %d control(s) where phase 1 returned %d.', ...
                kp,numel(ukt),nu);
        end
          uk(kt,:) = ukt(:).';
    end

%% Drop the duplicated first sample when continuing from a previous phase:
                k0 = 1;
    if kp > 1
                k0 = 2;
    end

%% Referencing every phase to its own tspan(1) keeps the cumulative clock right
%% for a phase that does not start at zero:
          tSeg{kp} = tOff + (tk(k0:end) - tk(1));
          xSeg{kp} = xk(k0:end,:);
          uSeg{kp} = uk(k0:end,:);
        idxSeg{kp} = kp*ones(numel(tk)-k0+1,1);

%% Carry the state forward, applying the optional link that maps this phase's
%% terminal state to the next phase's initial state. The link runs BEFORE the
%% junction is recorded, so junction(kp).x is the state as HANDED TO THE NEXT
%% PHASE -- the far side of a staging jump, not the near side:
             xCurr = xk(end,:)';
              tOff = tOff + (tk(end) - tk(1));
    if isfield(ph,'link') && ~isempty(ph.link)
             xCurr = ph.link(xCurr);
             xCurr = xCurr(:);
        if numel(xCurr) ~= nx
            error('coorbital:phaseRun:linkWidth', ...
                'The phase %d link returned %d state(s); the chain runs %d.', ...
                kp,numel(xCurr),nx);
        end
    end

%% Record the junction:
    if kp < nPh
    junction(kp).t = tOff;
    junction(kp).x = xCurr;
    end
end

%% Assemble:
            traj.t = vertcat(tSeg{:});
            traj.x = vertcat(xSeg{:});
            traj.u = vertcat(uSeg{:});
     traj.phaseIdx = vertcat(idxSeg{:});
     traj.junction = junction;
     traj.phaseEnd = phaseEnd;
 traj.endedOnEvent = endedOnEvent;
end
