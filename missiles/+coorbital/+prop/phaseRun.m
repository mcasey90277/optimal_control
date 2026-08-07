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
%
%  x0               [nx x 1]                    Initial state. Its length nx
%                                               sets the state dimension for
%                                               the whole call
%
%  veh              Struct                      Vehicle parameters
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
%                                               tolerances, default 1e-10 each
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
%
%% Notes:
%
%  Each phase is integrated on its own tspan and then referenced to that
%  tspan's start, so a phase given tspan = [10 50] contributes 40 s to the
%  cumulative clock and opens no gap at the junction. The guide, by contrast,
%  still sees the phase's own tspan-relative time, so a schedule written
%  against tspan = [10 50] is evaluated at 10 through 50.
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
%  before the phase loop, at (tspan(1), x0), purely to size traj.u. A
%  STATEFUL guide -- one that latches, counts calls, or otherwise remembers
%  what it has been asked -- will therefore see one evaluation more than the
%  integrator makes, and must tolerate it.
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
%  A call with no phases returns an empty traj and does not integrate.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Michael Casey  Any state dimension; optional phase link      08/07/2026
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
              odeF = @(t,x) ph.eom(t,x,ph.guide(t,x),veh,env);
              opts = odeset('RelTol',relTol,'AbsTol',absTol, ...
                            'Events',ph.terminate);
          [tk,xk] = ode45(odeF,ph.tspan,xCurr,opts);

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
end
