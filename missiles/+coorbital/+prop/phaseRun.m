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
%
%  x0               [6 x 1]                     Initial state
%
%  veh              Struct                      Vehicle parameters
%
%  env              Struct                      Environment model handles:
%                                               atmos, grav, aero, omegaE.
%                                               Optional odeRelTol and
%                                               odeAbsTol override the ode45
%                                               tolerances, default 1e-10 each
%
%% Outputs:
%
%  traj             Struct                      Trajectory:
%                                               t        [N x 1] (s), cumulative
%                                               x        [N x 6] state
%                                               u        [N x 2] control (rad)
%                                               phaseIdx [N x 1] 1-based phase
%                                               junction [(P-1) x 1] struct
%                                                        with fields t and x
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
%% Revision History:
%  Michael Casey                                                08/06/2026
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

%% Rebuild the control history on the returned grid:
                uk = zeros(numel(tk),2);
    for kt = 1:numel(tk)
         uk(kt,:) = ph.guide(tk(kt),xk(kt,:)').';
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

%% Record the junction and carry the state forward:
             xCurr = xk(end,:)';
              tOff = tOff + (tk(end) - tk(1));
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
