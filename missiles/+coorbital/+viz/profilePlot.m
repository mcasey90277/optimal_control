function hFig = profilePlot(traj,veh,env,opts)
%% Purpose:
%
%  Plot the five time histories that describe what a trajectory did to the
%  vehicle -- altitude, speed, Mach number, dynamic pressure and aerodynamic
%  load factor -- one panel each, against one shared time axis, with the phase
%  boundaries drawn across every panel. Replaces the two time-history figures
%  each of HGV/run_glide, BM/run_ballistic and HGV/run_boost_glide wrote out
%  by hand.
%
%  The library works in metres, metres per second, pascals and radians. Three
%  of the five panels display a different unit -- kilometres, kilopascals and
%  standard gravities -- and each conversion is applied once, explicitly, at
%  the point of plotting, with the displayed unit named in the axis label.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun. Reads
%                                               t [N x 1] (s), x(:,1) radius
%                                               (m), x(:,4) speed (m/s),
%                                               x(:,7) mass (kg) when the state
%                                               carries one, u(:,1) angle of
%                                               attack (rad), phaseIdx [N x 1]
%
%  veh              Struct                      Vehicle flown, for the load
%                                               factor: uses mass (kg) and
%                                               Sref (m^2). Chains that fly
%                                               more than one airframe pass
%                                               them through opts.VehPhase and
%                                               this argument is then unread
%
%  env              Struct                      Environment model handles.
%                                               Reads atmos, for density and
%                                               sound speed, and aero, for the
%                                               lift and drag coefficients.
%                                               Both are required
%
%  opts             Struct, optional            All fields optional:
%                                               Parent   [1 x 5] axes array to
%                                                        draw into, or a
%                                                        figure/panel to create
%                                                        five axes in; [] or
%                                                        absent makes a new
%                                                        figure
%                                               VehPhase cell array, one
%                                                        vehicle struct per
%                                                        phase; default is veh
%                                                        for every phase
%                                               Name     char, figure Name
%                                               Visible  'on' (default) or
%                                                        'off'; applies only to
%                                                        a figure this function
%                                                        creates
%
%% Outputs:
%
%  hFig             [1 x 1] figure              The figure the panels were
%                                               drawn in: the one created here,
%                                               or the one the Parent lives in
%
%% Notes:
%
%  THE LOAD FACTOR IS AERODYNAMIC AND EXCLUDES THRUST. It is the magnitude of
%  the lift and drag accelerations together, over standard gravity -- what the
%  airframe feels from the air. On a powered phase the accelerometer would
%  also read the thrust, and this panel does not show it. Two reasons. The
%  thrust term needs a propulsion model, an ambient pressure and a knowledge
%  of which phases are burning, none of which a trajectory carries; and the
%  aerodynamic load is the quantity that is comparable across every phase of
%  every chain in this library. An entry script that wants the sensed load
%  including thrust already computes and reports it in its own summary.
%
%  IT IS ALSO NOT A DECELERATION. With a constant-L/D aero model the load is
%  mostly lift, and only the drag part brakes the vehicle. The axis says
%  "aerodynamic load factor", not "deceleration", for exactly that reason.
%
%  WHY opts.VehPhase EXISTS. coorbital.prop.phaseRun carries one vehicle for a
%  whole chain, but a boosted flight is a stack under power and a separated
%  body afterwards, with different reference areas and different coefficients
%  -- see the "three vehicles, one chain" note in HGV/run_boost_glide. The
%  load factor is the only panel that depends on the airframe, so it is the
%  only one that needs to know. Given nothing, this function uses veh for
%  every phase, which is right for a single-airframe run and wrong, silently,
%  for a staged one; a staged caller must pass VehPhase.
%
%  THE MASS USED IS THE STATE MASS when the state carries one, and the phase
%  vehicle's mass otherwise. That is not a preference, it is what makes the
%  panel correct: during a burn the carried mass falls continuously while
%  veh.mass cannot, and after burnout the two agree exactly, which the entry
%  scripts assert.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A short real propagation, because every panel but altitude is a
%% derived quantity and a synthetic state would not exercise the derivation:
if nargin == 0
                 c = coorbital.util.missileConst();
            vehDem = coorbital.util.vehicleDefaults();
      envDem.atmos = @coorbital.atmos.expAtmos;
       envDem.grav = @coorbital.grav.sphereGrav;
       envDem.aero = @coorbital.aero.constLD;
     envDem.omegaE = 0;
             sched = struct('tGrid',[0 2000],'alpha',[0 0],'sigma',[0 0]);
           phD.eom = @coorbital.eom.glide3DOF;
         phD.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
     phD.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
         phD.tspan = [0 2000];
             x0Dem = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
            trjDem = coorbital.prop.phaseRun(phD,x0Dem,vehDem,envDem);
              hFig = coorbital.viz.profilePlot(trjDem,vehDem,envDem, ...
                         struct('Name','coorbital.viz.profilePlot self-demo'));
    return;
end

%% Options, every one of them defaulted before anything is computed:
    if nargin < 4 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
            parent = optionOf(opts,'Parent',[]);
          vehPhase = optionOf(opts,'VehPhase',{});
            namTxt = optionOf(opts,'Name','Trajectory profile');
            visTxt = optionOf(opts,'Visible','on');

    assert(isfield(env,'atmos') && isfield(env,'aero'), ...
        ['env must carry both atmos and aero: the Mach, dynamic pressure and ' ...
         'load panels are derived from them.']);

                 c = coorbital.util.missileConst();
                nS = numel(traj.t);
             phIdx = traj.phaseIdx;
            phList = unique(phIdx(:)).';

%% One vehicle per phase. Given none, the single vehicle flies the whole
%% chain, which is right for an unstaged run and is what the caller has
%% asserted by not supplying the list:
    if isempty(vehPhase)
          vehPhase = repmat({veh},1,numel(phList));
    end
    assert(iscell(vehPhase) && numel(vehPhase) == numel(phList), ...
        ['VehPhase must hold one vehicle per phase: %d given against %d ' ...
         'phases in the trajectory.'],numel(vehPhase),numel(phList));

%% Atmosphere along the flight, evaluated once:
               hM  = traj.x(:,1) - c.rE;
                 V = traj.x(:,4);
    [rho,~,~,aSnd] = env.atmos(hM);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;

%% Mass actually carried. The state mass when there is one -- it falls through
%% a burn, and veh.mass cannot -- and the phase vehicle's mass otherwise. The
%% branch is not an optimisation: a booster struct legitimately has NO mass
%% field at all, because its mass is split into massDry and massProp and only
%% the chain that assembles them knows the total. Reading vehPhase{kp}.mass
%% unconditionally would therefore fail on exactly the staged runs the mass
%% column exists to serve:
    if size(traj.x,2) >= 7
             massV = traj.x(:,7);
    else
             massV = zeros(nS,1);
        for kp = 1:numel(phList)
               sel = phIdx == phList(kp);
            assert(isfield(vehPhase{kp},'mass'), ...
                ['phase %d''s vehicle carries no mass field and the state ' ...
                 'carries no mass column, so the load factor has nothing to ' ...
                 'divide by.'],phList(kp));
        massV(sel) = vehPhase{kp}.mass;
        end
    end

%% Aerodynamic accelerations sample by sample, each with ITS OWN phase's
%% airframe, so a staged chain reports the stack under power and the separated
%% body afterwards rather than one of them throughout:
             aLift = zeros(nS,1);
             aDrag = zeros(nS,1);
    for ks = 1:nS
              kPh  = find(phList == phIdx(ks),1);
              vehK = vehPhase{kPh};
         [CLk,CDk] = env.aero(traj.u(ks,1),mach(ks),vehK);
         aLift(ks) = qbar(ks).*vehK.Sref.*CLk./massV(ks);
         aDrag(ks) = qbar(ks).*vehK.Sref.*CDk./massV(ks);
    end
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;

%% The five channels, each already in the unit its label will name. This is
%% the only place a unit conversion happens:
              chan = {traj.t, hM./1000       ,'altitude (km)'             ,'Altitude'; ...
                      traj.t, V              ,'speed (m/s)'               ,'Speed'; ...
                      traj.t, mach           ,'Mach number (-)'           ,'Mach number'; ...
                      traj.t, qbar./1000     ,'dynamic pressure (kPa)'    ,'Dynamic pressure'; ...
                      traj.t, nAero          ,'aerodynamic load factor (g)','Aerodynamic load factor'};
              nPan = size(chan,1);

%% Axes to draw into, and the only place a figure can be created:
        [hFig,hAx] = vizParent(parent,nPan,namTxt,visTxt);

%% Phase boundary times, drawn across every panel so a reader can tell which
%% phase a feature belongs to without a second figure:
             tBnd  = zeros(1,0);
    for kp = 2:numel(phList)
             kLast = find(phIdx == phList(kp-1),1,'last');
             tBnd  = [tBnd traj.t(kLast)];
    end

    for ka = 1:nPan
           wasHeld = ishold(hAx(ka));
        hold(hAx(ka),'on');
        for kb = 1:numel(tBnd)
            line(hAx(ka),[tBnd(kb) tBnd(kb)],[min(chan{ka,2}) max(chan{ka,2})], ...
                 'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1, ...
                 'Tag','phaseBoundary');
        end
        line(hAx(ka),chan{ka,1},chan{ka,2}, ...
             'Color',[0 0.298 0.6],'LineWidth',1.5,'Tag','profile');
        grid(hAx(ka),'on');
        xlabel(hAx(ka),'time (s)');
        ylabel(hAx(ka),chan{ka,3});
        title(hAx(ka),chan{ka,4});
        if ~wasHeld
            hold(hAx(ka),'off');
        end
    end
end

function val = optionOf(opts,name,val)
%% Purpose:
%
%  Return one option's value when the caller supplied it and the default
%  otherwise. Factored out so the option block above reads as one unbranched
%  line per option, matching the override blocks in the entry scripts.
%
%% Inputs:
%
%  opts             Struct                      Caller options
%
%  name             Char [1 x n]                Option name
%
%  val              Any                         Default value
%
%% Outputs:
%
%  val              Any                         opts.(name) when present,
%                                               otherwise the default
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if isfield(opts,name) && ~isempty(opts.(name))
               val = opts.(name);
    end
end
