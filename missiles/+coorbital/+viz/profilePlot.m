function hFig = profilePlot(traj,veh,env,opts)
%% Purpose:
%
%  Plot selected time histories of a trajectory, one panel each, against one
%  shared time axis, with the phase boundaries drawn across every panel.
%  Replaces the two time-history figures each of HGV/run_glide,
%  BM/run_ballistic and HGV/run_boost_glide wrote out by hand.
%
%  Seven channels are available and five are drawn by default. The choice is
%  made with ONE option, opts.Channels, naming them in the order they should
%  appear; see the pool below. A caller with a series this library cannot
%  derive -- the sensed load factor including thrust is the case that forced
%  the feature -- appends it through opts.Extra rather than asking for a new
%  option here.
%
%  The library works in metres, metres per second, pascals and radians.
%  Several panels display a different unit -- kilometres, kilopascals,
%  standard gravities, degrees -- and each conversion is applied once,
%  explicitly, at the point of plotting, with the displayed unit named in the
%  axis label.
%
%% The channel pool, by name:
%
%    'altitude'   radius above the reference sphere, km
%    'speed'      planet-relative speed, m/s
%    'mach'       Mach number from the atmosphere's sound speed, -
%    'q'          dynamic pressure, kPa
%    'nAero'      AERODYNAMIC load factor, g. See the note below
%    'mass'       the mass STATE x(:,7), kg. Requires a seven-state trajectory
%    'gamma'      flight path angle, deg
%
%  The default is the first five, which is what every caller wanted before
%  'mass' and 'gamma' were added and is what a six-state glide can draw.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun. Reads
%                                               t [N x 1] (s), x(:,1) radius
%                                               (m), x(:,4) speed (m/s),
%                                               x(:,5) flight path angle
%                                               (rad), x(:,7) mass (kg) when
%                                               the state carries one, u(:,1)
%                                               angle of attack (rad),
%                                               phaseIdx [N x 1]
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
%                                               sound speed, always; and aero,
%                                               for the lift and drag
%                                               coefficients, only when the
%                                               'nAero' channel is drawn
%
%  opts             Struct, optional            All fields optional:
%                                               Channels cellstr naming the
%                                                        channels to draw, in
%                                                        order. Default
%                                                        {'altitude','speed',
%                                                        'mach','q','nAero'}
%                                               Extra    {series,ylabel,title}
%                                                        one caller-supplied
%                                                        [N x 1] series drawn
%                                                        as a FINAL panel. The
%                                                        label must name its
%                                                        unit in parentheses
%                                               Parent   axes array with one
%                                                        axes per panel, or a
%                                                        figure/panel to create
%                                                        them in; [] or absent
%                                                        makes a new figure
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
%  THE PANEL COUNT IS numel(Channels) + ~isempty(Extra), and it is fixed
%  before anything is drawn. A caller passing its own axes must supply exactly
%  that many, and is told the number if it does not.
%
%  'nAero' IS AERODYNAMIC AND EXCLUDES THRUST. It is the magnitude of the lift
%  and drag accelerations together, over standard gravity -- what the airframe
%  feels from the air. On a powered phase an accelerometer would also read the
%  thrust, and this channel does not show it, because the thrust term needs a
%  propulsion model, an ambient pressure and a knowledge of which phases are
%  burning, none of which a trajectory carries.
%
%  THAT IS WHAT opts.Extra IS FOR. A caller that has already computed the
%  sensed load -- BM/run_ballistic computes it for its own summary -- hands
%  the series over and gets it as a panel, rather than this function growing a
%  propulsion model it has no business owning. It is a DISPLAY channel: this
%  function does not check what the series means, only that it is the right
%  length and that its label states a unit.
%
%  'nAero' IS ALSO NOT A DECELERATION. With a constant-L/D aero model the load
%  is mostly lift, and only the drag part brakes the vehicle. The axis says
%  "aerodynamic load factor", not "deceleration", for exactly that reason.
%
%  WHY opts.VehPhase EXISTS. coorbital.prop.phaseRun carries one vehicle for a
%  whole chain, but a boosted flight is a stack under power and a separated
%  body afterwards, with different reference areas and different coefficients
%  -- see the "three vehicles, one chain" note in HGV/run_boost_glide. 'nAero'
%  is the only channel that depends on the airframe, so it is the only one
%  that needs to know. Given nothing, this function uses veh for every phase,
%  which is right for a single-airframe run and wrong, silently, for a staged
%  one; a staged caller must pass VehPhase.
%
%  THE MASS USED IN 'nAero' IS THE STATE MASS when the state carries one, and
%  the phase vehicle's mass otherwise. That is not a preference, it is what
%  makes the channel correct: during a burn the carried mass falls
%  continuously while veh.mass cannot, and after burnout the two agree
%  exactly, which the entry scripts assert. The branch also matters because a
%  booster struct legitimately has NO mass field at all -- see
%  coorbital.util.boosterDefaults, whose mass is split into massDry and
%  massProp -- so reaching for vehPhase{k}.mass unconditionally would fail on
%  exactly the staged runs the mass column exists to serve.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Channel pool; caller-supplied Extra panel     08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A short real propagation, because every channel but altitude is
%% a derived quantity and a synthetic state would not exercise the derivation:
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
                         struct('Name','coorbital.viz.profilePlot self-demo', ...
                                'Channels',{{'altitude','speed','mach','q', ...
                                             'nAero','gamma'}}));
    return;
end

%% Options, every one of them defaulted before anything is computed:
    if nargin < 4 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
            parent = vizOption(opts,'Parent',[]);
          vehPhase = vizOption(opts,'VehPhase',{});
            namTxt = vizOption(opts,'Name','Trajectory profile');
            visTxt = vizOption(opts,'Visible','on');
             extra = vizOption(opts,'Extra',{});
              want = vizOption(opts,'Channels', ...
                         {'altitude','speed','mach','q','nAero'});
    if ischar(want)
              want = {want};
    end
    assert(iscellstr(want) && ~isempty(want), ...
        'Channels must be a non-empty cellstr naming channels to draw.');

    assert(isfield(env,'atmos'), ...
        'env must carry atmos: the Mach and dynamic pressure channels come from it.');

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
                hM = traj.x(:,1) - c.rE;
                 V = traj.x(:,4);
    [rho,~,~,aSnd] = env.atmos(hM);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;

%% The channel pool. Series, axis label WITH its unit, and panel title. Only
%% the aerodynamic load costs anything to build, so it is built only if asked
%% for, and only then is env.aero required:
     pool.altitude = {hM./1000            ,'altitude (km)'              ,'Altitude'};
        pool.speed = {V                   ,'speed (m/s)'                ,'Speed'};
         pool.mach = {mach                ,'Mach number (-)'            ,'Mach number'};
            pool.q = {qbar./1000          ,'dynamic pressure (kPa)'     ,'Dynamic pressure'};
        pool.gamma = {rad2deg(traj.x(:,5)),'flight path angle (deg)'    ,'Flight path angle'};
    if any(strcmp('mass',want))
        assert(size(traj.x,2) >= 7, ...
            ['the ''mass'' channel needs the mass state x(:,7), and this ' ...
             'trajectory carries %d state(s). A six-state unpowered chain ' ...
             'has no mass to plot; drop ''mass'' from Channels, or run the ' ...
             'chain seven-state through coorbital.eom.massConstant.'], ...
            size(traj.x,2));
     pool.mass     = {traj.x(:,7)         ,'mass (kg)'                  ,'Mass'};
    end
    if any(strcmp('nAero',want))
        assert(isfield(env,'aero'), ...
            ['env must carry aero for the ''nAero'' channel; drop it from ' ...
             'Channels if the environment has no aerodynamic model.']);
             nAero = aeroLoad(traj,vehPhase,phList,qbar,mach,env,c);
     pool.nAero    = {nAero               ,'aerodynamic load factor (g)','Aerodynamic load factor'};
    end

%% Resolve the requested names against the pool. An unknown name is an error
%% and not a silently skipped panel, for the same reason a misspelt USER
%% PARAMETERS override is an error in the entry scripts:
             known = {'altitude','speed','mach','q','nAero','mass','gamma'};
              chan = cell(numel(want),3);
    for kc = 1:numel(want)
        assert(any(strcmp(want{kc},known)), ...
            '"%s" is not a channel. The pool is: %s.', ...
            want{kc},strjoin(known,', '));
        chan(kc,:) = pool.(want{kc});
    end

%% The caller's own panel, appended last. Checked for length and for a stated
%% unit, and for nothing else: what the series MEANS is the caller's business,
%% which is exactly why it is not a channel:
    if ~isempty(extra)
        assert(iscell(extra) && numel(extra) == 3, ...
            'Extra must be {series, ylabel, title}; %d element(s) given.',numel(extra));
        assert(isnumeric(extra{1}) && numel(extra{1}) == nS, ...
            ['the Extra series has %d points against %d trajectory samples; ' ...
             'it must be sampled on traj.t.'],numel(extra{1}),nS);
        assert(ischar(extra{2}) && contains(extra{2},'(') && contains(extra{2},')'), ...
            ['the Extra label "%s" must state its unit in parentheses, as ' ...
             'every other axis in this package does. A panel whose unit a ' ...
             'reader has to guess is the defect the labels exist to prevent.'], ...
            extra{2});
        assert(ischar(extra{3}) && ~isempty(strtrim(extra{3})), ...
            'the Extra panel needs a non-empty title.');
              chan = [chan; {extra{1}(:), extra{2}, extra{3}}];
    end
              nPan = size(chan,1);

%% Axes to draw into, and the only place a figure can be created:
        [hFig,hAx] = vizParent(parent,nPan,namTxt,visTxt);

%% Phase boundary times, drawn across every panel so a reader can tell which
%% phase a feature belongs to without a second figure. This is also what
%% replaced HGV/run_boost_glide's hand-drawn handoff marker, so the POSITION
%% of these lines carries information and is asserted in tests/test_viz.m:
              tBnd = zeros(1,0);
    for kp = 2:numel(phList)
             kLast = find(phIdx == phList(kp-1),1,'last');
              tBnd = [tBnd traj.t(kLast)];
    end

    for ka = 1:nPan
           wasHeld = ishold(hAx(ka));
        hold(hAx(ka),'on');
        for kb = 1:numel(tBnd)
            line(hAx(ka),[tBnd(kb) tBnd(kb)],[min(chan{ka,1}) max(chan{ka,1})], ...
                 'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1, ...
                 'Tag','phaseBoundary');
        end
        line(hAx(ka),traj.t,chan{ka,1}, ...
             'Color',[0 0.298 0.6],'LineWidth',1.5,'Tag','profile');
        grid(hAx(ka),'on');
        xlabel(hAx(ka),'time (s)');
        ylabel(hAx(ka),chan{ka,2});
        title(hAx(ka),chan{ka,3});
        if ~wasHeld
            hold(hAx(ka),'off');
        end
    end
end

function nAero = aeroLoad(traj,vehPhase,phList,qbar,mach,env,c)
%% Purpose:
%
%  Aerodynamic load factor along the trajectory: the magnitude of the lift and
%  drag accelerations together, over standard gravity. Thrust is excluded --
%  see the notes in the caller's header. Each sample is evaluated with ITS OWN
%  phase's airframe, so a staged chain reports the stack under power and the
%  separated body afterwards rather than one of them throughout.
%
%% Inputs:
%
%  traj             Struct                      Trajectory; reads x, u and
%                                               phaseIdx
%
%  vehPhase         Cell [1 x P]                One vehicle struct per phase;
%                                               each needs Sref (m^2), and
%                                               mass (kg) only when the state
%                                               carries no mass column
%
%  phList           [1 x P]                     Phase indices present, sorted
%
%  qbar             [N x 1]                     Dynamic pressure (Pa)
%
%  mach             [N x 1]                     Mach number (-)
%
%  env              Struct                      Environment; reads aero
%
%  c                Struct                      Constants; uses g0 (m/s^2)
%
%% Outputs:
%
%  nAero            [N x 1]                     Aerodynamic load factor (g)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = numel(traj.t);
             phIdx = traj.phaseIdx;

%% Mass actually carried. The state mass when there is one -- it falls through
%% a burn, and veh.mass cannot -- and the phase vehicle's mass otherwise:
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

             aLift = zeros(nS,1);
             aDrag = zeros(nS,1);
    for ks = 1:nS
               kPh = find(phList == phIdx(ks),1);
              vehK = vehPhase{kPh};
         [CLk,CDk] = env.aero(traj.u(ks,1),mach(ks),vehK);
         aLift(ks) = qbar(ks).*vehK.Sref.*CLk./massV(ks);
         aDrag(ks) = qbar(ks).*vehK.Sref.*CDk./massV(ks);
    end
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;
end
