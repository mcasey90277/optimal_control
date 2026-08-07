function hFig = groundTrack(traj,veh,env,opts)
%% Purpose:
%
%  Plot the ground track of a trajectory on a latitude/longitude grid, one
%  coloured segment per flight phase, with the launch point and the impact
%  point marked and an aim point marked when one is supplied. Replaces the
%  ground-track figure that was written out by hand in HGV/run_glide,
%  BM/run_ballistic and HGV/run_boost_glide, which were the same plot three
%  times.
%
%  The state carries longitude and latitude in RADIANS. They are converted to
%  degrees once, here, because degrees are what a reader of a ground track
%  wants; both axis labels say so.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun. Reads
%                                               x(:,2) longitude (rad),
%                                               x(:,3) latitude (rad) and
%                                               phaseIdx [N x 1]
%
%  veh              Struct                      Vehicle parameters. NOT READ.
%                                               Present so that every function
%                                               in coorbital.viz takes the same
%                                               (traj,veh,env,opts) arguments
%                                               and a caller can hold one
%                                               argument list for all of them
%
%  env              Struct                      Environment model handles. NOT
%                                               READ, for the same reason. A
%                                               ground track is pure geometry:
%                                               no atmosphere, no gravity and
%                                               no vehicle enters it
%
%  opts             Struct, optional            All fields optional:
%                                               Parent    axes to draw into, or
%                                                         a figure/panel to
%                                                         create one in; [] or
%                                                         absent makes a new
%                                                         figure
%                                               Target    [2 x 1] aim point
%                                                         [lat; lon] (RAD),
%                                                         marked when given
%                                               PhaseName cellstr, one name per
%                                                         phase, for the legend
%                                               StartName char, legend text for
%                                                         the first sample;
%                                                         default 'launch'
%                                               EndName   char, legend text for
%                                                         the last sample;
%                                                         default 'impact'
%                                               Title     char, axes title
%                                               Visible   'on' (default) or
%                                                         'off'; applies only
%                                                         to a figure this
%                                                         function creates
%
%% Outputs:
%
%  hFig             [1 x 1] figure              The figure the track was drawn
%                                               in: the one created here, or
%                                               the one the Parent lives in
%
%% Notes:
%
%  LONGITUDE IS PLOTTED AS INTEGRATED, not wrapped into [-180,180]. The state
%  longitude is a continuous integral, so a flight crossing the antimeridian
%  simply continues past 180 deg. Wrapping it would put a full-width
%  horizontal streak across the figure at the crossing, which is an artefact
%  of the wrap and not something the vehicle did. The axis limits follow the
%  data, so the grid is still a latitude/longitude grid; it is just not always
%  centred on the prime meridian.
%
%  THE END POINTS ARE NAMED BY THE CALLER. They default to 'launch' and
%  'impact' because two of the three shipped entry scripts fly from a pad to
%  the ground, but HGV/run_glide begins at an entry interface and ends at a
%  terminal altitude, and calling either of those a launch or an impact would
%  be wrong. The MARKER TAGS stay launchMarker and impactMarker whatever the
%  legend says, because a tag is a programmatic handle and renaming it with
%  the label would break every caller that looks one up.
%
%  PHASE SEGMENTS AND THE JUNCTION SAMPLE. coorbital.prop.phaseRun records the
%  sample at a phase boundary ONCE and labels it with the OUTGOING phase, so
%  the first sample labelled phase p sits one integration step past the
%  junction. Drawing each phase over its own samples alone would therefore
%  leave a one-segment hole at every boundary. Each segment is instead
%  extended BACKWARD by the last sample of the preceding phase. That closes
%  the hole without inventing a point: the sample prepended is the junction
%  itself, which belongs to both segments as a matter of geometry even though
%  the phase index can only name one owner.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A synthetic two-phase great-circle arc, so the demo costs no
%% integration and still exercises the phase colouring and every marker:
if nargin == 0
                 c = coorbital.util.missileConst();
              nDem = 200;
              lam  = linspace(0,deg2rad(60),nDem)';
           demo.x  = [(c.rE + 50e3).*ones(nDem,1), ...
                      lam, ...
                      deg2rad(20).*sin(lam), ...
                      6000.*ones(nDem,1), ...
                      zeros(nDem,1), ...
                      zeros(nDem,1)];
           demo.t  = linspace(0,1800,nDem)';
           demo.u  = zeros(nDem,2);
    demo.phaseIdx  = 1 + (lam > deg2rad(30));
    demo.junction  = repmat(struct('t',[],'x',[]),1,1);
          demoOpt  = struct('Target',[deg2rad(19); deg2rad(61)], ...
                            'PhaseName',{{'glide','descent'}}, ...
                            'Title','coorbital.viz.groundTrack self-demo');
              hFig = coorbital.viz.groundTrack(demo, ...
                         coorbital.util.vehicleDefaults(),struct(),demoOpt);
    return;
end

%% Options, every one of them defaulted before anything is drawn:
    if nargin < 4 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
            parent = vizOption(opts,'Parent',[]);
            target = vizOption(opts,'Target',[]);
            phName = vizOption(opts,'PhaseName',{});
            titTxt = vizOption(opts,'Title','Ground track');
            visTxt = vizOption(opts,'Visible','on');
            begTxt = vizOption(opts,'StartName','launch');
            endTxt = vizOption(opts,'EndName','impact');

%% The plotted quantities. ONE conversion from radians to degrees, and the
%% axis labels below state the unit it produced:
            lonDeg = rad2deg(traj.x(:,2));
            latDeg = rad2deg(traj.x(:,3));
            phIdx  = traj.phaseIdx;
    assert(numel(phIdx) == numel(lonDeg), ...
        'phaseIdx has %d entries against %d samples.',numel(phIdx),numel(lonDeg));
            phList = unique(phIdx(:)).';
              nPh  = numel(phList);

%% Axes to draw into, and the only place a figure can be created:
        [hFig,hAx] = vizParent(parent,1,'Ground track',visTxt);
           wasHeld = ishold(hAx);
    hold(hAx,'on');

%% One line per phase, each extended backward to the junction so the track has
%% no holes in it. See the note in the header:
               col = lines(max(nPh,7));
            hTrack = gobjects(1,nPh);
            trkTxt = cell(1,nPh);
    for kp = 1:nPh
               sel = find(phIdx == phList(kp));
        if kp > 1
               sel = [find(phIdx == phList(kp-1),1,'last'); sel];
        end
        hTrack(kp) = line(hAx,lonDeg(sel),latDeg(sel), ...
                          'Color',col(kp,:), ...
                          'LineWidth',1.5, ...
                          'Tag','groundTrack', ...
                          'UserData',phList(kp));
         trkTxt{kp} = sprintf('phase %d',phList(kp));
        if numel(phName) >= kp
         trkTxt{kp} = phName{kp};
        end
    end

%% Launch and impact. Marked from the FIRST and LAST samples of the flown
%% trajectory, never from the user block that set them up, so a marker sitting
%% off the end of the track is evidence of a real disagreement:
             hLnch = line(hAx,lonDeg(1),latDeg(1), ...
                          'LineStyle','none','Marker','o','MarkerSize',9, ...
                          'LineWidth',1.5,'Color',[0 0.45 0], ...
                          'Tag','launchMarker');
              hImp = line(hAx,lonDeg(end),latDeg(end), ...
                          'LineStyle','none','Marker','s','MarkerSize',9, ...
                          'LineWidth',1.5,'Color',[0.75 0 0], ...
                          'Tag','impactMarker');
             hLeg  = [hTrack hLnch hImp];
            legTxt = [trkTxt {begTxt,endTxt}];

%% The aim point, when the caller has one. Distinct from the impact marker on
%% purpose: the gap between the two IS the miss distance, and a plot that drew
%% them the same way would hide it:
    if ~isempty(target)
        assert(numel(target) == 2, ...
            'Target must be [lat; lon] in radians, %d element(s) given.',numel(target));
            tgtLat = rad2deg(target(1));
            tgtLon = rad2deg(target(2));
             hTgt  = line(hAx,tgtLon,tgtLat, ...
                          'LineStyle','none','Marker','p','MarkerSize',13, ...
                          'LineWidth',1.5,'Color',[0 0 0.8], ...
                          'Tag','targetMarker');
             hLeg  = [hLeg hTgt];
            legTxt = [legTxt {'target'}];
    end

%% Labels. Both carry the unit the data was converted into:
    grid(hAx,'on');
    xlabel(hAx,'longitude (deg)');
    ylabel(hAx,'latitude (deg)');
    title(hAx,titTxt);
    legend(hAx,hLeg,legTxt,'Location','best');

%% Leave the axes' hold state as it was found, so composing into a caller's
%% axes has no side effect beyond the children added:
    if ~wasHeld
        hold(hAx,'off');
    end
end
