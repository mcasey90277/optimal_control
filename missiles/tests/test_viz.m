function test_viz()
%% Purpose:
%
%  Pin what the coorbital.viz figures SAY. A plot test that draws a figure and
%  checks that nothing threw is worth almost nothing: every mistake worth
%  catching -- plotting the wrong column, converting metres to kilometres
%  twice, labelling an axis without its unit, marking a launch point that is
%  not where the flight began, drawing a phase boundary in the wrong place --
%  produces a figure that draws perfectly and is wrong. So this file reads the
%  drawn objects back out and compares them against the trajectory the figure
%  claims to be showing.
%
%  Four things are asserted for every figure in the package:
%
%    1. The axes count is what the function promises.
%    2. Every axes carries a non-empty XLabel, YLabel and Title.
%    3. Every axis label names the UNIT that axis is displaying, which is the
%       only thing standing between a reader and a silent unit error.
%    4. The XData and YData of the primary line equal the trajectory arrays,
%       recomputed here from traj rather than read back from the function.
%
%  The fourth is the one that makes the file a test. The recomputations below
%  duplicate the arithmetic in the plotting routines ON PURPOSE, in the same
%  units, so that a wrong column, a stale array or a doubled conversion moves
%  one side and not the other.
%
%  THE PHASE BOUNDARY LINES ARE CHECKED FOR POSITION, not just counted.
%  HGV/run_boost_glide's hand-drawn handoff marker was DELETED in favour of
%  them, so where they sit is now load-bearing information rather than
%  decoration, and a boundary drawn at the wrong instant is a figure that lies
%  about which phase a feature belongs to. Counting them alone did not catch
%  a mutation that moved every one of them to t(1); section 9 does.
%
%  THE TRAJECTORIES ARE SYNTHETIC, and deliberately so. The figures do not
%  integrate anything, so a real propagation would add seconds to the suite
%  and buy nothing; what a synthetic state buys instead is control over
%  phaseIdx, which is how the phase-segment and phase-boundary contracts below
%  get exercised at all. The composition of viz with the real propagator is
%  covered where it belongs, in the entry-script tests.
%
%  Every figure is created with 'Visible','off' and closed again, and the
%  figure count is asserted back to its starting value at the end, so the
%  suite neither displays nor leaks.
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
%  Michael Casey                                                08/07/2026
%  Michael Casey  Channel pool, Extra panel, boundary positions 08/07/2026
%  Michael Casey  globeMovie: (traj,veh,env,opts), marker shell,
%                 phase width grading, handoff markers, legend  08/07/2026
%  Michael Casey  globeMovie: the altitude inset, the arc-frame
%                 camera, the spin centring and the two
%                 degenerate arcs -- all of which shipped with
%                 nothing at all asserted about them            08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
             nFig0 = numel(findall(groot,'Type','figure'));

%% The environment the derived panels are computed from, and two vehicles that
%% are genuinely different from one another. The difference is the point: a
%% load-factor panel that ignored the per-phase vehicle list and used the
%% single veh argument throughout would agree with the expectation below only
%% if the vehicles happened to match, so they are made not to.
%%
%% The boost vehicle is the real coorbital.util.boosterDefaults struct rather
%% than a doctored copy of the glide vehicle, and that is deliberate: a
%% booster has NO mass field -- its mass is split into massDry and massProp,
%% and only the chain that assembles them knows the total. A load-factor panel
%% that reached for veh.mass unconditionally would throw on it, which is
%% exactly the defect this argument caught during construction:
         env.atmos = @coorbital.atmos.expAtmos;
          env.aero = @coorbital.aero.constLD;
          env.grav = @coorbital.grav.sphereGrav;
        env.omegaE = 0;
            vehOne = coorbital.util.vehicleDefaults();
            vehBst = coorbital.util.boosterDefaults();
    assert(~isfield(vehBst,'mass'), ...
        ['boosterDefaults has grown a mass field; this test relied on its ' ...
         'absence to prove profilePlot does not need one.']);
          vehPhase = {vehBst,vehOne,vehOne};

%% A three-phase, seven-state trajectory, and a one-phase six-state one. Both
%% shapes ship in this library and they exercise different branches: the mass
%% column, the mass CHANNEL, and the per-phase vehicle list:
            trjCha = chainTraj(c.rE);
            trjOne = singleTraj(c.rE);
            phList = unique(trjCha.phaseIdx(:)).';

%% ---------------------------------------------------------------------
%% 1. groundTrack -- structure, labels and units
%% ---------------------------------------------------------------------
            target = [deg2rad(35.5); deg2rad(-77.0)];
            hFigGT = coorbital.viz.groundTrack(trjCha,vehOne,env, ...
                         struct('Visible','off','Target',target, ...
                                'PhaseName',{{'boost','glide','descent'}}, ...
                                'Title','Ground track, 7663 km great-circle range'));
              axGT = axesInOrder(hFigGT);
    assert(numel(axGT) == 1, ...
        'groundTrack drew %d axes; it promises exactly 1.',numel(axGT));
    assertLabelled(axGT,'groundTrack');
    assertUnit(axGT(1).XLabel.String,'deg','groundTrack x axis');
    assertUnit(axGT(1).YLabel.String,'deg','groundTrack y axis');
    assertUnit(axGT(1).XLabel.String,'longitude','groundTrack x axis');
    assertUnit(axGT(1).YLabel.String,'latitude','groundTrack y axis');

%% ---------------------------------------------------------------------
%% 2. groundTrack -- the plotted data IS the trajectory
%% ---------------------------------------------------------------------
%% One line per phase, each covering its own samples with the preceding
%% phase's last sample prepended, which is the contract the function's header
%% states and the thing that keeps the track from breaking at every junction.
%% Reconstructed here from phaseIdx alone, so a segment drawn over the wrong
%% index range fails even though the union of the segments still looks right:
            hTrack = lineByTag(axGT(1),'groundTrack');
    assert(numel(hTrack) == numel(phList), ...
        'groundTrack drew %d track segments for %d phases.', ...
        numel(hTrack),numel(phList));
    for kp = 1:numel(phList)
              hSeg = hTrack([hTrack.UserData] == phList(kp));
        assert(isscalar(hSeg), ...
            'expected exactly one track segment tagged for phase %d.',phList(kp));
               sel = phaseSpan(trjCha.phaseIdx,phList,kp);
        assertSeries(hSeg.XData(:),rad2deg(trjCha.x(sel,2)), ...
            sprintf('groundTrack phase %d longitude (deg)',phList(kp)));
        assertSeries(hSeg.YData(:),rad2deg(trjCha.x(sel,3)), ...
            sprintf('groundTrack phase %d latitude (deg)',phList(kp)));
    end

%% ...and the markers, each against the sample it claims to mark:
    assertMarker(axGT(1),'launchMarker', ...
        rad2deg(trjCha.x(1,2)),rad2deg(trjCha.x(1,3)),'launch');
    assertMarker(axGT(1),'impactMarker', ...
        rad2deg(trjCha.x(end,2)),rad2deg(trjCha.x(end,3)),'impact');
    assertMarker(axGT(1),'targetMarker', ...
        rad2deg(target(2)),rad2deg(target(1)),'target');

%% The launch and impact markers must not be the same point, or the two
%% assertions above would both pass on a figure that marked one thing twice:
             dMark = abs(rad2deg(trjCha.x(1,2)) - rad2deg(trjCha.x(end,2)));
    assert(dMark > 1, ...
        ['the synthetic track starts and ends within %.3f deg of longitude; ' ...
         'the marker checks cannot discriminate.'],dMark);

%% The default legend names the ends launch and impact:
    assertLegend(hFigGT,{'boost','glide','descent','launch','impact','target'});
    close(hFigGT);

%% ---------------------------------------------------------------------
%% 3. groundTrack -- optional target, and caller-named end points
%% ---------------------------------------------------------------------
%% The aim point is optional, and an optional marker that is drawn anyway --
%% at the origin, or at the impact point -- is worse than no marker at all,
%% because a reader would read the miss distance off it.
%%
%% The end-point NAMES are the caller's, because HGV/run_glide begins at an
%% entry interface and ends at a terminal altitude, and neither is a launch or
%% an impact. The marker TAGS must not follow the names, or every caller that
%% looks a marker up by tag would break the moment a script renamed one:
            hFigNT = coorbital.viz.groundTrack(trjCha,vehOne,env, ...
                         struct('Visible','off', ...
                                'StartName','entry','EndName','terminus'));
              axNT = axesInOrder(hFigNT);
    assert(isempty(findobj(axNT(1),'Tag','targetMarker')), ...
        'groundTrack drew a target marker when no Target was supplied.');
    assert(~isempty(findobj(axNT(1),'Tag','launchMarker')) && ...
           ~isempty(findobj(axNT(1),'Tag','impactMarker')), ...
        ['launch and impact must be marked, under their fixed tags, whether ' ...
         'or not a target is given and whatever the legend calls them.']);
    assertLegend(hFigNT,{'phase 1','phase 2','phase 3','entry','terminus'});
    close(hFigNT);

%% ---------------------------------------------------------------------
%% 4. profilePlot -- the DEFAULT five channels, structure, labels and units
%% ---------------------------------------------------------------------
            hFigPP = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Visible','off','VehPhase',{vehPhase}));
              axPP = axesInOrder(hFigPP);
    assert(numel(axPP) == 5, ...
        'profilePlot drew %d axes; its default is exactly 5 channels.',numel(axPP));
    assertLabelled(axPP,'profilePlot');

%% Every panel shares one time axis, and every panel's own unit must appear in
%% its own label. The list is written out rather than looped over a pattern so
%% that a panel silently relabelled into another panel's unit still fails:
           dfltChn = {'altitude','speed','mach','q','nAero'};
    for ka = 1:5
        assertUnit(axPP(ka).XLabel.String,'(s)', ...
            sprintf('profilePlot panel %d x axis',ka));
        assertUnit(axPP(ka).YLabel.String,unitOf(dfltChn{ka}), ...
            sprintf('profilePlot panel %d y axis',ka));
    end

%% ---------------------------------------------------------------------
%% 5. profilePlot -- every panel's data, recomputed from traj
%% ---------------------------------------------------------------------
             wantC = profileChannels(trjCha,vehPhase,c);
    assertPanels(axPP,dfltChn,wantC,trjCha.t,'default');

%% The five channels must be five DIFFERENT things. Without this, a function
%% that plotted altitude into all five panels could still pass every check
%% above if the expectations were ever built from the drawn data:
    assertDistinct(dfltChn,wantC);
    close(hFigPP);

%% ---------------------------------------------------------------------
%% 6. profilePlot -- the per-phase vehicle list is actually read
%% ---------------------------------------------------------------------
%% Run the same trajectory again with the single-vehicle default. The load
%% factor must MOVE, because the boost phase is then flown on the wrong
%% airframe. If it does not, VehPhase is being ignored and section 5 proved
%% nothing about it:
            hFigSV = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Visible','off'));
              axSV = axesInOrder(hFigSV);
             nSing = lineByTag(axSV(5),'profile');
    assert(max(abs(nSing.YData(:) - wantC.nAero)) > 1e-3, ...
        ['the load factor is unchanged when the per-phase vehicle list is ' ...
         'withheld; profilePlot is ignoring VehPhase.']);
    close(hFigSV);

%% ---------------------------------------------------------------------
%% 7. profilePlot -- a six-state, single-phase trajectory
%% ---------------------------------------------------------------------
%% No mass column and no staging: the mass must then come from the vehicle,
%% and there must be no phase boundary to draw:
             hFig6 = coorbital.viz.profilePlot(trjOne,vehOne,env, ...
                         struct('Visible','off'));
               ax6 = axesInOrder(hFig6);
    assert(numel(ax6) == 5,'profilePlot drew %d axes on a six-state run.',numel(ax6));
             want6 = profileChannels(trjOne,{vehOne},c);
    assertPanels(ax6,dfltChn,want6,trjOne.t,'six-state');
    assert(isempty(findobj(ax6(1),'Tag','phaseBoundary')), ...
        'a single-phase run has no phase boundary to draw.');
    close(hFig6);

%% ---------------------------------------------------------------------
%% 8. profilePlot -- the WHOLE channel pool, plus a caller-supplied panel
%% ---------------------------------------------------------------------
%% 'mass' and 'gamma' are the two channels the first cut of this package
%% dropped, and both were restored because the printed summary cannot show
%% what they show: mass carries the STAGING DISCONTINUITY, which is this
%% library's headline feature and is a jump rather than a number, and gamma
%% through a skip phugoid is a shape whose per-phase mean says nothing about
%% it.
%%
%% The Extra panel is the sensed load factor -- the thing profilePlot must not
%% compute, because thrust needs a propulsion model, an ambient pressure and a
%% list of which phases burn, and a trajectory carries none of them. It is
%% built here the way BM/run_ballistic builds it, from arbitrary but distinct
%% numbers, because what is under test is that the SERIES HANDED OVER is the
%% series drawn, not what it means:
            allChn = {'altitude','speed','mach','q','nAero','mass','gamma'};
            xtraSr = 3 + sin(linspace(0,9,numel(trjCha.t))');
           hFigAll = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Visible','off','VehPhase',{vehPhase}, ...
                                'Channels',{allChn}, ...
                                'Extra',{{xtraSr,'sensed load factor (g)', ...
                                          'Sensed load factor'}}));
             axAll = axesInOrder(hFigAll);
    assert(numel(axAll) == numel(allChn) + 1, ...
        ['profilePlot drew %d axes for %d channels plus an Extra panel; the ' ...
         'count must be numel(Channels) + ~isempty(Extra).'], ...
        numel(axAll),numel(allChn));
    assertLabelled(axAll,'profilePlot full pool');
    for ka = 1:numel(allChn)
        assertUnit(axAll(ka).YLabel.String,unitOf(allChn{ka}), ...
            sprintf('profilePlot ''%s'' panel',allChn{ka}));
    end
    assertPanels(axAll(1:numel(allChn)),allChn,wantC,trjCha.t,'full pool');
    assertDistinct(allChn,wantC);

%% The Extra panel is LAST, carries the caller's series unchanged, and wears
%% the caller's label and title:
              hXtr = lineByTag(axAll(end),'profile');
    assert(isscalar(hXtr),'the Extra panel holds %d primary lines.',numel(hXtr));
    assertSeries(hXtr.XData(:),trjCha.t,'Extra panel time (s)');
    assertSeries(hXtr.YData(:),xtraSr,'Extra panel series');
    assertUnit(axAll(end).YLabel.String,'(g)','Extra panel y axis');
    assert(contains(labelText(axAll(end).Title.String),'Sensed'), ...
        'the Extra panel does not wear the title it was given.');

%% ...and the Extra series must be nothing the pool already draws, or the
%% check above would pass on a panel that quietly plotted a channel:
    for kc = 1:numel(allChn)
        assert(max(abs(xtraSr - wantC.(allChn{kc}))) > 1e-6, ...
            'the Extra series coincides with the ''%s'' channel.',allChn{kc});
    end

%% ---------------------------------------------------------------------
%% 9. profilePlot -- WHERE the phase boundaries are, not just how many
%% ---------------------------------------------------------------------
%% This is the check that makes the dashed lines trustworthy. They replaced
%% HGV/run_boost_glide's deleted handoff marker, so their POSITION is the
%% information. The junction instants are recomputed here from phaseIdx alone:
              tBnd = zeros(1,0);
    for kp = 2:numel(phList)
              tBnd = [tBnd trjCha.t(find(trjCha.phaseIdx == phList(kp-1),1,'last'))];
    end
    assert(min(abs(tBnd - trjCha.t(1))) > 1, ...
        ['every junction sits within 1 s of the start of the trajectory; ' ...
         'this check cannot tell a correct boundary from one collapsed to ' ...
         't(1).']);
    for ka = 1:numel(axAll)
              hBnd = findobj(axAll(ka),'Tag','phaseBoundary');
        assert(numel(hBnd) == numel(tBnd), ...
            'panel %d drew %d phase boundaries for %d junctions.', ...
            ka,numel(hBnd),numel(tBnd));
               got = zeros(1,numel(hBnd));
        for kb = 1:numel(hBnd)
            assert(numel(hBnd(kb).XData) == 2 && ...
                   hBnd(kb).XData(1) == hBnd(kb).XData(2), ...
                'panel %d boundary %d is not a vertical line.',ka,kb);
           got(kb) = hBnd(kb).XData(1);
        end
        assertSeries(sort(got(:)),sort(tBnd(:)), ...
            sprintf('panel %d phase boundary times (s)',ka));
    end
    close(hFigAll);

%% ---------------------------------------------------------------------
%% 10. profilePlot -- what it must REFUSE
%% ---------------------------------------------------------------------
%% Each of these is a caller error that would otherwise produce a plausible
%% figure. A misspelt channel silently dropped, or an Extra panel whose unit
%% a reader has to guess, is the failure mode the whole package exists to
%% prevent:
    assertThrows(@() coorbital.viz.profilePlot(trjCha,vehOne,env, ...
        struct('Visible','off','Channels',{{'altitude','altitud'}})), ...
        'a misspelt channel name');
    assertThrows(@() coorbital.viz.profilePlot(trjOne,vehOne,env, ...
        struct('Visible','off','Channels',{{'altitude','mass'}})), ...
        'the mass channel on a six-state trajectory');
    assertThrows(@() coorbital.viz.profilePlot(trjCha,vehOne,env, ...
        struct('Visible','off','Extra',{{[1 2 3],'load (g)','Load'}})), ...
        'an Extra series of the wrong length');
    assertThrows(@() coorbital.viz.profilePlot(trjCha,vehOne,env, ...
        struct('Visible','off', ...
               'Extra',{{xtraSr,'sensed load','Load'}})), ...
        'an Extra label that states no unit');

%% ---------------------------------------------------------------------
%% 11. globe3D -- structure, labels, units and the arc
%% ---------------------------------------------------------------------
             altSc = 30;
            hFigG3 = coorbital.viz.globe3D(trjCha,vehOne,env, ...
                         struct('Visible','off','AltScale',altSc, ...
                                'Target',target));
              axG3 = axesInOrder(hFigG3);
    assert(numel(axG3) == 1,'globe3D drew %d axes; it promises exactly 1.',numel(axG3));
    assertLabelled(axG3,'globe3D');
    assertUnit(axG3(1).XLabel.String,'km','globe3D x axis');
    assertUnit(axG3(1).YLabel.String,'km','globe3D y axis');
    assertUnit(axG3(1).ZLabel.String,'km','globe3D z axis');

%% A non-unity exaggeration MUST be stated in the title. A vertical scale
%% thirty times the horizontal one, unannounced, is a lie told with a plot:
             titG3 = labelText(axG3(1).Title.String);
    assert(contains(titG3,'exaggerat') && contains(titG3,sprintf('%g',altSc)), ...
        ['globe3D was given AltScale = %g and its title does not say so. ' ...
         'Title was "%s".'],altSc,titG3);

%% The planet is drawn, once, by the shared renderer:
            hEarth = findobj(axG3(1),'Tag','earthSurface');
    assert(isscalar(hEarth),'globe3D drew %d Earth surfaces; expected 1.',numel(hEarth));
    assert(abs(max(hEarth.XData(:)) - c.rE./1000) < 1e-6, ...
        'the Earth surface has radius %.3f km against rE = %.3f km.', ...
        max(hEarth.XData(:)),c.rE./1000);

%% The arc, phase by phase, against the spherical-to-Cartesian conversion done
%% here independently. This is where an exaggeration applied to the sphere
%% instead of to the altitude, or a latitude/longitude transposition in the
%% conversion, is caught:
              hArc = lineByTag(axG3(1),'globeTrack');
    assert(numel(hArc) == numel(phList), ...
        'globe3D drew %d arc segments for %d phases.',numel(hArc),numel(phList));
              rPlt = c.rE + altSc.*(trjCha.x(:,1) - c.rE);
              lonC = trjCha.x(:,2);
              latC = trjCha.x(:,3);
             xWant = rPlt.*cos(latC).*cos(lonC)./1000;
             yWant = rPlt.*cos(latC).*sin(lonC)./1000;
             zWant = rPlt.*sin(latC)./1000;
    for kp = 1:numel(phList)
              hSeg = hArc([hArc.UserData] == phList(kp));
        assert(isscalar(hSeg),'expected one globe arc segment for phase %d.',phList(kp));
               sel = phaseSpan(trjCha.phaseIdx,phList,kp);
        assertSeries(hSeg.XData(:),xWant(sel),sprintf('globe3D phase %d x (km)',phList(kp)));
        assertSeries(hSeg.YData(:),yWant(sel),sprintf('globe3D phase %d y (km)',phList(kp)));
        assertSeries(hSeg.ZData(:),zWant(sel),sprintf('globe3D phase %d z (km)',phList(kp)));
    end

%% The arc really is lifted off the surface by the exaggeration, so the check
%% above is discriminating rather than comparing two things that coincide:
            liftKm = max(rPlt - c.rE)./1000;
    assert(liftKm > 100, ...
        'the exaggerated arc stands only %.3f km off the sphere.',liftKm);

%% ...and the surface markers are NOT lifted. They mark ground positions, and
%% an exaggerated launch marker sits somewhere the vehicle never was:
    assertSurfaceMarker(axG3(1),'launchMarker',c.rE,latC(1)  ,lonC(1)  ,'launch');
    assertSurfaceMarker(axG3(1),'impactMarker',c.rE,latC(end),lonC(end),'impact');
    assertSurfaceMarker(axG3(1),'targetMarker',c.rE,target(1),target(2),'target');
    close(hFigG3);

%% A unity exaggeration is the default and must NOT be announced, or every
%% honest figure would carry a caveat that does not apply to it:
            hFigG1 = coorbital.viz.globe3D(trjCha,vehOne,env,struct('Visible','off'));
              axG1 = axesInOrder(hFigG1);
    assert(~contains(labelText(axG1(1).Title.String),'exaggerat'), ...
        'globe3D announced an exaggeration it did not apply.');
    close(hFigG1);

%% ---------------------------------------------------------------------
%% 12. globeMovie -- a movie of the trajectory DEVELOPING
%% ---------------------------------------------------------------------
%% The one assertion that distinguishes a movie of a growing trajectory from a
%% static arc spun on a turntable: the track drawn in an INTERMEDIATE frame
%% must be the trajectory TRUNCATED at that frame's time, and the track drawn
%% in the FINAL frame must be the whole of it. Both are compared against the
%% same xWant/yWant/zWant arrays section 11 built by hand from traj, at the
%% same exaggeration, so the movie and the still are held to one standard.
%%
%% The frames are reached through the FrameFcn hook, because everything a movie
%% draws is deleted by the time it returns. grabFrame is nested inside this
%% function so it can write into frmSeen:
             nFrmT = 16;
             kMidF = 8;
            mvFile = [tempname '.mp4'];
            frmSeen = cell(nFrmT,1);
             nFigMv = numel(findall(groot,'Type','figure'));
                mv = coorbital.viz.globeMovie(trjCha,vehOne,env, ...
                         struct('File',mvFile,'NFrame',nFrmT, ...
                                'Size',[480 320],'AltScale',altSc, ...
                                'SpinDeg',0,'Texture','plain','Sky','black', ...
                                'PhaseName',{{'boost','glide','descent'}}, ...
                                'FrameFcn',@grabFrame));

%% The file exists, is not empty, and says what it is:
    assert(isfile(mv.file),'globeMovie reported "%s", which is not a file.',mv.file);
             mvInf = dir(mv.file);
    assert(mvInf.bytes > 0,'globeMovie wrote a zero-byte file.');
    assert(strcmp(mv.file,mvFile), ...
        'globeMovie wrote "%s" against the requested "%s".',mv.file,mvFile);

%% The frame count is the count that was asked for, and the frames were really
%% rendered rather than merely counted:
    assert(mv.nFrame == nFrmT, ...
        'globeMovie reports %d frames against the %d requested.',mv.nFrame,nFrmT);
    assert(numel(mv.frameIdx) == nFrmT, ...
        'mv.frameIdx holds %d entries for %d frames.',numel(mv.frameIdx),nFrmT);
    assert(sum(~cellfun(@isempty,frmSeen)) == nFrmT, ...
        'the per-frame hook fired %d times for %d frames.', ...
        sum(~cellfun(@isempty,frmSeen)),nFrmT);

%% MPEG-4 needs both dimensions even, and a frame size that is not is a file
%% that will not open. Pinned, because the cropping that guarantees it is easy
%% to lose in a refactor:
    assert(all(mod(mv.frameSize,2) == 0), ...
        'globeMovie wrote %d x %d frames; MPEG-4 needs both dimensions even.', ...
        mv.frameSize(1),mv.frameSize(2));

%% Which SAMPLE each frame shows, recomputed here from the documented rule --
%% the last sample at or before a uniform walk in TIME -- rather than read back
%% from the function:
              tFrm = linspace(trjCha.t(1),trjCha.t(end),nFrmT)';
            idxWnt = zeros(nFrmT,1);
    for kf = 1:nFrmT
         idxWnt(kf) = find(trjCha.t <= tFrm(kf),1,'last');
    end
    assertSeries(mv.frameIdx,idxWnt,'globeMovie frame sample indices');
    assert(idxWnt(1) == 1 && idxWnt(end) == numel(trjCha.t), ...
        ['the movie must open on the first sample and close on the last; it ' ...
         'opens on %d and closes on %d of %d.'], ...
        idxWnt(1),idxWnt(end),numel(trjCha.t));

%% The intermediate frame must be a PROPER truncation, or the comparison below
%% cannot tell a growing track from a static one:
              kEndM = idxWnt(kMidF);
    assert(kEndM > 1 && kEndM < numel(trjCha.t), ...
        ['frame %d shows sample %d of %d; an intermediate frame that is the ' ...
         'first or the last sample proves nothing about growth.'], ...
        kMidF,kEndM,numel(trjCha.t));

%% THE ASSERTION. Final frame against the whole trajectory, intermediate frame
%% against the trajectory truncated at its own time:
    assertGrownTrack(frmSeen{nFrmT},trjCha.phaseIdx,phList, ...
        numel(trjCha.t),xWant,yWant,zWant,'final');
    assertGrownTrack(frmSeen{kMidF},trjCha.phaseIdx,phList, ...
        kEndM,xWant,yWant,zWant,'intermediate');

%% ...and the two frames really are different pictures. Without this a function
%% that drew the full arc every frame would have to fail the truncation check
%% above, but a function that drew NOTHING would sail through both:
             nPtsM = sum(cellfun(@(m) size(m,1),frmSeen{kMidF}.xyz));
             nPtsF = sum(cellfun(@(m) size(m,1),frmSeen{nFrmT}.xyz));
    assert(nPtsM > 0 && nPtsM < nPtsF, ...
        ['the intermediate frame holds %d track points against the final ' ...
         'frame''s %d; the track is not growing.'],nPtsM,nPtsF);

%% The MARKERS ride a shell 2 per cent above the sphere, because a marker drawn
%% exactly on the surface is bisected by the depth buffer and renders as a
%% half-disc. Two per cent, not the 0.2 that looks sufficient: the marker is a
%% screen-aligned quad spanning of order 115 km of ground, and 0.2 per cent is
%% only 12.8 km of lift, which was measured to change the rendered ring not at
%% all (41 lit pixels at shell 1.000 and 41 at 1.002; 58 and a closed annulus
%% at 1.02). Do not tighten it back. The shell is recomputed here rather than read back, in the same
%% form the function documents -- max(plotted radius, 1.02 rE), radial only --
%% so a shell applied to the wrong quantity, or applied to the TRACK, fails:
              rMkW = max(rPlt,1.02.*c.rE);
             xMkW  = rMkW.*cos(latC).*cos(lonC)./1000;
             yMkW  = rMkW.*cos(latC).*sin(lonC)./1000;
             zMkW  = rMkW.*sin(latC)./1000;

%% The vehicle is marked at the sample the frame shows, and it MOVES:
    assertSeries(frmSeen{kMidF}.veh(:),[xMkW(kEndM);yMkW(kEndM);zMkW(kEndM)], ...
        'globeMovie vehicle marker, intermediate frame (km)');
    assertSeries(frmSeen{nFrmT}.veh(:),[xMkW(end);yMkW(end);zMkW(end)], ...
        'globeMovie vehicle marker, final frame (km)');

%% ...and the shell really is doing something at the end of this trajectory, or
%% the two assertions above would be indistinguishable from no shell at all.
%% The synthetic chain ends at 0.1 km, so the final sample is a case where the
%% marker must be LIFTED off the track it is riding:
             liftM = sqrt((xMkW(end) - xWant(end)).^2 + ...
                          (yMkW(end) - yWant(end)).^2 + ...
                          (zMkW(end) - zWant(end)).^2);
    assert(liftM > 1, ...
        ['the marker shell lifts the final vehicle marker by only %.4f km; ' ...
         'this check cannot tell a shell from no shell.'],liftM);

%% The launch marker sits on the shell too, over the ground the flight began
%% on. It is NOT exaggerated: it marks a ground position:
    assertSeries(frmSeen{1}.lau(:), ...
        1.02.*(c.rE./1000).*[cos(latC(1)).*cos(lonC(1)); ...
                              cos(latC(1)).*sin(lonC(1)); ...
                              sin(latC(1))], ...
        'globeMovie launch marker (km)');

%% ---------------------------------------------------------------------
%% 12b. globeMovie -- the picture reads as THREE phases, not two
%% ---------------------------------------------------------------------
%% A boost-glide chain puts two orders of magnitude more length into the glide
%% than into the terminal descent, and drawn at one width on a whole globe the
%% descent measured FIVE coloured pixels in a 1280 x 720 frame -- three phases
%% in the data and two in the picture. None of that is visible to a test that
%% reads properties, so what is pinned here is the MECHANISM: the width grading,
%% the junction markers and the legend. The frame itself still has to be looked
%% at, and was.
%%
%% Width first. The phase drawn over the least distance must be drawn WIDER
%% than the phase drawn over the most, and the two must not be equal:
              lenP = zeros(1,numel(phList));
    for kp = 1:numel(phList)
               sel = phaseSpan(trjCha.phaseIdx,phList,kp);
          lenP(kp) = sum(sqrt(diff(xWant(sel)).^2 + ...
                              diff(yWant(sel)).^2 + ...
                              diff(zWant(sel)).^2));
    end
        [~,kShort] = min(lenP);
         [~,kLong] = max(lenP);
    assert(lenP(kShort) < 0.5.*lenP(kLong), ...
        ['the shortest synthetic phase is %.1f km against the longest %.1f km; ' ...
         'the width grading cannot be told from a constant width.'], ...
        lenP(kShort),lenP(kLong));
              lwF  = frmSeen{nFrmT}.lw;
              phF  = frmSeen{nFrmT}.ph;
    assert(lwF(phF == phList(kShort)) > lwF(phF == phList(kLong)) + 0.5, ...
        ['the shortest phase is drawn %.2f pt wide against the longest phase''s ' ...
         '%.2f pt; a short leg drawn no bolder than a long one is a leg a ' ...
         'reader cannot see.'], ...
        lwF(phF == phList(kShort)),lwF(phF == phList(kLong)));

%% A junction marker per handoff, in the right place, and -- the part that
%% matters -- NOT drawn before the flight reaches it. Frame 8 of this
%% trajectory has passed the first junction and not the second, so one marker
%% must be placed and one must be off-screen:
              jWnt = zeros(1,numel(phList) - 1);
    for kp = 2:numel(phList)
        jWnt(kp-1) = find(trjCha.phaseIdx == phList(kp-1),1,'last');
    end
    assert(numel(frmSeen{nFrmT}.jctP) == numel(jWnt), ...
        'globeMovie drew %d handoff markers for %d junctions.', ...
        numel(frmSeen{nFrmT}.jctP),numel(jWnt));
    assert(jWnt(1) <= kEndM && jWnt(end) > kEndM, ...
        ['frame %d must sit BETWEEN two junctions for the growth check below ' ...
         'to mean anything; junctions are at samples %s and it shows %d.'], ...
        kMidF,mat2str(jWnt),kEndM);
    for kj = 1:numel(jWnt)
              kRow = find(frmSeen{nFrmT}.jctP == phList(kj+1),1);
        assert(~isempty(kRow), ...
            'no handoff marker carries the index of phase %d.',phList(kj+1));
        assertSeries(frmSeen{nFrmT}.jctXyz(kRow,:)', ...
            [xMkW(jWnt(kj));yMkW(jWnt(kj));zMkW(jWnt(kj))], ...
            sprintf('globeMovie handoff marker %d (km)',kj));
              kRwM = find(frmSeen{kMidF}.jctP == phList(kj+1),1);
               gotM = frmSeen{kMidF}.jctXyz(kRwM,:);
        if jWnt(kj) <= kEndM
            assertSeries(gotM(:), ...
                [xMkW(jWnt(kj));yMkW(jWnt(kj));zMkW(jWnt(kj))], ...
                sprintf('globeMovie handoff marker %d, intermediate frame (km)',kj));
        else
            assert(all(isnan(gotM)), ...
                ['handoff marker %d is drawn at %s in frame %d, which has only ' ...
                 'reached sample %d of the junction''s %d. A handoff drawn ' ...
                 'before the vehicle gets there announces the ending.'], ...
                kj,mat2str(gotM,4),kMidF,kEndM,jWnt(kj));
        end
    end

%% The legend names the phases, in the caller's own words. Without it the
%% colours are a code the frame never explains:
    assert(isequal(reshape(frmSeen{nFrmT}.leg,1,[]),{'boost','glide','descent'}), ...
        'globeMovie''s legend reads {%s}; it must name the phases it was given.', ...
        strjoin(frmSeen{nFrmT}.leg,', '));

%% The readout updates per frame, states its units, and names the phase the
%% vehicle is actually in. A frozen annotation on a moving vehicle is worse
%% than none, because a reader would believe it:
            phNames = {'boost','glide','descent'};
            altKmT  = (trjCha.x(:,1) - c.rE)./1000;
    for kf = [kMidF nFrmT]
              kEndK = idxWnt(kf);
              hudTx = frmSeen{kf}.hud;
        assert(contains(hudTx,sprintf('%.1f s',trjCha.t(kEndK))), ...
            'frame %d readout "%s" does not state its time %.1f s, with the unit.', ...
            kf,hudTx,trjCha.t(kEndK));
        assert(contains(hudTx,sprintf('%.2f km',altKmT(kEndK))), ...
            'frame %d readout "%s" does not state its altitude %.2f km, with the unit.', ...
            kf,hudTx,altKmT(kEndK));
        assertUnit(hudTx,phNames{trjCha.phaseIdx(kEndK)}, ...
            sprintf('frame %d readout',kf));
    end
    assert(~strcmp(frmSeen{kMidF}.hud,frmSeen{nFrmT}.hud), ...
        'the readout is identical in frames %d and %d; it is not updating.', ...
        kMidF,nFrmT);

%% The caption states the exaggeration. A movie whose vertical scale is thirty
%% times its horizontal one and does not say so is the same lie a still figure
%% would be telling, and section 11 holds globe3D to exactly this:
             titMv = frmSeen{nFrmT}.tit;
    assert(contains(titMv,'exaggerat') && contains(titMv,sprintf('%g',altSc)), ...
        ['globeMovie was given AltScale = %g and its caption does not say so. ' ...
         'Caption was "%s".'],altSc,titMv);

%% The phases are drawn in DISTINGUISHABLE colours -- the whole reason the
%% track is split into one line per phase rather than drawn as one polyline:
              colF = frmSeen{nFrmT}.col;
    for ka = 1:numel(colF)
        for kb = ka+1:numel(colF)
            assert(max(abs(colF{ka} - colF{kb})) > 0.1, ...
                'globeMovie drew two phases in the same colour %s.', ...
                mat2str(colF{ka},3));
        end
    end

%% The planet is the SHARED renderer's, at the right radius. A movie that grew
%% its own Earth would drift out of step with coorbital.viz.globe3D:
    assert(frmSeen{nFrmT}.nEarth == 1, ...
        'globeMovie drew %d Earth surfaces; expected exactly 1.', ...
        frmSeen{nFrmT}.nEarth);
    assert(abs(frmSeen{nFrmT}.rEarth - c.rE./1000) < 1e-6, ...
        'the movie''s Earth has radius %.3f km against rE = %.3f km.', ...
        frmSeen{nFrmT}.rEarth,c.rE./1000);

%% No figure survives the call. A movie renders into a figure and must take it
%% away with it, or a batch that writes twenty movies ends with twenty windows:
    assert(numel(findall(groot,'Type','figure')) == nFigMv, ...
        'globeMovie left %d figure(s) open.', ...
        numel(findall(groot,'Type','figure')) - nFigMv);

%% The test cleans up after itself: the repository must stay free of binaries:
    delete(mvFile);
    assert(~isfile(mvFile),'the test movie "%s" could not be deleted.',mvFile);

%% ---------------------------------------------------------------------
%% 12c. globeMovie -- the altitude inset says what the globe cannot
%% ---------------------------------------------------------------------
%% The globe cannot show height: at the sub-camera point the altitude vector
%% points at the lens. The inset is the answer, and it shipped with NOTHING
%% asserted about it -- the render above drew it every frame and looked away.
%% Four mutations lived in it: drawing the whole flight in every frame, putting
%% downrange in metres under a label that says kilometres, exaggerating the
%% profile that calls itself true scale, and a fixed marker size. Each of the
%% checks below kills one.
%%
%% Downrange and altitude are recomputed HERE from traj, in the units the
%% panel's own labels claim, so a unit slip moves one side and not the other:
             dRngW = coorbital.util.greatCircle(latC(1),lonC(1),latC,lonC) ...
                     .*(c.rE./1000);
            altKmW = (trjCha.x(:,1) - c.rE)./1000;
    assert(frmSeen{nFrmT}.nIns == 1, ...
        ['globeMovie drew %d axes tagged insetAxes; Inset defaults to true ' ...
         'and must produce exactly one.'],frmSeen{nFrmT}.nIns);

%% The labels, each naming its quantity AND its unit. "downrange" without
%% "(km)" is the same silent unit error every other panel in this file is
%% checked for, and "true scale" is the claim that distinguishes this panel
%% from the exaggerated arc above it:
    assertUnit(frmSeen{nFrmT}.insLbl{1},'downrange','inset x axis');
    assertUnit(frmSeen{nFrmT}.insLbl{1},'km'       ,'inset x axis');
    assertUnit(frmSeen{nFrmT}.insLbl{2},'altitude' ,'inset y axis');
    assertUnit(frmSeen{nFrmT}.insLbl{2},'km'       ,'inset y axis');
    assertUnit(frmSeen{nFrmT}.insLbl{3},'true scale','inset title');

%% THE PROFILE IS THE TRAJECTORY, phase by phase, and it GROWS. The final
%% frame carries the whole flight and the intermediate frame carries it
%% truncated at its own sample -- the same contract the globe track is held to
%% in section 12, applied to the panel that was never checked:
    for kp = 1:numel(phList)
        for kf = [kMidF nFrmT]
              kCut = idxWnt(kf);
              kRow = find(frmSeen{kf}.insPh == phList(kp),1);
            assert(~isempty(kRow), ...
                'the inset has no line tagged for phase %d.',phList(kp));
               sel = phaseSpan(trjCha.phaseIdx,phList,kp);
               sel = sel(sel <= kCut);
               got = frmSeen{kf}.insXy{kRow};
            assert(size(got,1) == numel(sel), ...
                ['inset frame %d, phase %d: %d points drawn against the %d ' ...
                 'samples the flight has reached. The panel is not growing ' ...
                 'with the flight.'],kf,phList(kp),size(got,1),numel(sel));
            if isempty(sel)
                continue;
            end
            assertSeries(got(:,1),dRngW(sel), ...
                sprintf('inset frame %d, phase %d downrange (km)',kf,phList(kp)));
            assertSeries(got(:,2),altKmW(sel), ...
                sprintf('inset frame %d, phase %d altitude (km)',kf,phList(kp)));
        end
    end

%% ...and the two frames are different pictures, so a panel that drew nothing
%% cannot pass the truncation check by drawing nothing twice:
             nInsM = sum(cellfun(@(m) size(m,1),frmSeen{kMidF}.insXy));
             nInsF = sum(cellfun(@(m) size(m,1),frmSeen{nFrmT}.insXy));
    assert(nInsM > 0 && nInsM < nInsF, ...
        ['the inset holds %d points in the intermediate frame against %d in ' ...
         'the final one; the profile is not growing.'],nInsM,nInsF);

%% TRUE SCALE, AND THE CHECK IS DISCRIMINATING. The globe arc above is drawn
%% at AltScale; the panel must not be. The assertion above compares against
%% altKmW, so an exaggerated panel fails it -- but only because the two differ,
%% which is stated here rather than assumed:
    assert(altSc ~= 1 && max(altKmW) > 0, ...
        ['the movie was rendered at AltScale = %g on a flight peaking at ' ...
         '%.3f km; the true-scale check cannot tell the panel from the arc.'], ...
        altSc,max(altKmW));

%% The moving point rides the profile at the sample the frame shows:
    assertSeries(frmSeen{kMidF}.insVeh(:),[dRngW(kEndM);altKmW(kEndM)], ...
        'inset vehicle point, intermediate frame (km)');
    assertSeries(frmSeen{nFrmT}.insVeh(:),[dRngW(end);altKmW(end)], ...
        'inset vehicle point, final frame (km)');

%% The panel's colours are the globe's colours. A profile in colours of its own
%% would make a reader match legs by shape:
    for kp = 1:numel(phList)
              kGlb = find(frmSeen{nFrmT}.ph    == phList(kp),1);
              kIns = find(frmSeen{nFrmT}.insPh == phList(kp),1);
        assert(max(abs(frmSeen{nFrmT}.col{kGlb} - ...
                       frmSeen{nFrmT}.insCl{kIns})) < 1e-9, ...
            ['phase %d is drawn %s on the globe and %s in the inset; the two ' ...
             'must be one colour or the panel cannot be read against the arc.'], ...
            phList(kp),mat2str(frmSeen{nFrmT}.col{kGlb},3), ...
            mat2str(frmSeen{nFrmT}.insCl{kIns},3));
    end

%% The limits are set from the WHOLE flight, so the curve grows into a frame
%% that never rescales, and they are FLOORED so a flight with no vertical
%% extent still produces a legal axis -- section 12f flies one:
             insLm = frmSeen{nFrmT}.insLim;
    assert(insLm(2) >= max(dRngW) && insLm(4) >= max(altKmW), ...
        ['the inset is limited to %.3f km downrange by %.3f km altitude on a ' ...
         'flight that reaches %.3f by %.3f; the curve would run off the panel.'], ...
        insLm(2),insLm(4),max(dRngW),max(altKmW));
    assert(insLm(2) > insLm(1) && insLm(4) > insLm(3), ...
        'the inset limits [%g %g] x [%g %g] are not increasing.', ...
        insLm(1),insLm(2),insLm(3),insLm(4));

%% THE PANEL IS TRANSPARENT. An opaque panel over a globe hides whatever the
%% flight put behind it, and on the shipped run_target case it hid part of the
%% glide -- the track descended behind its left edge and re-emerged above its
%% top. Pinned as a property because no property check can see the picture:
    assert(strcmp(frmSeen{nFrmT}.insCol,'none'), ...
        ['the inset axes Color is %s; an opaque panel over the planet hides ' ...
         'the track behind it.'],mat2str(frmSeen{nFrmT}.insCol));

%% ...and its x label is ON THE PAGE. The label under the panel is an absolute
%% stack of tick text roughly 26 px tall, so a panel whose origin is a fixed
%% FRACTION of the figure loses the label at every small frame size -- which
%% is both the 480 x 320 rendered here and the 640 x 360 the self-demo uses.
%% Measured in pixels, because pixels are what the label occupies:
             yInPx = frmSeen{nFrmT}.insPos(2).*320;
    assert(yInPx > 25.99, ...
        ['the inset sits %.1f px off the bottom of a 320 px frame; its x ' ...
         'label needs 26 and is clipped.'],yInPx);

%% ---------------------------------------------------------------------
%% 12d. globeMovie -- the camera aims where it says it aims
%% ---------------------------------------------------------------------
%% Everything the arc-frame construction does ends in view(az,el), and nothing
%% read it. The camera aim, the plane normal, the along-arc offset and the
%% out-of-plane tilt were therefore all unpinned at once: negating the tilt,
%% or swapping the endpoints the midpoint is built from, changed every frame of
%% every movie and failed no assertion.
%%
%% The render above used SpinDeg = 0, so its camera is the same in every frame
%% and is exactly the arc-frame direction. Recomputed here from the two
%% endpoints and the DOCUMENTED defaults, restated in the test's own hand:
            offDef = 0;
            tltDef = 32;
               uvT = @(a,b) [cos(a).*cos(b); cos(a).*sin(b); sin(a)];
               u1T = uvT(latC(1)  ,lonC(1));
               u2T = uvT(latC(end),lonC(end));
               nPT = unitOf3(cross(u1T,u2T));
               mMT = unitOf3(u1T + u2T);
               tTT = unitOf3(cross(nPT,mMT));
               cDT = cosd(tltDef).*(cosd(offDef).*mMT + sind(offDef).*tTT) ...
                     + sind(tltDef).*nPT;
               cDT = unitOf3(cDT);
    for kf = [1 kMidF nFrmT]
        assertCamera(frmSeen{kf}.view,cDT, ...
            sprintf('globeMovie camera, frame %d',kf));
    end

%% ...and the assertion is DISCRIMINATING, which is the half that was missing.
%% Each of the three ways this construction has been got wrong must move the
%% aim by far more than the tolerance above, or the check is passing on a
%% coincidence rather than on correctness:
               cNeg = unitOf3(cosd(tltDef).*mMT - sind(tltDef).*nPT);
               cSwp = unitOf3(cosd(tltDef).*mMT ...
                              + sind(tltDef).*unitOf3(cross(u2T,u1T)));
               kMed = ceil(numel(trjCha.t)./2);
               cMed = unitOf3(cosd(tltDef).*uvT(latC(kMed),lonC(kMed)) ...
                              + sind(tltDef).*nPT);
    assert(angBetween(cDT,cNeg) > 5, ...
        ['a negated ViewTiltDeg moves the aim only %.3f deg; the camera check ' ...
         'cannot see the sign of the tilt.'],angBetween(cDT,cNeg));
    assert(angBetween(cDT,cSwp) > 5, ...
        ['swapping the two endpoints moves the aim only %.3f deg; the camera ' ...
         'check cannot see which way the arc frame was built.'], ...
        angBetween(cDT,cSwp));
%% The median-sample threshold is ONE degree and not five, and the reason is
%% worth stating: this synthetic trajectory is uniform in index, so its median
%% SAMPLE sits close to its arc midpoint and the two aims differ by a measured
%% 4.85 deg. The defect this replaced was much worse than that in the field --
%% ode45 clusters samples where the dynamics are fast, and on a 97 deg
%% intercontinental arc the median sample hid the launch point entirely -- so
%% a synthetic case cannot reproduce its full size. It does not need to:
%% assertCamera's budget is 1e-6 deg, so even 4.85 is a margin of five million:
    assert(angBetween(cDT,cMed) > 1, ...
        ['aiming at the MEDIAN SAMPLE instead of the arc midpoint moves the ' ...
         'aim only %.3f deg on this trajectory; the check cannot see the ' ...
         'defect it exists for.'],angBetween(cDT,cMed));

%% ---------------------------------------------------------------------
%% 12e. globeMovie -- the midpoint is the great-circle midpoint, and the
%%      spin is CENTRED on it
%% ---------------------------------------------------------------------
%% Section 12d recomputes the function's own arithmetic. This one does not: it
%% flies the camera flat, at ViewTiltDeg = 0 and ViewOffsetDeg = 0, where the
%% aim must be the great-circle midpoint of the endpoints and nothing else, and
%% compares it against the standard midpoint formula -- a different expression
%% of the same geometry, so an error common to both sides is not possible.
%%
%% The frame is rendered TALL on purpose. The inset marker is sized from the
%% frame like everything else in the file, and at 320 or 720 px of height the
%% rule and the 5 pt constant it replaced give the same answer; only above 720
%% do they part company, so only a tall frame can tell them apart:
            nFrmV = 4;
            spinV = 40;
            sizeV = [320 1080];
          viewSeen = zeros(nFrmV,2);
         insMkSeen = NaN;
            mvFilV = [tempname '.mp4'];
    coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',mvFilV,'NFrame',nFrmV,'Size',sizeV, ...
               'AltScale',altSc,'SpinDeg',spinV, ...
               'ViewOffsetDeg',0,'ViewTiltDeg',0, ...
               'Texture','plain','Sky','black','FrameFcn',@grabView));
    delete(mvFilV);

%% The great-circle midpoint, from the formula, in the test's own hand:
      [latMd,lonMd] = gcMidpoint(latC(1),lonC(1),latC(end),lonC(end));
               cMid = uvT(latMd,lonMd);

%% The sweep is CENTRED on it: the first and the last frame sit half the sweep
%% either side, so their mean is the aim and their difference is SpinDeg. A
%% one-sided sweep -- the defect this replaced -- puts the aim at frame one and
%% walks a long arc off the limb by the end:
               azMd = rad2deg(atan2(cMid(2),cMid(1))) + 90;
               elMd = rad2deg(asin(cMid(3)));
    assertAngle(viewSeen(1    ,1),azMd - spinV./2,1e-6, ...
        'first-frame azimuth against the midpoint less half the sweep (deg)');
    assertAngle(viewSeen(nFrmV,1),azMd + spinV./2,1e-6, ...
        'last-frame azimuth against the midpoint plus half the sweep (deg)');
    for kf = 1:nFrmV
        assertAngle(viewSeen(kf,2),elMd,1e-6, ...
            sprintf('frame %d elevation at zero tilt against the midpoint (deg)',kf));
    end

%% ...and the midpoint is not the launch point, the impact point or the median
%% sample, so the three assertions above are statements about the MIDPOINT and
%% not about whichever sample happens to sit near it:
    for kEnds = [1 kMed numel(trjCha.t)]
        assert(angBetween(cMid,uvT(latC(kEnds),lonC(kEnds))) > 5, ...
            ['the arc midpoint is only %.3f deg from sample %d; this ' ...
             'trajectory cannot tell the midpoint from that sample.'], ...
            angBetween(cMid,uvT(latC(kEnds),lonC(kEnds))),kEnds);
    end

%% The inset marker grew with the frame. mkIns is half the globe marker,
%% floored at 5 pt, so a 1080 px frame must carry more than the 5 pt a 320 px
%% frame does; a hard-coded 5 gives the same speck at every size:
    assert(insMkSeen > frmSeen{nFrmT}.insMkS, ...
        ['the inset marker is %.1f pt at %d px of frame height and %.1f pt at ' ...
         '320; a marker that does not grow with the frame is a speck at ' ...
         '1080.'],insMkSeen,sizeV(2),frmSeen{nFrmT}.insMkS);

%% ---------------------------------------------------------------------
%% 12f. globeMovie -- the two degenerate arcs, and a flight with no height
%% ---------------------------------------------------------------------
%% COINCIDENT and ANTIPODAL endpoints both collapse the trajectory plane, and
%% both used to be handled by one guard that was wrong twice over. The fallback
%% normal was the polar axis whatever the aim was, so the frame was NOT
%% orthonormal -- m . n came out at 0.5 on a 30 deg case -- and ViewTiltDeg was
%% then not the angle out of anything. And the MIDPOINT was left unguarded: for
%% antipodal endpoints u1 + u2 is pure roundoff, measured at 2.0e-16, and
%% normalising it aimed the camera in a direction decided by floating point.
%%
%% Both are checked by ONE property that needs no knowledge of the fallback:
%% the camera must sit exactly ViewTiltDeg away from the point it is aimed at.
%% That is what the option means, it holds in every branch of an orthonormal
%% frame, and it fails for both defects -- the old polar fallback puts the
%% camera 25.8 deg from a 20 N aim point when asked for 32, and the roundoff
%% midpoint puts it 90 deg from anywhere in particular:
            tltDeg = 32;
    for kCase = 1:2
        if kCase == 1
              trDgn = degenTraj(c.rE,'coincident');
              aimWt = uvT(trDgn.x(1,3),trDgn.x(1,2));
              whatT = 'coincident endpoints';
        else
              trDgn = degenTraj(c.rE,'antipodal');
              aimWt = uvT(trDgn.x(1,3),trDgn.x(1,2));
              whatT = 'antipodal endpoints';
        end
          viewSeen = zeros(2,2);
            mvFilD = [tempname '.mp4'];
            nFigDg = numel(findall(groot,'Type','figure'));
        coorbital.viz.globeMovie(trDgn,vehOne,env, ...
            struct('File',mvFilD,'NFrame',2,'Size',[320 240], ...
                   'SpinDeg',0,'ViewOffsetDeg',0,'ViewTiltDeg',tltDeg, ...
                   'Texture','plain','Sky','black','FrameFcn',@grabView));
        delete(mvFilD);
        assert(numel(findall(groot,'Type','figure')) == nFigDg, ...
            'globeMovie leaked a figure on %s.',whatT);
               cGot = camDirOf(viewSeen(1,:));
        assert(all(isfinite(cGot)), ...
            'the camera came back %s on %s.',mat2str(viewSeen(1,:),4),whatT);
               angT = angBetween(cGot,aimWt);
        assert(abs(angT - tltDeg) < 1e-6, ...
            ['the camera sits %.6f deg from the aim point on %s against the ' ...
             '%g deg tilt it was given. Either the fallback frame is not ' ...
             'orthonormal, so the tilt is not the angle it claims, or the aim ' ...
             'itself came out of roundoff.'],angT,whatT,tltDeg);
    end

%% ...and a flight with NO VERTICAL EXTENT AT ALL must still render with the
%% inset on. The panel takes its limits from the flight, and a zero span is an
%% axis limit MATLAB refuses; that error used to be thrown with the video
%% writer already open, so the run left a figure on screen and a part-written
%% MP4 on disk. With Inset false the very same trajectory rendered perfectly,
%% which is what made it an inset defect rather than a trajectory one:
            trFlat = degenTraj(c.rE,'flat');
            mvFilF = [tempname '.mp4'];
            nFigFl = numel(findall(groot,'Type','figure'));
              mvFl = coorbital.viz.globeMovie(trFlat,vehOne,env, ...
                         struct('File',mvFilF,'NFrame',4,'Size',[320 240], ...
                                'Texture','plain','Sky','black'));
    assert(isfile(mvFl.file),'the flat-trajectory movie was not written.');
             infFl = dir(mvFl.file);
    assert(infFl.bytes > 0,'the flat-trajectory movie is empty.');
    delete(mvFilF);
    assert(numel(findall(groot,'Type','figure')) == nFigFl, ...
        'globeMovie leaked a figure on a flat trajectory.');

%% ---------------------------------------------------------------------
%% 13. globeMovie -- what it must REFUSE
%% ---------------------------------------------------------------------
%% All four are caught BEFORE a frame is rendered, which is the point: a movie
%% that renders for ten minutes and then cannot be written has wasted ten
%% minutes, and a one-frame movie is a still with a codec:
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.avi'],'NFrame',4)),'an output that is not .mp4');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',fullfile(tempdir,'no_such_folder_here','m.mp4'),'NFrame',4)), ...
        'an output folder that does not exist');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',1)),'a one-frame movie');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'AltScale',-3)), ...
        'a negative altitude exaggeration');

%% ...and a misspelt Texture or Sky, because silently downgrading one would be
%% indistinguishable from pumpkyn not being installed, which is the one thing
%% the function's degradation note promises a reader can tell apart:
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'Texture','atuo')), ...
        'a misspelt Texture');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'Sky','starz')), ...
        'a misspelt Sky');

%% ...and the three camera-and-inset options, which shipped with no check at
%% all while every other option in the file had one. The NaN cases are the
%% reason this matters: view() takes a NaN angle without complaint and leaves
%% the camera on MATLAB's default, so a NaN tilt produced a movie shot from
%% the wrong place with nothing anywhere saying the option had been thrown
%% away. A refusal is the only way a caller finds out:
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'ViewTiltDeg',NaN)), ...
        'a NaN ViewTiltDeg');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'ViewOffsetDeg',NaN)), ...
        'a NaN ViewOffsetDeg');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'ViewTiltDeg',[10 20])), ...
        'a non-scalar ViewTiltDeg');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'Inset','yes')), ...
        'a character Inset');
    assertThrows(@() coorbital.viz.globeMovie(trjCha,vehOne,env, ...
        struct('File',[tempname '.mp4'],'NFrame',4,'Inset',[true true])), ...
        'a non-scalar Inset');

%% ---------------------------------------------------------------------
%% 14. 'Parent' composes instead of creating
%% ---------------------------------------------------------------------
%% The whole point of the option: handed somewhere to draw, none of these may
%% open a figure of their own. Counted, not assumed:
             hHost = figure('Visible','off');
               axG = axes('Parent',hHost);
           nBefore = numel(findall(groot,'Type','figure'));
               fGT = coorbital.viz.groundTrack(trjCha,vehOne,env, ...
                         struct('Parent',axG));
    assert(numel(findall(groot,'Type','figure')) == nBefore, ...
        'groundTrack opened a figure although it was handed an axes.');
    assert(isequal(fGT,hHost), ...
        'groundTrack returned a figure other than the one its Parent lives in.');
    assert(~isempty(lineByTag(axG,'groundTrack')), ...
        'groundTrack drew nothing into the axes it was given.');

%% The globe goes into its own axes in the SAME figure -- a lat/lon grid and a
%% sphere sharing one set of limits would be unreadable, and the point being
%% made is only that no new figure appears:
              axG2 = axes('Parent',hHost);
               fG3 = coorbital.viz.globe3D(trjCha,vehOne,env, ...
                         struct('Parent',axG2));
    assert(numel(findall(groot,'Type','figure')) == nBefore, ...
        'globe3D opened a figure although it was handed an axes.');
    assert(isequal(fG3,hHost), ...
        'globe3D returned a figure other than the one its Parent lives in.');
    assert(~isempty(lineByTag(axG2,'globeTrack')), ...
        'globe3D drew nothing into the axes it was given.');

%% profilePlot needs one axes per panel, and the panel count is the caller's
%% own choice, so it is handed exactly as many as it asked for:
            hHost5 = figure('Visible','off');
               ax5 = gobjects(1,5);
    for ka = 1:5
           ax5(ka) = axes('Parent',hHost5, ...
                          'OuterPosition',[0 (5-ka)./5 1 1./5]);
    end
          nBefore5 = numel(findall(groot,'Type','figure'));
               fPP = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Parent',ax5,'VehPhase',{vehPhase}));
    assert(numel(findall(groot,'Type','figure')) == nBefore5, ...
        'profilePlot opened a figure although it was handed five axes.');
    assert(isequal(fPP,hHost5), ...
        'profilePlot returned a figure other than the one its Parent lives in.');
    for ka = 1:5
        assert(isscalar(lineByTag(ax5(ka),'profile')), ...
            'profilePlot did not draw into supplied axes %d.',ka);
    end

%% ...and the wrong number of axes must be refused rather than half-filled.
%% Five axes are the right number for the DEFAULT channels and the wrong
%% number the moment a sixth channel or an Extra panel is asked for, which is
%% the mistake this refusal exists to catch:
    assertThrows(@() coorbital.viz.profilePlot(trjCha,vehOne,env, ...
        struct('Parent',ax5(1:3))),'three axes for a five-panel figure');
    assertThrows(@() coorbital.viz.profilePlot(trjCha,vehOne,env, ...
        struct('Parent',ax5,'VehPhase',{vehPhase},'Channels',{allChn})), ...
        'five axes for a seven-channel figure');

    close(hHost);
    close(hHost5);

%% ---------------------------------------------------------------------
%% 15. No leaked figures
%% ---------------------------------------------------------------------
             nFig1 = numel(findall(groot,'Type','figure'));
    assert(nFig1 == nFig0, ...
        ['test_viz leaked %d figure(s); the suite must end with the %d it ' ...
         'started with.'],nFig1 - nFig0,nFig0);

    function grabFrame(hAxF,kF)
%% Purpose:
%
%  The FrameFcn hook handed to coorbital.viz.globeMovie: read one rendered
%  frame's drawn objects out of the axes and keep them. NESTED inside test_viz
%  rather than written as a subfunction, because it has to write into frmSeen,
%  and an anonymous function cannot assign.
%
%  This is the only way an intermediate frame can be inspected at all: a movie
%  deletes its figure before it returns, so by the time the call is over there
%  is nothing left to read.
%
%% Inputs:
%
%  hAxF             [1 x 1] axes                The movie's axes, mid-render
%
%  kF               [1 x 1]                     Frame number (-)
%
%% Outputs:
%
%  none                                         Writes frmSeen{kF}
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hTr = findobj(hAxF,'Type','line','Tag','globeTrack');
               got = struct();
           got.ph  = zeros(1,numel(hTr));
           got.xyz = cell(1,numel(hTr));
           got.col = cell(1,numel(hTr));
           got.lw  = zeros(1,numel(hTr));
        for kr = 1:numel(hTr)
        got.ph(kr) = hTr(kr).UserData;
       got.xyz{kr} = [hTr(kr).XData(:),hTr(kr).YData(:),hTr(kr).ZData(:)];
       got.col{kr} = hTr(kr).Color;
        got.lw(kr) = hTr(kr).LineWidth;
        end
               hVe = findobj(hAxF,'Type','line','Tag','vehicleMarker');
           got.veh = [hVe.XData,hVe.YData,hVe.ZData];
               hLa = findobj(hAxF,'Type','line','Tag','launchMarker');
           got.lau = [hLa.XData,hLa.YData,hLa.ZData];

%% The junction markers, phase index and position together, so a marker drawn
%% at the wrong junction is a different thing from one drawn at the right one:
               hHd = findobj(hAxF,'Type','line','Tag','handoffMarker');
          got.jctP = zeros(1,numel(hHd));
        got.jctXyz = zeros(numel(hHd),3);
        for kr = 1:numel(hHd)
      got.jctP(kr) = hHd(kr).UserData;
   got.jctXyz(kr,:) = [hHd(kr).XData,hHd(kr).YData,hHd(kr).ZData];
        end

%% The caption and the readout are FIGURE annotations, not axes children, so
%% they are looked up in the figure. findall rather than findobj, because an
%% annotation lives in a hidden overlay axes that findobj will not descend
%% into:
               hFg = ancestor(hAxF,'figure');
               hHu = findall(hFg,'Tag','hudText');
           got.hud = labelText(hHu.String);
               hTi = findall(hFg,'Tag','titleText');
           got.tit = labelText(hTi.String);
               hLg = findobj(hFg,'Type','legend');
           got.leg = {};
        if isscalar(hLg)
           got.leg = cellstr(string(hLg.String));
        end
               hEa = findobj(hAxF,'Tag','earthSurface');
        got.nEarth = numel(hEa);
        got.rEarth = NaN;
        if isscalar(hEa)
        got.rEarth = max(hEa.XData(:));
        end

%% THE CAMERA. view(az,el) is the entire visible output of the arc-frame
%% construction -- the midpoint, the plane normal, the offset and the tilt all
%% end in these two numbers and nowhere else -- so a frame that does not record
%% them cannot see any of that code work:
          got.view = hAxF.View;

%% THE ALTITUDE INSET, everything about it that can be compared against the
%% trajectory: the drawn series per phase, the moving point, the limits, the
%% labels, the marker size, the panel geometry and whether the panel is
%% transparent:
          got.nIns = 0;
         got.insXy = {};
         got.insPh = [];
         got.insCl = {};
        got.insVeh = [];
        got.insLim = [];
        got.insPos = [];
        got.insMkS = NaN;
        got.insCol = '';
        got.insLbl = {'','',''};
               hIn = findobj(hFg,'Type','axes','Tag','insetAxes');
          got.nIns = numel(hIn);
        if isscalar(hIn)
               hIl = findobj(hIn,'Type','line','Tag','insetTrack');
         got.insPh = zeros(1,numel(hIl));
         got.insXy = cell(1,numel(hIl));
         got.insCl = cell(1,numel(hIl));
            for kr = 1:numel(hIl)
     got.insPh(kr) = hIl(kr).UserData;
    got.insXy{kr}  = [hIl(kr).XData(:),hIl(kr).YData(:)];
    got.insCl{kr}  = hIl(kr).Color;
            end
               hIv = findobj(hIn,'Type','line','Tag','insetVehicle');
        got.insVeh = [hIv.XData,hIv.YData];
        got.insMkS = hIv.MarkerSize;
        got.insLim = [hIn.XLim,hIn.YLim];
        got.insPos = hIn.Position;
        got.insCol = hIn.Color;
        got.insLbl = {labelText(hIn.XLabel.String), ...
                      labelText(hIn.YLabel.String), ...
                      labelText(hIn.Title.String)};
        end
      frmSeen{kF}  = got;
    end

    function grabView(hAxF,kF)
%% Purpose:
%
%  A lighter FrameFcn than grabFrame: record only the camera angles of every
%  frame, and the inset marker size. Used by the renders below that exist to
%  check the camera and the frame-size scaling rather than the drawn data, so
%  they cost one property read per frame instead of a whole scene walk.
%
%% Inputs:
%
%  hAxF             [1 x 1] axes                The movie's axes, mid-render
%
%  kF               [1 x 1]                     Frame number (-)
%
%% Outputs:
%
%  none                                         Writes viewSeen and insMkSeen
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    viewSeen(kF,:) = hAxF.View;
               hFv = ancestor(hAxF,'figure');
               hIv = findobj(hFv,'Type','axes','Tag','insetAxes');
        if isscalar(hIv)
               hMv = findobj(hIv,'Type','line','Tag','insetVehicle');
         insMkSeen = hMv.MarkerSize;
        end
    end
end

function traj = chainTraj(rE)
%% Purpose:
%
%  A synthetic three-phase, seven-state trajectory: a climb, a long shallow
%  glide and a dive, with the mass falling through phase 1 and constant
%  afterwards. It is not a solution of anything, and it does not need to be --
%  what a plotting test needs is a state history whose every column is
%  distinct and whose phase structure is known exactly.
%
%  The mass column carries a real staging DROP at the phase 1/2 boundary, so
%  the 'mass' channel is exercised on the discontinuity it exists to show and
%  not on a smooth curve.
%
%% Inputs:
%
%  rE               [1 x 1]                     Reference sphere radius (m)
%
%% Outputs:
%
%  traj             Struct                      Fields t [N x 1] (s),
%                                               x [N x 7], u [N x 2] (rad),
%                                               phaseIdx [N x 1], junction
%                                               [2 x 1] struct
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                n1 = 20;
                n2 = 40;
                n3 = 15;
                nS = n1 + n2 + n3;
            traj.t = linspace(0,2100,nS)';

%% Altitude: up to 47 km through boost, a decaying skip through the glide,
%% then down to the deck:
               hKm = [linspace(0.2,47,n1), ...
                      40 + 8.*exp(-linspace(0,2.2,n2)).*cos(linspace(0,7,n2)), ...
                      linspace(15,0.1,n3)]';
               spd = [linspace(10,4600,n1), ...
                      linspace(4600,900,n2), ...
                      linspace(900,400,n3)]';
               lon = deg2rad(linspace(-155,-77.6,nS))';
               lat = deg2rad(linspace(20,34,nS))';
               gam = deg2rad([linspace(89,2,n1), ...
                              2.*cos(linspace(0,7,n2)), ...
                              linspace(-5,-32,n3)])';
               psi = deg2rad(linspace(60,80,nS))';
               mas = [linspace(9400,4300,n1), 900.*ones(1,n2 + n3)]';

            traj.x = [rE + 1000.*hKm, lon, lat, spd, gam, psi, mas];
            traj.u = zeros(nS,2);
     traj.phaseIdx = [ones(n1,1); 2.*ones(n2,1); 3.*ones(n3,1)];
     traj.junction = repmat(struct('t',[],'x',[]),2,1);
    for kj = 1:2
        traj.junction(kj).t = traj.t(find(traj.phaseIdx == kj,1,'last'));
        traj.junction(kj).x = traj.x(find(traj.phaseIdx == kj,1,'last'),:)';
    end
end

function traj = singleTraj(rE)
%% Purpose:
%
%  A synthetic one-phase, six-state trajectory -- the shape HGV/run_glide
%  returns. Exists so the mass-column, mass-channel and per-phase-vehicle
%  branches of coorbital.viz.profilePlot are exercised on the side where none
%  of them is present.
%
%% Inputs:
%
%  rE               [1 x 1]                     Reference sphere radius (m)
%
%% Outputs:
%
%  traj             Struct                      Fields t [N x 1] (s),
%                                               x [N x 6], u [N x 2] (rad),
%                                               phaseIdx [N x 1] all ones,
%                                               junction [0 x 1] struct
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = 50;
            traj.t = linspace(0,2073,nS)';
               hKm = linspace(60,5,nS)';
               spd = linspace(6000,321,nS)';
               lon = deg2rad(linspace(0,62.76,nS))';
               lat = zeros(nS,1);
               gam = deg2rad(linspace(-1,-12.6,nS))';
               psi = deg2rad(90).*ones(nS,1);
            traj.x = [rE + 1000.*hKm, lon, lat, spd, gam, psi];
            traj.u = zeros(nS,2);
     traj.phaseIdx = ones(nS,1);
     traj.junction = repmat(struct('t',[],'x',[]),0,1);
end

function chn = profileChannels(traj,vehPhase,c)
%% Purpose:
%
%  Recompute every channel coorbital.viz.profilePlot can draw, in the units it
%  draws them in, from the trajectory and the per-phase vehicles. The
%  duplication of the plotting routine's arithmetic is the whole mechanism:
%  the two are compared, so a wrong column or a doubled conversion moves one
%  and not the other.
%
%% Inputs:
%
%  traj             Struct                      Trajectory to recompute from
%
%  vehPhase         Cell [1 x P]                One vehicle per phase
%
%  c                Struct                      Constants from missileConst;
%                                               uses rE (m) and g0 (m/s^2)
%
%% Outputs:
%
%  chn              Struct                      One [N x 1] field per channel:
%                                               altitude (km), speed (m/s),
%                                               mach (-), q (kPa), nAero (g),
%                                               gamma (deg), and mass (kg)
%                                               when the state carries one
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = numel(traj.t);
                hM = traj.x(:,1) - c.rE;
                 V = traj.x(:,4);
    [rho,~,~,aSnd] = coorbital.atmos.expAtmos(hM);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;
      chn.altitude = hM./1000;
         chn.speed = V;
          chn.mach = mach;
             chn.q = qbar./1000;
         chn.gamma = rad2deg(traj.x(:,5));
            phList = unique(traj.phaseIdx(:)).';

    if size(traj.x,2) >= 7
          chn.mass = traj.x(:,7);
             massV = traj.x(:,7);
    else
             massV = zeros(nS,1);
        for kp = 1:numel(phList)
               sel = traj.phaseIdx == phList(kp);
        massV(sel) = vehPhase{kp}.mass;
        end
    end

             aLift = zeros(nS,1);
             aDrag = zeros(nS,1);
    for ks = 1:nS
               kPh = find(phList == traj.phaseIdx(ks),1);
              vehK = vehPhase{kPh};
         [CLk,CDk] = coorbital.aero.constLD(traj.u(ks,1),mach(ks),vehK);
         aLift(ks) = qbar(ks).*vehK.Sref.*CLk./massV(ks);
         aDrag(ks) = qbar(ks).*vehK.Sref.*CDk./massV(ks);
    end
         chn.nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;
end

function txt = unitOf(name)
%% Purpose:
%
%  The unit token that must appear in a given channel's axis label. Kept as a
%  lookup rather than derived from the drawn label, because a token read back
%  out of the figure would agree with itself no matter what it said.
%
%% Inputs:
%
%  name             Char [1 x n]                Channel name
%
%% Outputs:
%
%  txt              Char [1 x m]                Required substring of the label
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             known = {'altitude','speed','mach','q','nAero','mass','gamma'};
             units = {'km'      ,'m/s'  ,'(-)' ,'kPa','(g)' ,'(kg)','(deg)'};
               idx = find(strcmp(name,known),1);
    assert(~isempty(idx),'no unit recorded for channel "%s".',name);
               txt = units{idx};
end

function assertPanels(hAx,names,chn,tWant,what)
%% Purpose:
%
%  Assert that a run of profile panels carries, in order, the channels named,
%  each on the trajectory's own time base. This is the assertion that makes
%  the file a test rather than a smoke check.
%
%% Inputs:
%
%  hAx              [1 x M] axes                Panels, in drawing order
%
%  names            Cellstr [1 x M]             Channel expected in each panel
%
%  chn              Struct                      Recomputed channels
%
%  tWant            [N x 1]                     Expected time base (s)
%
%  what             Char [1 x n]                Case name, for the message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(numel(hAx) == numel(names), ...
        '%s: %d panels against %d expected channels.',what,numel(hAx),numel(names));
    for ka = 1:numel(names)
               hLn = lineByTag(hAx(ka),'profile');
        assert(isscalar(hLn), ...
            '%s: panel %d holds %d primary lines; expected 1.',what,ka,numel(hLn));
        assertSeries(hLn.XData(:),tWant, ...
            sprintf('%s panel %d time (s)',what,ka));
        assertSeries(hLn.YData(:),chn.(names{ka}), ...
            sprintf('%s panel %d, channel ''%s''',what,ka,names{ka}));
    end
end

function assertDistinct(names,chn)
%% Purpose:
%
%  Assert that the named channels are pairwise different series. Without this,
%  a function that drew one channel into every panel could satisfy every other
%  check in the file the moment an expectation was ever built from drawn data
%  rather than from the trajectory.
%
%% Inputs:
%
%  names            Cellstr [1 x M]             Channel names to compare
%
%  chn              Struct                      Recomputed channels
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    for ka = 1:numel(names)
        for kb = ka+1:numel(names)
               dif = max(abs(chn.(names{ka}) - chn.(names{kb})));
            assert(dif > 1e-6, ...
                'channels ''%s'' and ''%s'' carry the same series; the checks cannot discriminate.', ...
                names{ka},names{kb});
        end
    end
end

function sel = phaseSpan(phIdx,phList,kp)
%% Purpose:
%
%  The sample indices one phase segment is drawn over: that phase's own
%  samples, with the preceding phase's LAST sample prepended. That is the
%  contract coorbital.viz.groundTrack and coorbital.viz.globe3D both document
%  -- phaseRun labels a junction sample with the outgoing phase, so drawing a
%  phase over its own samples alone leaves a one-segment hole at every
%  boundary. Reconstructed here from phaseIdx alone, so a segment drawn over
%  the wrong index range fails even though the union still looks right.
%
%% Inputs:
%
%  phIdx            [N x 1]                     Phase index per sample
%
%  phList           [1 x P]                     Phase indices present, sorted
%
%  kp               [1 x 1]                     Which entry of phList
%
%% Outputs:
%
%  sel              [M x 1]                     Sample indices, ascending
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               sel = find(phIdx == phList(kp));
    if kp > 1
               sel = [find(phIdx == phList(kp-1),1,'last'); sel];
    end
end

function assertGrownTrack(frm,phIdx,phList,kEnd,xWant,yWant,zWant,what)
%% Purpose:
%
%  Assert that ONE captured frame of coorbital.viz.globeMovie holds the
%  trajectory truncated at sample kEnd -- every phase's own samples up to that
%  point and not one beyond it, in the same phase-by-phase segmentation
%  coorbital.viz.globe3D uses.
%
%  This is what separates a movie of a trajectory developing from a static arc
%  spun on a turntable. Called with kEnd = N it asserts the whole trajectory;
%  called with an intermediate kEnd it asserts a proper prefix, and a function
%  that drew the full arc in every frame fails the second call.
%
%  The expected coordinates are the arrays section 11 built by hand from traj,
%  passed in whole and indexed here, so the movie is compared against the same
%  independent arithmetic the still figure is.
%
%% Inputs:
%
%  frm              Struct                      One captured frame: fields ph
%                                               [1 x P] phase index per drawn
%                                               segment and xyz [1 x P] cell of
%                                               [M x 3] drawn points (km)
%
%  phIdx            [N x 1]                     Phase index per sample
%
%  phList           [1 x P]                     Phase indices present, sorted
%
%  kEnd             [1 x 1]                     Last sample this frame may show
%
%  xWant            [N x 1]                     Expected x of every sample (km)
%
%  yWant            [N x 1]                     Expected y of every sample (km)
%
%  zWant            [N x 1]                     Expected z of every sample (km)
%
%  what             Char [1 x n]                Frame name, for the message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(numel(frm.xyz) == numel(phList), ...
        '%s frame holds %d track segments for %d phases.', ...
        what,numel(frm.xyz),numel(phList));
    for kp = 1:numel(phList)
               kSg = find(frm.ph == phList(kp),1);
        assert(~isempty(kSg), ...
            '%s frame has no track segment tagged for phase %d.',what,phList(kp));
               sel = phaseSpan(phIdx,phList,kp);
               sel = sel(sel <= kEnd);
               got = frm.xyz{kSg};
        assert(size(got,1) == numel(sel), ...
            ['%s frame, phase %d: %d points drawn against the %d samples the ' ...
             'trajectory has up to sample %d. The track is not growing with ' ...
             'the flight.'],what,phList(kp),size(got,1),numel(sel),kEnd);
        if isempty(sel)
            continue;
        end
        assertSeries(got(:,1),xWant(sel), ...
            sprintf('%s frame, phase %d x (km)',what,phList(kp)));
        assertSeries(got(:,2),yWant(sel), ...
            sprintf('%s frame, phase %d y (km)',what,phList(kp)));
        assertSeries(got(:,3),zWant(sel), ...
            sprintf('%s frame, phase %d z (km)',what,phList(kp)));
    end
end

function traj = degenTraj(rE,kind)
%% Purpose:
%
%  One of the three trajectories that break the movie's geometry rather than
%  its arithmetic. None is a solution of anything; each is the smallest input
%  that reaches one degenerate branch:
%
%    'coincident'  Launch and impact at the SAME point, so the trajectory
%                  plane is undefined and the inset's downrange span is zero.
%    'antipodal'   Endpoints exactly opposite on the sphere, so the plane
%                  normal AND the endpoint sum both vanish -- the sum into
%                  roundoff rather than into zero, which is the harder case.
%    'flat'        A flight that moves but never leaves the surface, so the
%                  inset's altitude span is zero.
%
%% Inputs:
%
%  rE               [1 x 1]                     Reference sphere radius (m)
%
%  kind             Char [1 x n]                'coincident', 'antipodal' or
%                                               'flat'
%
%% Outputs:
%
%  traj             Struct                      Fields t [N x 1] (s),
%                                               x [N x 6], u [N x 2] (rad),
%                                               phaseIdx [N x 1] all ones,
%                                               junction [0 x 1] struct
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = 12;
            traj.t = linspace(0,600,nS)';
    switch lower(kind)
        case 'coincident'
               lat = deg2rad(20).*ones(nS,1);
               lon = deg2rad(-155).*ones(nS,1);
               hKm = linspace(40,10,nS)';
        case 'antipodal'
               lat = deg2rad(linspace(10,-10,nS))';
               lon = deg2rad(linspace(20,-160,nS))';
               hKm = linspace(40,10,nS)';
        case 'flat'
               lat = zeros(nS,1);
               lon = deg2rad(linspace(0,20,nS))';
               hKm = zeros(nS,1);
        otherwise
            error('degenTraj: unknown kind "%s".',kind);
    end
            traj.x = [rE + 1000.*hKm, lon, lat, ...
                      3000.*ones(nS,1), zeros(nS,1), zeros(nS,1)];
            traj.u = zeros(nS,2);
     traj.phaseIdx = ones(nS,1);
     traj.junction = repmat(struct('t',[],'x',[]),0,1);
end

function u = unitOf3(v)
%% Purpose:
%
%  Normalise a three-vector. Written out rather than taken from norm() because
%  this library does not use norm anywhere.
%
%% Inputs:
%
%  v                [3 x 1]                     Any non-zero vector
%
%% Outputs:
%
%  u                [3 x 1]                     v scaled to unit length
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               mag = sqrt(sum(v(:).^2));
    assert(mag > 0,'cannot normalise a zero vector.');
                 u = v(:)./mag;
end

function ang = angBetween(a,b)
%% Purpose:
%
%  The angle between two three-vectors, in degrees, by atan2 of the cross and
%  the dot rather than by acos of the dot: acos loses all its precision
%  exactly where these checks need it, on nearly parallel directions.
%
%% Inputs:
%
%  a                [3 x 1]                     First vector
%
%  b                [3 x 1]                     Second vector
%
%% Outputs:
%
%  ang              [1 x 1]                     Angle between them (deg)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                cr = cross(a(:),b(:));
               ang = atan2d(sqrt(sum(cr.^2)),dot(a(:),b(:)));
end

function d = camDirOf(azEl)
%% Purpose:
%
%  The unit vector a MATLAB view(az,el) looks FROM, inverting the convention
%  coorbital.viz.globeMovie uses to set it: a camera over longitude L and
%  latitude B is view(L + 90, B). This is how a rendered frame's camera is
%  turned back into geometry the trajectory can be compared against.
%
%% Inputs:
%
%  azEl             [1 x 2]                     Azimuth and elevation as read
%                                               from an axes View (deg)
%
%% Outputs:
%
%  d                [3 x 1]                     Unit vector from the centre of
%                                               the sphere towards the camera
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               azR = deg2rad(azEl(1) - 90);
               elR = deg2rad(azEl(2));
                 d = [cos(elR).*cos(azR); cos(elR).*sin(azR); sin(elR)];
end

function assertCamera(azEl,cWant,what)
%% Purpose:
%
%  Assert that a rendered frame's camera points along an expected direction.
%  Compared as an ANGLE between two unit vectors rather than as two numbers,
%  because azimuth is periodic and degenerate at the poles, and a test that
%  compared the raw pair would fail on wrapping and pass on nonsense.
%
%% Inputs:
%
%  azEl             [1 x 2]                     Axes View as rendered (deg)
%
%  cWant            [3 x 1]                     Expected camera direction
%
%  what             Char [1 x n]                Name of the case, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               got = camDirOf(azEl);
               ang = angBetween(got,cWant);
    assert(ang < 1e-6, ...
        ['%s aims %.6f deg away from the direction the arc frame requires. ' ...
         'Rendered view was [%.4f %.4f] deg.'],what,ang,azEl(1),azEl(2));
end

function assertAngle(got,want,tol,what)
%% Purpose:
%
%  Assert two angles are equal in degrees, comparing them wrapped into
%  (-180,180]. MATLAB stores an axes azimuth wrapped, so a camera correctly
%  placed at 370 deg reads back as 10 and a raw difference would call that a
%  360 degree error.
%
%% Inputs:
%
%  got              [1 x 1]                     Angle read from the figure (deg)
%
%  want             [1 x 1]                     Expected angle (deg)
%
%  tol              [1 x 1]                     Absolute tolerance (deg)
%
%  what             Char [1 x n]                Name of the quantity, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               dif = mod(got - want + 180,360) - 180;
    assert(abs(dif) < tol, ...
        '%s: got %.9g against %.9g, %.3e deg apart, budget %.1e', ...
        what,got,want,abs(dif),tol);
end

function [latM,lonM] = gcMidpoint(lat1,lon1,lat2,lon2)
%% Purpose:
%
%  The great-circle midpoint of two points on a sphere, by the standard
%  navigation formula. Written this way ON PURPOSE: coorbital.viz.globeMovie
%  reaches the same point by normalising the SUM of the two Cartesian unit
%  vectors, so an error common to both sides is not available, and comparing
%  the two is a genuine check rather than a restatement.
%
%% Inputs:
%
%  lat1, lon1       [1 x 1]                     First point (rad)
%
%  lat2, lon2       [1 x 1]                     Second point (rad)
%
%% Outputs:
%
%  latM             [1 x 1]                     Midpoint latitude (rad)
%
%  lonM             [1 x 1]                     Midpoint longitude (rad)
%
%% REFERENCES:
%   [1] Veness, C., "Calculating distance, bearing and more between
%       Latitude/Longitude points," Movable Type Scripts, midpoint formula.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               dLn = lon2 - lon1;
                bX = cos(lat2).*cos(dLn);
                bY = cos(lat2).*sin(dLn);
              latM = atan2(sin(lat1) + sin(lat2), ...
                           sqrt((cos(lat1) + bX).^2 + bY.^2));
              lonM = lon1 + atan2(bY,cos(lat1) + bX);
end

function hAx = axesInOrder(hFig)
%% Purpose:
%
%  Return a figure's axes ordered top to bottom on the page. findobj hands
%  them back in reverse creation order, which is an implementation detail of
%  the graphics system rather than anything a test should depend on; the
%  panels of a stacked profile figure have to be identified by WHERE they are.
%  Legends and colorbars are not axes and do not appear here.
%
%% Inputs:
%
%  hFig             [1 x 1] figure              Figure to search
%
%% Outputs:
%
%  hAx              [1 x M] axes                Axes, topmost first
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hAx = findobj(hFig,'Type','axes');
    assert(~isempty(hAx),'the figure carries no axes at all.');
              yTop = zeros(1,numel(hAx));
    for ka = 1:numel(hAx)
              oPos = hAx(ka).OuterPosition;
          yTop(ka) = oPos(2);
    end
           [~,ord] = sort(yTop,'descend');
               hAx = reshape(hAx(ord),1,[]);
end

function hLn = lineByTag(hAx,tagTxt)
%% Purpose:
%
%  All line objects in one axes carrying a given tag, as a row.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to search
%
%  tagTxt           Char [1 x n]                Tag to match
%
%% Outputs:
%
%  hLn              [1 x M] line                Matching lines, possibly empty
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hLn = findobj(hAx,'Type','line','Tag',tagTxt);
               hLn = reshape(hLn,1,[]);
end

function txt = labelText(str)
%% Purpose:
%
%  Flatten an axis label or title into one char row. MATLAB allows a label to
%  be char, string or a cellstr of stacked lines, and a test that assumed one
%  of the three would fail for the wrong reason.
%
%% Inputs:
%
%  str              Char / string / cellstr     Label as read from the object
%
%% Outputs:
%
%  txt              Char [1 x n]                One-line rendering
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if iscell(str) || (isstring(str) && ~isscalar(str))
               txt = char(strjoin(string(str),' '));
    else
               txt = char(str);
    end
end

function assertLabelled(hAx,what)
%% Purpose:
%
%  Every axes in a figure must carry a non-empty XLabel, YLabel and Title. An
%  unlabelled axis is not a stylistic lapse; it is a plot whose reader cannot
%  tell what is being shown or in what unit.
%
%% Inputs:
%
%  hAx              [1 x M] axes                Axes to check
%
%  what             Char [1 x n]                Name of the function under
%                                               test, for the message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    for ka = 1:numel(hAx)
              want = {'XLabel','YLabel','Title'};
        for kw = 1:numel(want)
               txt = labelText(hAx(ka).(want{kw}).String);
            assert(~isempty(strtrim(txt)), ...
                '%s axes %d has an empty %s.',what,ka,want{kw});
        end
    end
end

function assertUnit(str,tok,what)
%% Purpose:
%
%  Assert that an axis label contains a required token -- normally the unit
%  the axis is displaying. This is the check that catches a panel showing
%  kilometres under a label that says only "altitude", which is the defect a
%  reader has no way of detecting from the figure.
%
%% Inputs:
%
%  str              Char / cellstr              Label as read from the object
%
%  tok              Char [1 x n]                Token that must appear
%
%  what             Char [1 x m]                Name of the axis, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               txt = labelText(str);
    assert(contains(txt,tok), ...
        '%s is labelled "%s", which does not state "%s".',what,txt,tok);
end

function assertLegend(hFig,want)
%% Purpose:
%
%  Assert that a figure's legend carries exactly the entries expected, in
%  order. The end-point entries are caller-settable, so a legend that ignored
%  StartName and EndName would draw a correct picture under wrong names --
%  and "launch" on a glide that begins at an entry interface is wrong.
%
%% Inputs:
%
%  hFig             [1 x 1] figure              Figure to search
%
%  want             Cellstr [1 x M]             Expected legend entries
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hLg = findobj(hFig,'Type','legend');
    assert(isscalar(hLg),'expected one legend, found %d.',numel(hLg));
               got = cellstr(string(hLg.String));
    assert(isequal(reshape(got,1,[]),reshape(want,1,[])), ...
        'legend reads {%s} against the expected {%s}.', ...
        strjoin(got,', '),strjoin(want,', '));
end

function assertMarker(hAx,tagTxt,xWant,yWant,what)
%% Purpose:
%
%  Assert that a tagged single-point marker exists in an axes and sits at the
%  expected coordinates. Existence alone is not enough: a launch marker drawn
%  at the impact point, or at the origin, is present and wrong.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to search
%
%  tagTxt           Char [1 x n]                Tag of the marker
%
%  xWant            [1 x 1]                     Expected x coordinate, in the
%                                               axis's displayed unit
%
%  yWant            [1 x 1]                     Expected y coordinate
%
%  what             Char [1 x m]                Name of the marker, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hMk = findobj(hAx,'Tag',tagTxt);
    assert(isscalar(hMk),'expected exactly one %s marker, found %d.',what,numel(hMk));
    assert(isscalar(hMk.XData) && isscalar(hMk.YData), ...
        'the %s marker holds %d points; a marker is one point.',what,numel(hMk.XData));
    assert(abs(hMk.XData - xWant) < 1e-9, ...
        'the %s marker is at x = %.9f against the trajectory''s %.9f', ...
        what,hMk.XData,xWant);
    assert(abs(hMk.YData - yWant) < 1e-9, ...
        'the %s marker is at y = %.9f against the trajectory''s %.9f', ...
        what,hMk.YData,yWant);
end

function assertSurfaceMarker(hAx,tagTxt,rE,latRad,lonRad,what)
%% Purpose:
%
%  Assert that a tagged marker in a three-dimensional globe figure sits on the
%  sphere at the given latitude and longitude, in kilometres. Written against
%  the spherical-to-Cartesian conversion done here rather than read back from
%  the figure, so a transposition of latitude and longitude fails.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to search
%
%  tagTxt           Char [1 x n]                Tag of the marker
%
%  rE               [1 x 1]                     Sphere radius (m)
%
%  latRad           [1 x 1]                     Expected latitude (rad)
%
%  lonRad           [1 x 1]                     Expected longitude (rad)
%
%  what             Char [1 x m]                Name of the marker, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hMk = findobj(hAx,'Tag',tagTxt);
    assert(isscalar(hMk),'expected exactly one %s marker, found %d.',what,numel(hMk));
             xWant = rE.*cos(latRad).*cos(lonRad)./1000;
             yWant = rE.*cos(latRad).*sin(lonRad)./1000;
             zWant = rE.*sin(latRad)./1000;
               got = [hMk.XData hMk.YData hMk.ZData];
              want = [xWant yWant zWant];
               err = max(abs(got - want));
    assert(err < 1e-6, ...
        ['the %s marker is at [%.4f %.4f %.4f] km against [%.4f %.4f %.4f] km ' ...
         'from its latitude and longitude; %.3e km apart.'], ...
        what,got(1),got(2),got(3),want(1),want(2),want(3),err);
end

function assertThrows(fn,what)
%% Purpose:
%
%  Assert that a call raises. Used for the caller errors that would otherwise
%  produce a plausible but wrong figure -- a misspelt channel name silently
%  dropped, a mass channel on a state that has no mass, an axes array of the
%  wrong length half-filled, an unlabelled extra panel.
%
%% Inputs:
%
%  fn               Function handle             Zero-argument call to make
%
%  what             Char [1 x n]                Description of the bad input,
%                                               for the message
%
%% Outputs:
%
%  none                                         Throws when fn does NOT
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             threw = false;
    try
        fn();
    catch
             threw = true;
    end
    assert(threw,'%s must be refused, not accepted.',what);
end

function assertSeries(got,want,what)
%% Purpose:
%
%  Assert that a series read back out of a figure equals the series recomputed
%  from the trajectory. The budget is relative to the largest magnitude in the
%  reference, so it is meaningful for a channel that passes through zero, and
%  it is tight enough that no real defect can hide under it: the loosest thing
%  it has to tolerate is a reassociation of the same arithmetic, and the
%  smallest thing it has to catch -- a doubled kilometre conversion -- is off
%  by three orders of magnitude.
%
%% Inputs:
%
%  got              [N x 1]                     Series read from the figure
%
%  want             [N x 1]                     Series recomputed from traj
%
%  what             Char [1 x n]                Name of the series, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(numel(got) == numel(want), ...
        '%s: the figure holds %d points against %d in the trajectory.', ...
        what,numel(got),numel(want));
             scale = max(abs(want));
    if scale == 0
             scale = 1;
    end
               err = max(abs(got(:) - want(:)))./scale;
    assert(err < 1e-12, ...
        '%s: the plotted series differs from the trajectory by %.3e relative.', ...
        what,err);
end
