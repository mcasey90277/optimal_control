function phase_quad_movie()
% PHASE_QUAD_MOVIE  Four-panel synced movie: the lunar-phase comparison at 10 N.
%
% One movie, four panels (phi0 = 0, pi/2, pi, 3pi/2), synchronized in
% PHYSICAL time. Per panel: Earth fixed at origin, GEO ring, the transfer
% trail colored red (burn) / blue (coast), and -- new for this movie -- the
% MOON ORBITING: its true direction angle n_M t + phi0, drawn as a marker on
% a clamped dashed circle at the frame edge (the real distance, 9.12 LU, is
% far outside the trajectory zone; the clamp radius is honest-labeled). Each
% panel carries its running Delta-V and its final Moon-effect verdict
% (dmf vs the same-chain control). At 10 N the Moon sweeps ~70 deg during
% the transfer -- visibly different tidal geometry per panel, same
% spacecraft, diverging outcomes.
%
% Writes results/phase_quad_movie.mp4 (300 frames, 1280x720, H.264-safe
% divisible-by-16 dims) and .gif (150 frames).
%
% INPUTS:  none        OUTPUTS: none (files written)
% REFERENCES:
%   [1] doc/cr3bp_geo_phase1_note.tex sec 'Lunar-phase dependence' (the
%       physics this movie animates); fig_phi_sweep.m (static companion).
%   [2] memory: matlab-movie-diagonal-streaks (div-16 frame rule).
here   = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup_paths.m'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'earth_elliptic_to_geo', 'direct', 'viz'));
resDir = fullfile(here, '..', 'results');
runs = { 0,      'cr3bp_T10N_phi0_fuel.mat',    '\phi_0 = 0';
         pi/2,   'cr3bp_T10N_phiPi2_fuel.mat',  '\phi_0 = \pi/2';
         pi,     'cr3bp_T10N_phiPi_fuel.mat',   '\phi_0 = \pi';
         3*pi/2, 'cr3bp_T10N_phi3Pi2_fuel.mat', '\phi_0 = 3\pi/2' };
C0 = load(fullfile(resDir, 'cr3bp_T10N_2body_control.mat'), 'm_f_kg');
nP = size(runs,1);  P = cell(nP,1);
for k = 1:nP
    S   = load(fullfile(resDir, runs{k,2}));
    par = kepler_lt_params(S.fp.thrustN, S.fp.m0kg, S.fp.ispS);
    cart = mee_res_to_cart_res(S.X, S.U, S.dL, S.sigma, S.fp.thrustN, 1.5, 1, 6);
    pk.r    = cart.fuel.X(1:3,:);
    pk.thr  = cart.fuel.U(4,:);
    pk.tTU  = cart.fuel.X(8,:);
    pk.m    = cart.fuel.X(7,:);
    pk.dV   = par.c*log(1./pk.m)*par.VU_kms;
    pk.nM   = S.fp.nM;  pk.phi0 = runs{k,1};
    pk.dmfG = 1000*(S.m_f_kg - C0.m_f_kg);
    pk.lbl  = runs{k,3};
    pk.tDays= pk.tTU*par.TU_s/86400;
    P{k} = pk;
end
tfTU = P{1}.tTU(end);  TUd = P{1}.tDays(end)/tfTU;
nFrames = 300;  Rcl = 1.45;  thG = linspace(0,2*pi,361);
fig = figure('Visible','off','Position',[40 40 1280 720],'Color','w');
vw = VideoWriter(fullfile(resDir,'phase_quad_movie.mp4'), 'MPEG-4');
vw.FrameRate = 24;  open(vw);
gifFile = fullfile(resDir,'phase_quad_movie.gif');
ax = gobjects(nP,1);
for k = 1:nP, ax(k) = subplot(2,2,k); end
for f = 1:nFrames
    tNow = tfTU * f/nFrames;
    for k = 1:nP
        pk = P{k};  cla(ax(k));  hold(ax(k),'on');
        idx = pk.tTU <= tNow;
        burn = idx & (pk.thr > 0.05);  coast = idx & ~ (pk.thr > 0.05);
        plot(ax(k), cos(thG), sin(thG), '-', 'Color',[.2 .6 .2 .5]);
        plot(ax(k), Rcl*cos(thG), Rcl*sin(thG), 'k:', 'Color',[.5 .5 .5]);
        plot(ax(k), pk.r(1,coast), pk.r(2,coast), '.', 'Color',[0.2 0.35 0.8], 'MarkerSize',3);
        plot(ax(k), pk.r(1,burn),  pk.r(2,burn),  '.', 'Color',[0.85 0.2 0.15], 'MarkerSize',3);
        plot(ax(k), 0,0,'o','MarkerFaceColor',[.2 .4 .9],'MarkerSize',8);
        % current spacecraft position
        ii = find(idx, 1, 'last');
        if ~isempty(ii), plot(ax(k), pk.r(1,ii), pk.r(2,ii), 'ko', 'MarkerFaceColor','k', 'MarkerSize',4); end
        % the orbiting Moon (direction on the clamped circle)
        aM = pk.nM*tNow + pk.phi0;
        plot(ax(k), Rcl*cos(aM), Rcl*sin(aM), 'o', 'MarkerSize',10, ...
             'MarkerFaceColor',[.6 .6 .6], 'Color',[.3 .3 .3]);
        dvNow = pk.dV(max(ii,1));
        sgn = 'helps';  if pk.dmfG < 0, sgn = 'hurts'; end
        title(ax(k), sprintf('%s   \\DeltaV = %.4f km/s   (Moon %s: %+.1f g)', ...
              pk.lbl, dvNow, sgn, pk.dmfG), 'FontSize', 10);
        axis(ax(k), 'equal');  xlim(ax(k), [-1.6 1.6]);  ylim(ax(k), [-1.6 1.6]);
        set(ax(k), 'XTick', [], 'YTick', []);
    end
    sgtitle(fig, sprintf(['lunar-phase comparison -- 10 N min-fuel (CR3BP), t = %5.2f d / %.2f d   ' ...
        '(gray = Moon direction, clamped at r=%.2f; true distance 9.12 LU)'], tNow*TUd, tfTU*TUd, Rcl), 'FontSize', 11);
    drawnow;
    fr = getframe(fig);
    img = fr.cdata;
    h16 = 16*floor(size(img,1)/16);  w16 = 16*floor(size(img,2)/16);
    img = img(1:h16, 1:w16, :);
    writeVideo(vw, img);
    if mod(f,2)==0
        [A,map] = rgb2ind(imresize(img, 0.55), 128);
        if f==2, imwrite(A,map,gifFile,'gif','LoopCount',Inf,'DelayTime',1/15);
        else,    imwrite(A,map,gifFile,'gif','WriteMode','append','DelayTime',1/15); end
    end
end
close(vw);  close(fig);
fprintf('WROTE %s (+.gif)\n', fullfile(resDir,'phase_quad_movie.mp4'));
end
