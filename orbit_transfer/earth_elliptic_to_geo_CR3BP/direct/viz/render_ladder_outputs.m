function render_ladder_outputs()
% RENDER_LADDER_OUTPUTS  The usual plots + movie for every thrust-ladder rung.
%
% For each certified rung product cr3bp_T<tag>N_phi0_fuel.mat, renders the
% campaign's standard trio into results/:
%   <run>_traj.png      - top-down 2D + 3D trajectory (red burn / blue coast)
%   <run>_throttle.png  - throttle + mass vs time
%   <run>_movie.mp4/gif - house transfer movie (honest CR3BP min-fuel label)
% Rungs whose movie already exists are skipped (delete to re-render). Deep
% rungs render with a reduced densify factor (nDense 3) to keep frame cost
% sane at N ~ 3500 nodes / ~700 revs.
%
% INPUTS:  none        OUTPUTS: none (files written; summary printed)
% REFERENCES:
%   [1] run_cr3bp_geo.m sections 4-5 (the plot/movie blocks mirrored here).
%   [2] ../../earth_elliptic_to_geo/direct/viz/transfer_movie.m (renderer).
here   = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup_paths.m'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'earth_elliptic_to_geo', 'direct', 'viz'));
resDir = fullfile(here, '..', 'results');
rungs  = {10,'T10N',8; 5,'T5N',8; 2.5,'T2p5N',6; 1,'T1N',4; 0.5,'T0p5N',4; 0.2,'T0p2N',3; 0.1,'T0p1N',3};
for k = 1:size(rungs,1)
    tag  = sprintf('cr3bp_%s_phi0_fuel', rungs{k,2});
    pf   = fullfile(resDir, [tag '.mat']);
    if ~isfile(pf), fprintf('SKIP %s (no product)\n', tag); continue; end
    if isfile(fullfile(resDir, [tag '_movie.mp4']))
        fprintf('SKIP %s (movie exists)\n', tag); continue;
    end
    S   = load(pf);
    par = kepler_lt_params(S.fp.thrustN, S.fp.m0kg, S.fp.ispS);
    cert = table3_certified(S.fp.thrustN);
    ctfEff = S.fp.tfTarget / cert.tfmin;
    cart = mee_res_to_cart_res(S.X, S.U, S.dL, S.sigma, S.fp.thrustN, ctfEff, 1, rungs{k,3});
    r = cart.fuel.X(1:3,:);  sTh = cart.fuel.U(4,:);  burn = sTh > 0.05;
    tDays = S.t_days;
    % --- traj png ---
    fig = figure('Visible','off','Position',[60 60 1280 560]);
    thG = linspace(0,2*pi,361);
    for sp = 1:2
        ax = subplot(1,2,sp);  hold(ax,'on');
        plot3(ax, cos(thG), sin(thG), 0*thG, '-', 'Color',[.2 .6 .2]);
        plot3(ax, r(1,burn), r(2,burn), r(3,burn), '.', 'Color',[0.85 0.2 0.15], 'MarkerSize',3);
        plot3(ax, r(1,~burn), r(2,~burn), r(3,~burn), '.', 'Color',[0.2 0.35 0.8], 'MarkerSize',3);
        plot3(ax, 0,0,0,'o','MarkerFaceColor',[.2 .4 .9],'MarkerSize',9);
        axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'x (ND)'); ylabel(ax,'y (ND)');
        if sp==1, view(ax,2); title(ax,'top-down'); else, view(ax,3); zlabel(ax,'z'); title(ax,'3D'); end
    end
    sgtitle(fig, sprintf('%s:  T=%g N, min-fuel (CR3BP), %d switches, m_f=%.4f kg', ...
        strrep(tag,'_','\_'), S.fp.thrustN, S.switches, S.m_f_kg));
    exportgraphics(fig, fullfile(resDir,[tag '_traj.png']), 'Resolution', 150); close(fig);
    % --- throttle png ---
    fig2 = figure('Visible','off','Position',[60 60 1100 500]);
    subplot(2,1,1); stairs(tDays, S.throttle, 'k','LineWidth',0.8); ylim([-0.05 1.05]); grid on
    ylabel('throttle \delta'); title(sprintf('T=%g N thrust profile (%d switches, nodal)', S.fp.thrustN, S.switches));
    subplot(2,1,2); plot(tDays, S.fp.m0kg*S.X(6,:), 'r','LineWidth',1.1); grid on
    xlabel('time (days)'); ylabel('mass (kg)');
    exportgraphics(fig2, fullfile(resDir,[tag '_throttle.png']), 'Resolution', 150); close(fig2);
    % --- movie ---
    cart.cfg.label = sprintf('min-fuel (CR3BP)');
    transfer_movie(cart, fullfile(resDir, [tag '_movie']));
    fprintf('DONE %s\n', tag);
end
fprintf('RENDER_LADDER_OUTPUTS complete\n');
end
