% MAKE_FIGS  Publication figures for the min-fuel method note (a) + briefing
% slides: transfer geometry (tulip + ELFO, with GTO/insertion velocity arrows)
% and the dV-vs-transfer-time fronts, plus the certified trajectory+throttle.
% Clean white-background PDFs/PNGs at 300 dpi.
here = fileparts(mfilename('fullpath'));  cd(here);
addpath('../../direct/sundman_minfuel'); addpath('../../../GTO_ELFO/direct/elfo'); addpath('../../../cr3bp_common');
run('../../direct/sundman_minfuel/setup_paths.m');
p = cr3bp_lt_params(25e-3, 15, 2100);  mu = p.muStar;
Earth = [-mu, 0];  Moon = [1-mu, 0];  outdir = here;

C = load('../../sundman_minfuel/sundman_minfuel_certified.mat');   % certified tulip min-fuel
E = load('../../elfo/results/mintime_elfo.mat');                    % an ELFO transfer
[~,~,trT] = gto_tulip_endpoints(p);
[~,~,trE] = gto_elfo_endpoints(p, struct('point','apolune'));
[rv0, rvfT] = insertion_states('tulip','campaign');
[~,   rvfE] = insertion_states('elfo','nearest');

% ---- Figures 1/2: transfer geometry -----------------------------------------
G(1) = struct('name','tulip','spiral',C.out.X(1:2,:),'orbit',trT(:,1:2),'rvf',rvfT,'fname','fig_geometry_tulip');
G(2) = struct('name','ELFO', 'spiral',E.X(1:2,:),    'orbit',trE(:,1:2),'rvf',rvfE,'fname','fig_geometry_elfo');
for k = 1:2
    g = G(k);  rvf = g.rvf;  spiral = g.spiral;  orbit = g.orbit;
    f = figure('Color','w','Position',[100 100 760 620]); ax=axes(f); hold(ax,'on'); box(ax,'on');
    hS = plot(ax, spiral(1,:), spiral(2,:), '-', 'Color',[0.55 0.75 0.95],'LineWidth',0.6);
    hO = plot(ax, orbit(:,1), orbit(:,2), '-', 'Color',[0.85 0.55 0.20],'LineWidth',1.6);
    plot(ax, Earth(1),Earth(2),'o','MarkerFaceColor',[0.20 0.45 0.85],'MarkerEdgeColor','k','MarkerSize',13);
    plot(ax, Moon(1), Moon(2), 'o','MarkerFaceColor',[0.55 0.55 0.55],'MarkerEdgeColor','k','MarkerSize',9);
    text(ax, Earth(1), Earth(2)-0.07,'Earth','HorizontalAlignment','center','FontSize',11);
    text(ax, Moon(1),  Moon(2)-0.07, 'Moon', 'HorizontalAlignment','center','FontSize',11);
    L = 0.11;  vg = rv0(4:5)/norm(rv0(4:5));  vf = rvf(4:5)/norm(rvf(4:5));
    plot(ax, rv0(1),rv0(2),'^','MarkerFaceColor',[0.10 0.55 0.20],'MarkerEdgeColor','k','MarkerSize',9);
    quiver(ax, rv0(1),rv0(2), L*vg(1), L*vg(2), 0,'Color',[0.10 0.55 0.20],'LineWidth',1.8,'MaxHeadSize',2);
    text(ax, rv0(1)+0.02, rv0(2)+0.06, sprintf('GTO departure (v=%.1f ND)',norm(rv0(4:6))),'FontSize',9,'Color',[0.10 0.45 0.15]);
    plot(ax, rvf(1),rvf(2),'p','MarkerFaceColor',[0.75 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',13);
    quiver(ax, rvf(1),rvf(2), L*vf(1), L*vf(2), 0,'Color',[0.75 0.15 0.15],'LineWidth',1.8,'MaxHeadSize',2);
    text(ax, rvf(1)+0.02, rvf(2)+0.07, sprintf('%s insertion (v=%.2f ND)',g.name,norm(rvf(4:6))),'FontSize',9,'Color',[0.6 0.1 0.1]);
    axis(ax,'equal'); xlim(ax,[-0.15 1.25]);
    ylo = min([spiral(2,:), orbit(:,2)'])-0.12;  yhi = max([spiral(2,:), orbit(:,2)'])+0.18;
    ylim(ax,[ylo yhi]);
    xlabel(ax,'x (rotating frame, ND)'); ylabel(ax,'y (ND)');
    title(ax, sprintf('GTO \\rightarrow %s transfer geometry (rotating frame, x{-}y)', g.name),'FontSize',13);
    legend(ax,[hS hO],{'transfer spiral (\sim40 rev)','target orbit'},'Location','northwest','FontSize',9);
    set(ax,'FontSize',11);
    exportgraphics(f, fullfile(outdir,[g.fname '.pdf']),'ContentType','image','Resolution',300);
    exportgraphics(f, fullfile(outdir,[g.fname '.png']),'Resolution',300); close(f);
    fprintf('wrote %s\n', g.fname);
end

% ---- Figure 3: tulip dV-tf front (from the campaign record) ------------------
fT  = [1.00 1.12 1.14 1.15 1.20 1.25 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80];
dVT = [4.47 3.83 3.49 3.37 3.24 3.14 2.782 2.647 2.520 2.444 2.434 2.466 2.523 2.591];
certT = [1.12 1.14 1.15 1.20 1.25];
f3 = figure('Color','w','Position',[100 100 780 500]); ax=axes(f3); hold(ax,'on'); box(ax,'on');
plot(ax, fT, dVT, '-','Color',[0.4 0.4 0.4],'LineWidth',1.2);
h1=plot(ax, fT, dVT,'o','MarkerFaceColor',[0.55 0.55 0.55],'MarkerEdgeColor','k','MarkerSize',6);
ic=ismember(fT,certT);
h2=plot(ax, fT(ic),dVT(ic),'o','MarkerFaceColor',[0.20 0.65 0.30],'MarkerEdgeColor','k','MarkerSize',7);
h3=plot(ax, 1.00,4.47,'s','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',9);
[~,im]=min(dVT); h4=plot(ax, fT(im),dVT(im),'v','MarkerFaceColor',[0.75 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',10);
text(ax, fT(im),dVT(im)-0.12,sprintf('best \\Delta V = %.3f km/s (%.2f\\times)',dVT(im),fT(im)),'HorizontalAlignment','center','FontSize',9,'Color',[0.6 0.1 0.1]);
xlabel(ax,'transfer time  t_f / t_f^{min}'); ylabel(ax,'\Delta V  (km/s)');
title(ax,'GTO \rightarrow tulip minimum-fuel front','FontSize',13);
legend(ax,[h1 h2 h3 h4],{'feasible-envelope','PMP-certified band','min-time (0 sw)','front minimum'},'Location','northeast','FontSize',9);
set(ax,'FontSize',11); grid(ax,'on'); ax.GridAlpha=0.15;
exportgraphics(f3, fullfile(outdir,'fig_front_tulip.pdf'),'ContentType','image','Resolution',300);
exportgraphics(f3, fullfile(outdir,'fig_front_tulip.png'),'Resolution',300); close(f3);
fprintf('wrote fig_front_tulip\n');

% ---- Figure 4: ELFO dV-tf front (from the batch summary) --------------------
S = load('../../elfo/results/elfo_batch_summary_minEps0.mat'); res=S.res;
ok=[res.ok] & [res.epsReached]; fE=[res(ok).factor]; dVE=[res(ok).dV]; [fE,ord]=sort(fE); dVE=dVE(ord);
f4 = figure('Color','w','Position',[100 100 780 500]); ax=axes(f4); hold(ax,'on'); box(ax,'on');
plot(ax, fE, dVE, '-','Color',[0.4 0.4 0.4],'LineWidth',1.2);
h1=plot(ax, fE, dVE,'o','MarkerFaceColor',[0.30 0.55 0.85],'MarkerEdgeColor','k','MarkerSize',7);
[~,im]=min(dVE); h2=plot(ax, fE(im),dVE(im),'v','MarkerFaceColor',[0.75 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',10);
text(ax, fE(im),dVE(im)-0.05,sprintf('min \\Delta V = %.3f km/s (%.2f\\times)',dVE(im),fE(im)),'HorizontalAlignment','center','FontSize',9,'Color',[0.6 0.1 0.1]);
xlabel(ax,'transfer time  t_f / t_f^{min} (tulip scale)'); ylabel(ax,'\Delta V  (km/s)');
title(ax,'GTO \rightarrow ELFO minimum-fuel front','FontSize',13);
legend(ax,[h1 h2],{'\epsilon=0 bang-bang','front minimum'},'Location','northeast','FontSize',9);
set(ax,'FontSize',11); grid(ax,'on'); ax.GridAlpha=0.15;
exportgraphics(f4, fullfile(outdir,'fig_front_elfo.pdf'),'ContentType','image','Resolution',300);
exportgraphics(f4, fullfile(outdir,'fig_front_elfo.png'),'Resolution',300); close(f4);
fprintf('wrote fig_front_elfo\n');

% ---- Figure 5: certified tulip trajectory colored by throttle + profile ------
X=C.out.X; U=C.out.U; s=U(4,:); t=X(8,:);
f5 = figure('Color','w','Position',[100 100 1050 460]);
ax1=subplot(1,2,1); hold(ax1,'on'); box(ax1,'on');
scatter(ax1, X(1,:),X(2,:),6,s,'filled'); colormap(ax1,flipud(gray)); cb=colorbar(ax1); cb.Label.String='throttle s'; clim(ax1,[0 1]);
plot(ax1, Earth(1),Earth(2),'o','MarkerFaceColor',[0.20 0.45 0.85],'MarkerEdgeColor','k','MarkerSize',11);
plot(ax1, Moon(1),Moon(2),'o','MarkerFaceColor',[0.55 0.55 0.55],'MarkerEdgeColor','k','MarkerSize',8);
axis(ax1,'equal'); xlabel(ax1,'x (ND)'); ylabel(ax1,'y (ND)');
title(ax1,'Burn near perigee, coast near apogee','FontSize',11); set(ax1,'FontSize',10);
ax2=subplot(1,2,2); hold(ax2,'on'); box(ax2,'on');
plot(ax2, t*p.tStar/86400, s,'-','Color',[0.15 0.15 0.15],'LineWidth',0.8);
xlabel(ax2,'time (days)'); ylabel(ax2,'throttle s'); ylim(ax2,[-0.05 1.05]);
title(ax2,sprintf('Throttle history (%d switches, 99.4%% bang-bang)',C.out.switches),'FontSize',11); set(ax2,'FontSize',10); grid(ax2,'on'); ax2.GridAlpha=0.15;
exportgraphics(f5, fullfile(outdir,'fig_traj_throttle.pdf'),'ContentType','image','Resolution',300);
exportgraphics(f5, fullfile(outdir,'fig_traj_throttle.png'),'Resolution',300); close(f5);
fprintf('wrote fig_traj_throttle\nMAKE_FIGS DONE\n');
