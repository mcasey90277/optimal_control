function plot_branch_map(outPng)
%% Purpose:
%
%   The extremal BRANCH MAP of the DRO -> 7-petal tulip minimum-time
%   transfer in the (thrust, t_f) plane, coloured by the sign of the
%   objective multiplier rho (rho > 0: a min-time PMP candidate; rho < 0: an
%   abnormal connector, not a candidate). FINDINGS section 34. Reads the
%   pseudo-arclength arcs in results/ (normal-chart diagnostic, homogeneous
%   arc from 75.5 mN, and the certified 70 mN branch in both directions).
%
%% Inputs:
%
%  outPng                   char (optional)         output file
%                                                   [results/mintime_branch_map.png]
%
%% Outputs: none (writes the figure)
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
R = fullfile(here, 'results');
if nargin < 1 || isempty(outPng), outPng = fullfile(R, 'mintime_branch_map.png'); end
fig=figure('Position',[100 100 1100 700],'Color','w'); hold on; grid on;
arcs={'mintime_arclength_hom.mat','fast family (from 75.5 mN)','o'; 'mintime_arclength_70mN_dn.mat','70 mN family, down','s'; 'mintime_arclength_70mN_up.mat','70 mN family, up','^'};
for k=1:size(arcs,1)
  f=fullfile(R, arcs{k,1}); if ~isfile(f), continue, end
  L=load(f); A=L.Aout; T=A.T_N*1000; tf=A.tf_days; rho=A.rho;
  pos=rho>=0; neg=rho<0;
  plot(T(pos),tf(pos),[arcs{k,3} '-'],'Color',[0.1 0.4 0.9],'MarkerSize',4,'MarkerFaceColor',[0.1 0.4 0.9],'DisplayName',[arcs{k,2} ' (\rho > 0)']);
  if any(neg), plot(T(neg),tf(neg),[arcs{k,3} '-'],'Color',[0.9 0.3 0.2],'MarkerSize',4,'DisplayName',[arcs{k,2} ' (\rho < 0, abnormal)']); end
  plot(T(1),tf(1),'kp','MarkerSize',12,'MarkerFaceColor','y','HandleVisibility','off');
end
% the normal-chart diagnostic arc and the certified points
L=load(fullfile(R,'mintime_arclength_diag.mat')); A=L.Aout; plot(A.T_N*1000,A.tf_days,'-','Color',[0.5 0.5 0.5],'LineWidth',1.5,'DisplayName','normal chart (\rho = 1), asymptotes at 72.0 mN');
plot(70,26.436,'r*','MarkerSize',16,'LineWidth',2,'DisplayName','certified 70 mN / 26.44 d (direct + ms + tfMin)');
plot(75.495,15.146,'g*','MarkerSize',14,'LineWidth',2,'DisplayName','75.5 mN / 15.15 d (deepest thrust-stepping root)');
xline(70,'k--','70 mN target','LabelVerticalAlignment','bottom','HandleVisibility','off');
yline(18,'k:','abstract: 18 d','HandleVisibility','off');
xlabel('thrust  [mN]'); ylabel('minimum time  t_f  [days]');
title('DRO \rightarrow 7-petal tulip, Isp 900 s, 150 kg: extremal branches in (T, t_f)');
legend('Location','northeast','FontSize',9); set(gca,'FontSize',11);
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('branch map -> %s\n', outPng);
end
