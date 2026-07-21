function animate_energy_solo(mode)
% ANIMATE_ENERGY_SOLO  Movie of the direct min-ENERGY full-spiral solution.
%
% Rotating CR3BP frame (Earth, Moon, tulip fixed). The full GTO -> tulip
% min-energy trajectory is animated with a moving spacecraft marker whose
% color tracks the instantaneous throttle, over a synced throttle subplot
% that shows the smooth saturated RAMP (the min-energy signature -- contrast
% the min-fuel bang-bang).
%
% INPUTS:
%   mode - 'preview' saves 3 stills; 'movie' renders MP4 + GIF [default movie]

if nargin < 1, mode = 'movie'; end
sp = [fileparts(mfilename('fullpath')) filesep];
S  = load([sp 'compare_data.mat']);
tStar = S.tStar;

r  = S.rME;  s = S.sME;  t = S.tME;
tDays = t*tStar/86400;
nNode = numel(t);

% axis limits over the whole spiral
allP = r.';
pad  = 0.05*(max(allP)-min(allP)+eps);
xl=[min(allP(:,1))-pad(1) max(allP(:,1))+pad(1)];
yl=[min(allP(:,2))-pad(2) max(allP(:,2))+pad(2)];
zl=[min(allP(:,3))-pad(3) max(allP(:,3))+pad(3)];

if strcmp(mode,'preview')
    frames=[round(0.15*nNode) round(0.6*nNode) nNode]; tags={'early','mid','late'};
else
    frames=unique(round(linspace(2,nNode,260)));
end

fig=figure('Color','w','Position',[100 100 1000 780],'Visible','off');
try, theme(fig,'light'); catch, end     % R2025b -batch defaults dark
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
axT=nexttile(tl,1,[3 1]); hold(axT,'on'); grid(axT,'on'); box(axT,'on');
if ~isempty(S.yTul)
    plot3(axT,S.yTul(:,1),S.yTul(:,2),S.yTul(:,3),'-','Color',[0.55 0.75 0.55],'LineWidth',0.8);
end
plot3(axT,r(1,:),r(2,:),r(3,:),'-','Color',[0.8 0.85 0.95],'LineWidth',0.5);   % full path faint
plot3(axT,S.earth(1),S.earth(2),S.earth(3),'o','MarkerFaceColor',[0.1 0.35 0.8],'MarkerEdgeColor','k','MarkerSize',11);
plot3(axT,S.moon(1),S.moon(2),S.moon(3),'o','MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','k','MarkerSize',8);
text(axT,S.earth(1),S.earth(2),S.earth(3)+0.04,'Earth','FontSize',9);
text(axT,S.moon(1),S.moon(2),S.moon(3)+0.04,'Moon','FontSize',9);
xlim(axT,xl); ylim(axT,yl); zlim(axT,zl); view(axT,-37,22); daspect(axT,[1 1 1]);
xlabel(axT,'x (rot, ND)'); ylabel(axT,'y'); zlabel(axT,'z');
title(axT,'Minimum-energy full spiral (direct NLP), rotating CR3BP frame');
hTraj=plot3(axT,nan,nan,nan,'-','Color',[0.15 0.35 0.85],'LineWidth',1.8);
hSC  =plot3(axT,nan,nan,nan,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.15 0.35 0.85],'MarkerSize',7);
hTxt =text(axT,xl(1)+0.03*diff(xl),yl(2)-0.05*diff(yl),zl(2),'','FontSize',10,'FontName','Menlo','VerticalAlignment','top');

axS=nexttile(tl,4); hold(axS,'on'); grid(axS,'on'); box(axS,'on');
plot(axS,tDays,s,'-','Color',[0.15 0.35 0.85],'LineWidth',1.2);
ylim(axS,[-0.05 1.08]); xlim(axS,[0 tDays(end)]);
xlabel(axS,'time (days)'); ylabel(axS,'throttle s');
title(axS,'min-energy throttle: continuous saturated ramp (never bang-bang)');
hCur=plot(axS,[0 0],[-0.05 1.08],'k-','LineWidth',1.0);
hDot=plot(axS,nan,nan,'o','MarkerFaceColor',[0.85 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',6);

if ~strcmp(mode,'preview')
    vw=VideoWriter([sp 'minenergy_solution'],'MPEG-4'); vw.FrameRate=24; vw.Quality=95; open(vw);
    gifFile=[sp 'minenergy_solution.gif']; gifStride=2; gifW=640; gifDelay=1/12; gifMap=[];
end

fc=0; vidHW=[];
for f=frames
    fc=fc+1;
    set(hTraj,'XData',r(1,1:f),'YData',r(2,1:f),'ZData',r(3,1:f));
    set(hSC,'XData',r(1,f),'YData',r(2,f),'ZData',r(3,f));
    set(hTxt,'String',sprintf(' t = %5.2f d  (%.0f%%)\n throttle s = %.3f', tDays(f),100*f/nNode,s(f)));
    set(hCur,'XData',[tDays(f) tDays(f)]); set(hDot,'XData',tDays(f),'YData',s(f));
    if strcmp(mode,'preview')
        exportgraphics(fig,sprintf('%spreview_energy_%s.png',sp,tags{fc}),'Resolution',130);
    else
        tmp=[sp 'frame_tmp.png']; exportgraphics(fig,tmp,'Resolution',110); img=imread(tmp);
        if isempty(vidHW), vidHW=2*floor([size(img,1) size(img,2)]/2); end
        if ~isequal([size(img,1) size(img,2)],vidHW), img=imresize(img,vidHW); end
        writeVideo(vw,img);
        if mod(fc-1,gifStride)==0
            gr=round(gifW*size(img,1)/size(img,2)); gimg=imresize(img,[gr gifW]);
            if isempty(gifMap)
                [gi,gifMap]=rgb2ind(gimg,256,'nodither'); imwrite(gi,gifMap,gifFile,'gif','LoopCount',Inf,'DelayTime',gifDelay);
            else
                gi=rgb2ind(gimg,gifMap,'nodither'); imwrite(gi,gifMap,gifFile,'gif','WriteMode','append','DelayTime',gifDelay);
            end
        end
    end
end
if ~strcmp(mode,'preview'), close(vw); fprintf('WROTE minenergy_solution.mp4 + .gif (%d frames)\n',numel(frames));
else fprintf('WROTE energy preview PNGs\n'); end
close(fig);
end
