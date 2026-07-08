function animate_coarse_minfuel(mode)
% ANIMATE_COARSE_MINFUEL  Movie of the coarse tf-continuation min-fuel control.
%
% The full-spiral min-fuel control found by the coarse tf-continuation at
% tf = 1.10 x tfMin (6 switches). Rotating CR3BP frame: the trajectory is
% colored by throttle state (red = burn, blue = coast) so the three brief
% coast arcs are visible, with primer-direction thrust arrows during the
% burn, over a synced throttle subplot (the 6-switch bang-bang).
%
% INPUTS:
%   mode - 'preview' saves 3 stills; 'movie' renders MP4 + GIF [default movie]

if nargin < 1, mode = 'movie'; end
sp = [fileparts(mfilename('fullpath')) filesep];
S  = load([sp 'coarse_minfuel_control_solution.mat']);   % r, s, U, tMesh, tf, ...
muStar = 0.012150585609624;
tStar  = S.tStar;

r  = S.r;  s = S.s;  w = S.U(1:3,:);
tDays = S.tMesh(:).'*tStar/86400;
nNode = numel(s);
burn  = s > 0.5;
earth = [-muStar,0,0];  moon = [1-muStar,0,0];

% tulip reference (from the min-fuel movie bundle, if present)
yTul = [];
try, mf = load([sp 'minfuel_movie_data.mat']); yTul = mf.yTul; catch, end

allP = r.';
pad  = 0.05*(max(allP)-min(allP)+eps);
xl=[min(allP(:,1))-pad(1) max(allP(:,1))+pad(1)];
yl=[min(allP(:,2))-pad(2) max(allP(:,2))+pad(2)];
zl=[min(allP(:,3))-pad(3) max(allP(:,3))+pad(3)];

if strcmp(mode,'preview'), frames=[round(0.3*nNode) round(0.65*nNode) nNode]; tags={'early','mid','late'};
else, frames=unique(round(linspace(2,nNode,260))); end

fig=figure('Color','w','Position',[100 100 1000 780],'Visible','off');
try, theme(fig,'light'); catch, end
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
axT=nexttile(tl,1,[3 1]); hold(axT,'on'); grid(axT,'on'); box(axT,'on');
if ~isempty(yTul), plot3(axT,yTul(:,1),yTul(:,2),yTul(:,3),'-','Color',[0.55 0.75 0.55],'LineWidth',0.8); end
plot3(axT,r(1,:),r(2,:),r(3,:),'-','Color',[0.82 0.85 0.9],'LineWidth',0.4);
plot3(axT,earth(1),earth(2),earth(3),'o','MarkerFaceColor',[0.1 0.35 0.8],'MarkerEdgeColor','k','MarkerSize',11);
plot3(axT,moon(1),moon(2),moon(3),'o','MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','k','MarkerSize',8);
text(axT,earth(1),earth(2),earth(3)+0.04,'Earth','FontSize',9);
text(axT,moon(1),moon(2),moon(3)+0.04,'Moon','FontSize',9);
xlim(axT,xl); ylim(axT,yl); zlim(axT,zl); view(axT,-37,22); daspect(axT,[1 1 1]);
xlabel(axT,'x (rot, ND)'); ylabel(axT,'y'); zlabel(axT,'z');
title(axT,sprintf('Min-fuel control, full spiral (tf=1.10x min, %d switches) -- rotating CR3BP frame', S.switches));
hBurn =plot3(axT,nan,nan,nan,'-','Color',[0.85 0.15 0.15],'LineWidth',2.0);
hCoast=plot3(axT,nan,nan,nan,'-','Color',[0.15 0.35 0.85],'LineWidth',2.6);
hSC   =plot3(axT,nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
hArr  =quiver3(axT,nan,nan,nan,nan,nan,nan,0,'Color',[0.85 0.15 0.15],'LineWidth',1.5,'MaxHeadSize',2);
hTxt  =text(axT,xl(1)+0.03*diff(xl),yl(2)-0.05*diff(yl),zl(2),'','FontSize',10,'FontName','Menlo','VerticalAlignment','top');

axS=nexttile(tl,4); hold(axS,'on'); grid(axS,'on'); box(axS,'on');
stairs(axS,tDays,s,'-','Color',[0.4 0.4 0.4],'LineWidth',1.0);
ylim(axS,[-0.05 1.08]); xlim(axS,[0 tDays(end)]);
xlabel(axS,'time (days)'); ylabel(axS,'throttle s');
title(axS,'min-fuel throttle: bang-bang (brief coasts at 3 switch pairs)');
hCur=plot(axS,[0 0],[-0.05 1.08],'k-','LineWidth',1.0);
hDot=plot(axS,nan,nan,'o','MarkerFaceColor',[0.85 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',6);
arrScale=0.06*diff(xl);

if ~strcmp(mode,'preview')
    vw=VideoWriter([sp 'coarse_minfuel_solution'],'MPEG-4'); vw.FrameRate=24; vw.Quality=95; open(vw);
    gifFile=[sp 'coarse_minfuel_solution.gif']; gifStride=2; gifW=640; gifDelay=1/12; gifMap=[];
end

fc=0; vidHW=[];
for f=frames
    fc=fc+1;
    bM=burn(1:f); cM=~bM;
    xb=r(1,1:f); yb=r(2,1:f); zb=r(3,1:f);
    set(hBurn ,'XData',maskv(xb,bM),'YData',maskv(yb,bM),'ZData',maskv(zb,bM));
    set(hCoast,'XData',maskv(xb,cM),'YData',maskv(yb,cM),'ZData',maskv(zb,cM));
    set(hSC,'XData',r(1,f),'YData',r(2,f),'ZData',r(3,f));
    if burn(f)
        a=w(:,f); an=a/max(sqrt(sum(a.^2)),eps)*arrScale;
        set(hArr,'XData',r(1,f),'YData',r(2,f),'ZData',r(3,f),'UData',an(1),'VData',an(2),'WData',an(3),'Visible','on');
        state='BURN';
    else
        set(hArr,'Visible','off'); state='COAST';
    end
    set(hTxt,'String',sprintf(' t = %5.2f d  (%.0f%%)\n throttle s = %.2f\n %s', tDays(f),100*f/nNode,s(f),state));
    set(hCur,'XData',[tDays(f) tDays(f)]); set(hDot,'XData',tDays(f),'YData',s(f));
    if strcmp(mode,'preview')
        exportgraphics(fig,sprintf('%spreview_coarse_%s.png',sp,tags{fc}),'Resolution',130);
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
if ~strcmp(mode,'preview'), close(vw); fprintf('WROTE coarse_minfuel_solution.mp4 + .gif (%d frames)\n',numel(frames));
else fprintf('WROTE coarse preview PNGs\n'); end
close(fig);
end

function v=maskv(x,m), v=x; v(~m)=nan; end
