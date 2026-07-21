function animate_three_way(mode)
% ANIMATE_THREE_WAY  Compare min-time, min-energy, min-fuel trajectories.
%
% Rotating CR3BP frame. Three solutions overlaid, animated by NORMALIZED
% progress (each swept 0->1 of its own arc so they finish together despite
% different durations):
%   min-time   (orange) : full GTO->tulip, throttle always 1 (always burns)
%   min-energy (blue)   : full GTO->tulip, smooth saturated-ramp throttle
%   min-fuel   (green)  : ARRIVAL LEG ONLY (tau>=4), bang-bang throttle
% The throttle subplot (vs normalized progress) is the story: flat-on vs
% ramp vs bang-bang. Min-fuel is a partial arc by design -- the full 40-rev
% min-fuel spiral defeats both methods; only the arrival leg is solved.
%
% INPUTS:
%   mode - 'preview' saves 3 stills; 'movie' renders MP4 + GIF [default movie]

if nargin < 1, mode = 'movie'; end
sp = [fileparts(mfilename('fullpath')) filesep];
S  = load([sp 'compare_data.mat']);

cMT=[0.90 0.50 0.10]; cME=[0.15 0.35 0.85]; cMF=[0.15 0.60 0.25];
sets = { struct('r',S.rMT,'s',S.sMT,'c',cMT,'name','min-time (always burn)'), ...
         struct('r',S.rME,'s',S.sME,'c',cME,'name','min-energy (ramp)'), ...
         struct('r',S.rMF,'s',S.sMF,'c',cMF,'name','min-fuel (bang-bang, arrival leg)') };

allP=[S.rMT S.rME S.rMF].';
pad=0.05*(max(allP)-min(allP)+eps);
xl=[min(allP(:,1))-pad(1) max(allP(:,1))+pad(1)];
yl=[min(allP(:,2))-pad(2) max(allP(:,2))+pad(2)];
zl=[min(allP(:,3))-pad(3) max(allP(:,3))+pad(3)];

if strcmp(mode,'preview'), pFrac=[0.2 0.6 1.0]; tags={'early','mid','late'};
else, pFrac=linspace(0.01,1,260); end

fig=figure('Color','w','Position',[100 100 1000 800],'Visible','off');
try, theme(fig,'light'); catch, end     % R2025b -batch defaults dark
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
axT=nexttile(tl,1,[3 1]); hold(axT,'on'); grid(axT,'on'); box(axT,'on');
if ~isempty(S.yTul)
    plot3(axT,S.yTul(:,1),S.yTul(:,2),S.yTul(:,3),'-','Color',[0.6 0.78 0.6],'LineWidth',0.7);
end
plot3(axT,S.earth(1),S.earth(2),S.earth(3),'o','MarkerFaceColor',[0.1 0.35 0.8],'MarkerEdgeColor','k','MarkerSize',11);
plot3(axT,S.moon(1),S.moon(2),S.moon(3),'o','MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','k','MarkerSize',8);
text(axT,S.earth(1),S.earth(2),S.earth(3)+0.04,'Earth','FontSize',9);
text(axT,S.moon(1),S.moon(2),S.moon(3)+0.04,'Moon','FontSize',9);
xlim(axT,xl); ylim(axT,yl); zlim(axT,zl); view(axT,-37,22); daspect(axT,[1 1 1]);
xlabel(axT,'x (rot, ND)'); ylabel(axT,'y'); zlabel(axT,'z');
title(axT,'Three low-thrust solutions, rotating CR3BP frame');

hTraj=gobjects(1,3); hSC=gobjects(1,3);
for k=1:3
    plot3(axT,sets{k}.r(1,:),sets{k}.r(2,:),sets{k}.r(3,:),'-','Color',[sets{k}.c 0.25],'LineWidth',0.5); %#ok<*NASGU>
    hTraj(k)=plot3(axT,nan,nan,nan,'-','Color',sets{k}.c,'LineWidth',1.8);
    hSC(k)=plot3(axT,nan,nan,nan,'o','MarkerFaceColor',sets{k}.c,'MarkerEdgeColor','k','MarkerSize',6);
end
legend(hTraj,{sets{1}.name,sets{2}.name,sets{3}.name},'Location','northeast','FontSize',8);

axS=nexttile(tl,4); hold(axS,'on'); grid(axS,'on'); box(axS,'on');
for k=1:3
    pk=linspace(0,1,numel(sets{k}.s));
    plot(axS,pk,sets{k}.s,'-','Color',sets{k}.c,'LineWidth',1.3);
end
ylim(axS,[-0.05 1.12]); xlim(axS,[0 1]);
xlabel(axS,'normalized progress along transfer'); ylabel(axS,'throttle s');
title(axS,'throttle: always-on (min-time) vs ramp (min-energy) vs bang-bang (min-fuel)');
hCur=plot(axS,[0 0],[-0.05 1.12],'k-','LineWidth',1.0);

if ~strcmp(mode,'preview')
    vw=VideoWriter([sp 'three_way_comparison'],'MPEG-4'); vw.FrameRate=24; vw.Quality=95; open(vw);
    gifFile=[sp 'three_way_comparison.gif']; gifStride=2; gifW=640; gifDelay=1/12; gifMap=[];
end

fc=0; vidHW=[];
for p=pFrac
    fc=fc+1;
    for k=1:3
        n=size(sets{k}.r,2); idx=max(2,round(p*(n-1))+1);
        set(hTraj(k),'XData',sets{k}.r(1,1:idx),'YData',sets{k}.r(2,1:idx),'ZData',sets{k}.r(3,1:idx));
        set(hSC(k),'XData',sets{k}.r(1,idx),'YData',sets{k}.r(2,idx),'ZData',sets{k}.r(3,idx));
    end
    set(hCur,'XData',[p p]);
    if strcmp(mode,'preview')
        exportgraphics(fig,sprintf('%spreview_three_%s.png',sp,tags{fc}),'Resolution',130);
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
if ~strcmp(mode,'preview'), close(vw); fprintf('WROTE three_way_comparison.mp4 + .gif (%d frames)\n',numel(pFrac));
else fprintf('WROTE three-way preview PNGs\n'); end
close(fig);
end
