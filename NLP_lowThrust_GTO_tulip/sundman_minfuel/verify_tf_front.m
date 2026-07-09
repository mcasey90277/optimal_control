function pmp = verify_tf_front(resultsFile, makePlot)
% VERIFY_TF_FRONT  Per-t_f first-order-PMP verification + reliability plot for a
% Delta-V vs transfer-time front.
%
% For each saved solution, recovers the switching function from the discrete
% costates (out.lamDef, the NLP's own KKT duals) and tests the Pontryagin
% first-order conditions -- the VALIDATED empirical-beta route from
% OPTIMALITY_VERIFICATION_PLAN.md sec D: S = 1 - beta*W with a single positive
% scale beta pinned by least squares at the switch intervals (this absorbs the
% trapezoid-weight / tau_f / kappa covector scaling, so no de-scaling is
% needed). A point is a certified first-order PMP EXTREMAL when it meets all of:
%   burn-sign agreement (S<0 on burns) >= 99%,  coast-sign (S>0) >= 99%,
%   per-switch beta spread <= 5%,  primer alignment <= 0.2 deg,
%   |terminal mass costate| <= 1e-3.
% The primer condition is scale-invariant and typically passes everywhere; the
% switching law is the discriminating test that exposes non-extremal points.
%
% INPUTS:
%   resultsFile - .mat with a `results` struct array carrying per-point
%                 .factor .tf_days .dV .switches .edge .X .U .lamDef
%                 .primerAlignDeg  [default tf_front_results.mat]
%   makePlot    - true (default): write a Delta-V-vs-t_f plot colored by PMP
%                 pass/fail (green = certified extremal, grey = not)
%
% OUTPUTS:
%   pmp - struct array per point: .factor .dV .switches .edge .primerDeg
%         .lamMend .beta .spreadPct .burnPct .coastPct .pmpPass
%         (also saved next to resultsFile as *_pmp.mat)
%
% REFERENCES:
%   [1] OPTIMALITY_VERIFICATION_PLAN.md sec D (empirical-beta switching law).
%   [2] Pontryagin minimum principle; primer-vector theory (Lawden).

here = fileparts(mfilename('fullpath'));  addpath(here);
if nargin<1||isempty(resultsFile), resultsFile=fullfile(here,'tf_front_results.mat'); end
if nargin<2||isempty(makePlot), makePlot=true; end
p = cr3bp_lt_params(0.025,15,2100);  c = p.c;
S = load(resultsFile);  R = S.results;  [~,ix]=sort([R.factor]);  R=R(ix);

pmp = struct('factor',{},'dV',{},'switches',{},'edge',{},'primerDeg',{}, ...
             'lamMend',{},'beta',{},'spreadPct',{},'burnPct',{},'coastPct',{},'pmpPass',{});
fprintf('%-6s %-8s %-4s %-6s %-7s %-8s %-6s %-6s %-5s\n', ...
    'fac','dV','sw','edge%','primer','betaSpr','burn%','coast%','PMP');
for e=1:numel(R)
    X=R(e).X; U=R(e).U; lamDef=R(e).lamDef; s=U(4,:); N=size(X,2)-1;
    mMid=0.5*(X(7,1:end-1)+X(7,2:end));
    lamV=lamDef(4:6,:); nlamV=sqrt(sum(lamV.^2,1)); lamM=lamDef(7,:);
    burnI=(s(1:end-1)>0.5)&(s(2:end)>0.5); coastI=(s(1:end-1)<0.5)&(s(2:end)<0.5);
    swI=find(diff(double(s>0.5))~=0); swI=min(swI,N); swI(swI<1)=[];
    W=nlamV.*c./mMid - lamM;                        % global sign -1 (verified)
    if isempty(swI)||sum(W(swI).^2)==0, beta=NaN; spreadPct=Inf;
    else, beta=sum(W(swI))/sum(W(swI).^2); spreadPct=100*std(1./W(swI))/abs(mean(1./W(swI))); end
    Sfun=1-beta*W;
    burnPct=100*mean(Sfun(burnI)<0); coastPct=100*mean(Sfun(coastI)>0);
    primerDeg=R(e).primerAlignDeg; lamMend=lamDef(7,end);
    pass=(burnPct>=99)&&(coastPct>=99)&&(spreadPct<=5)&&(primerDeg<=0.2)&&(abs(lamMend)<=1e-3);
    fprintf('%-6.2f %-8.4f %-4d %-6.1f %-7.3f %-8.1f %-6.1f %-6.1f %-5d\n', ...
        R(e).factor,R(e).dV,R(e).switches,100*R(e).edge,primerDeg,spreadPct,burnPct,coastPct,pass);
    pmp(end+1)=struct('factor',R(e).factor,'dV',R(e).dV,'switches',R(e).switches, ...
        'edge',R(e).edge,'primerDeg',primerDeg,'lamMend',lamMend,'beta',beta, ...
        'spreadPct',spreadPct,'burnPct',burnPct,'coastPct',coastPct,'pmpPass',pass); %#ok<AGROW>
end
[pp,ff]=fileparts(resultsFile); outF=fullfile(pp,[ff '_pmp.mat']); save(outF,'pmp');
fprintf('\n%d points, %d certified first-order PMP extremals. WROTE %s\n', ...
    numel(pmp), sum([pmp.pmpPass]), outF);

if makePlot
    tStar=382981.289129055; tfMin=6.290694;
    d=[R.tf_days]; v=[pmp.dV]; ps=[pmp.pmpPass]==1;
    fig=figure('Color','w','Position',[100 100 840 500],'Visible','off');
    try
        theme(fig,'light');
    catch
    end
    hold on; grid on; box on;
    plot(d(~ps), v(~ps), 'o','Color',[0.6 0.6 0.65],'MarkerFaceColor',[0.86 0.86 0.9],'MarkerSize',8);
    plot(d(ps),  v(ps),  'o','Color',[0.10 0.45 0.15],'MarkerFaceColor',[0.20 0.65 0.25], ...
         'MarkerSize',9,'LineWidth',1.2);
    plot(tfMin*tStar/86400,4.4665,'ks','MarkerFaceColor','k','MarkerSize',9);
    text(tfMin*tStar/86400+0.4,4.4665,'min-time (4.4665, 0 sw)','FontSize',9,'Color',[0.2 0.2 0.2]);
    xlabel('transfer time t_f (days)'); ylabel('\DeltaV (km/s)');
    title('Min-fuel \DeltaV vs t_f  (green = PMP-certified extremal, grey = not certified)');
    legend({'not certified','PMP-certified','min-time'},'Location','northwest','Box','off');
    exportgraphics(fig,fullfile(here,'tf_front_verified.png'),'Resolution',150); close(fig);
    fprintf('WROTE %s\n', fullfile(here,'tf_front_verified.png'));
end
end
