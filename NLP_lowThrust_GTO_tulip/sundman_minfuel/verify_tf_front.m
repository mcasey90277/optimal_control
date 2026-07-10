function pmp = verify_tf_front(resultsFile, makePlot)
% VERIFY_TF_FRONT  Per-t_f first-order-PMP verification + reliability plot for a
% Delta-V vs transfer-time front.
%
% For each saved solution, recovers the switching function from the discrete
% costates (out.lamDef, the NLP's own KKT duals) and tests the Pontryagin
% first-order conditions -- the empirical-beta route from
% OPTIMALITY_VERIFICATION_PLAN.md sec D: S = 1 - beta*W with a single positive
% scale beta (absorbs the trapezoid-weight / tau_f / kappa covector scaling).
%
% beta is estimated ROBUSTLY: beta = median over switches of the implied
% per-switch scale 1/W_k. The original W^2-weighted least-squares estimator is
% fragile to ONE outlier switch -- diag_beta_checker.m (2026-07-09) showed the
% FIRST switch carries a real, smoothly-t_f-growing anomaly (implied scale
% 0.12-0.93 vs a dead-flat ~1.0 bulk, MAD <0.8%) that dragged the LS beta and
% cascaded fake mid-arc violations across the whole upper band.
%
% TWO-TIER verdict (pmpPass):
%   2 = FULL certified     : all gates below AND first-switch scale within 10%
%   1 = INTERIOR certified : all gates below, first switch anomalous (the
%       open first-switch question -- mesh under-resolution vs genuine local
%       non-extremality -- is flagged, not hidden; ms_band arbiter to settle)
%   0 = not certified
% Gates (with robust beta): burn-sign >= 99%, coast-sign >= 99%, robust
% per-switch scale MAD <= 5%, primer alignment <= 0.2 deg, RELATIVE
% transversality |lam_m(tau_f)|/max|lam_m| <= 1e-3 (scale-invariant).
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
if nargin<1||isempty(resultsFile), resultsFile=fullfile(here,'results','fronts','tf_front_results.mat'); end
if nargin<2||isempty(makePlot), makePlot=true; end
p = cr3bp_lt_params(0.025,15,2100);  c = p.c;
S = load(resultsFile);  R = S.results;  [~,ix]=sort([R.factor]);  R=R(ix);

pmp = struct('factor',{},'dV',{},'switches',{},'edge',{},'primerDeg',{}, ...
             'lamMend',{},'relTrans',{},'beta',{},'madPct',{},'firstSwScale',{}, ...
             'burnPct',{},'coastPct',{},'pmpPass',{});
fprintf('%-6s %-8s %-4s %-6s %-7s %-7s %-7s %-6s %-6s %-5s\n', ...
    'fac','dV','sw','edge%','primer','MAD%','1stSw','burn%','coast%','PMP');
for e=1:numel(R)
    X=R(e).X; U=R(e).U; lamDef=R(e).lamDef; s=U(4,:); N=size(X,2)-1;
    mMid=0.5*(X(7,1:end-1)+X(7,2:end));
    lamV=lamDef(4:6,:); nlamV=sqrt(sum(lamV.^2,1)); lamM=lamDef(7,:);
    burnI=(s(1:end-1)>0.5)&(s(2:end)>0.5); coastI=(s(1:end-1)<0.5)&(s(2:end)<0.5);
    swI=find(diff(double(s>0.5))~=0); swI=min(swI,N); swI(swI<1)=[];
    W=nlamV.*c./mMid - lamM;                        % global sign -1 (verified)
    if isempty(swI)||any(W(swI)==0)
        beta=NaN; madPct=Inf; firstSw=NaN;
    else
        invW=1./W(swI);
        beta=median(invW);                          % ROBUST single scale
        madPct=100*mad(invW,1)/abs(beta);           % robust spread (median-based)
        firstSw=invW(1)/beta;                       % first-switch implied scale
    end
    Sfun=1-beta*W;
    burnPct=100*mean(Sfun(burnI)<0); coastPct=100*mean(Sfun(coastI)>0);
    primerDeg=R(e).primerAlignDeg; lamMend=lamDef(7,end);
    % transversality lam_m(tau_f)=0 checked RELATIVE to the mass-costate's own
    % magnitude (scale-invariant; absolute gate is scale-dependent).
    relTrans=abs(lamMend)/max(max(abs(lamDef(7,:))),eps);
    interior=(burnPct>=99)&&(coastPct>=99)&&(madPct<=5)&&(primerDeg<=0.2)&&(relTrans<=1e-3);
    pass=double(interior)*(1+double(abs(firstSw-1)<=0.10));   % 0 / 1 interior / 2 full
    fprintf('%-6.2f %-8.4f %-4d %-6.1f %-7.3f %-7.2f %-7.2f %-6.1f %-6.1f %-5d\n', ...
        R(e).factor,R(e).dV,R(e).switches,100*R(e).edge,primerDeg,madPct,firstSw,burnPct,coastPct,pass);
    pmp(end+1)=struct('factor',R(e).factor,'dV',R(e).dV,'switches',R(e).switches, ...
        'edge',R(e).edge,'primerDeg',primerDeg,'lamMend',lamMend,'relTrans',relTrans,'beta',beta, ...
        'madPct',madPct,'firstSwScale',firstSw,'burnPct',burnPct,'coastPct',coastPct,'pmpPass',pass); %#ok<AGROW>
end
[pp,ff]=fileparts(resultsFile); outF=fullfile(pp,[ff '_pmp.mat']); save(outF,'pmp');
fprintf('\n%d points: %d FULL certified, %d interior-certified, %d not. WROTE %s\n', ...
    numel(pmp), sum([pmp.pmpPass]==2), sum([pmp.pmpPass]==1), sum([pmp.pmpPass]==0), outF);

if makePlot
    tStar=382981.289129055; tfMin=6.290694;
    d=[R.tf_days]; v=[pmp.dV]; tier=[pmp.pmpPass];
    fig=figure('Color','w','Position',[100 100 840 500],'Visible','off');
    try
        theme(fig,'light');
    catch
    end
    hold on; grid on; box on;
    hh=[]; lbl={};
    if any(tier==0)
        hh(end+1)=plot(d(tier==0), v(tier==0), 'o','Color',[0.6 0.6 0.65],'MarkerFaceColor',[0.86 0.86 0.9],'MarkerSize',8);
        lbl{end+1}='not certified';
    end
    if any(tier==1)
        hh(end+1)=plot(d(tier==1), v(tier==1), 'o','Color',[0.35 0.60 0.30],'MarkerFaceColor',[0.62 0.82 0.55],'MarkerSize',9);
        lbl{end+1}='interior-certified (1st switch flagged)';
    end
    if any(tier==2)
        hh(end+1)=plot(d(tier==2), v(tier==2), 'o','Color',[0.10 0.45 0.15],'MarkerFaceColor',[0.20 0.65 0.25], ...
             'MarkerSize',9,'LineWidth',1.2);
        lbl{end+1}='FULL PMP-certified';
    end
    hh(end+1)=plot(tfMin*tStar/86400,4.4665,'ks','MarkerFaceColor','k','MarkerSize',9);
    lbl{end+1}='min-time';
    text(tfMin*tStar/86400+0.4,4.4665,'min-time (4.4665, 0 sw)','FontSize',9,'Color',[0.2 0.2 0.2]);
    xlabel('transfer time t_f (days)'); ylabel('\DeltaV (km/s)');
    title('Min-fuel \DeltaV vs t_f  (two-tier PMP certification, robust \beta)');
    legend(hh,lbl,'Location','northwest','Box','off');
    exportgraphics(fig,fullfile(here,'tf_front_verified.png'),'Resolution',150); close(fig);
    fprintf('WROTE %s\n', fullfile(here,'tf_front_verified.png'));
end
end
