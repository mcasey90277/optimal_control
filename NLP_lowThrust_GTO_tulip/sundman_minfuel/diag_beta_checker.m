function diag = diag_beta_checker()
% DIAG_BETA_CHECKER  Diagnose WHY the empirical-beta switching-law check fails
% above 1.25x while passing at 1.12-1.25x and 1.85x.
%
% Competing hypotheses (HONEST_EVALUATION / campaign "checker under test"):
%   H1 gate artifact  - violations live in the +/-1 intervals AROUND switches,
%                       where the trapezoid midpoint straddles S=0 and the
%                       sign is genuinely ambiguous at mesh resolution.
%   H2 wrong pattern  - violations fill WHOLE arcs (an arc that should not
%                       burn does, or vice versa): the switch pattern is not
%                       extremal-supported. Signature: arcs >50% violating.
%   H3 scale drift    - the single-beta model breaks: the implied per-switch
%                       scale 1/W_k drifts systematically along the
%                       trajectory instead of scattering around a constant.
% Discriminators, computed per case on the SAVED KKT duals (no solves):
%   - violation distance-to-nearest-switch histogram (H1: mass at <=1)
%   - per-arc violation fraction + count of flipped arcs (H2: flipped arcs)
%   - 1/W_k vs switch index + correlation with time (H3: monotone drift)
%
% Cases: FAIL 1.65x dn (burn ~67%), FAIL 1.35x en (92.6%), PASS 1.20x
% (99.7%), PASS 1.85x (100%).
%
% INPUTS:  none (loads stored results)
% OUTPUTS: diag - struct array per case with the discriminator numbers
%          (also writes results/plots/beta_checker_diag.png)
%
% REFERENCES:
%   [1] verify_tf_front.m (the gate being diagnosed).
%   [2] OPTIMALITY_VERIFICATION_PLAN.md sec D (empirical-beta route).

here = fileparts(mfilename('fullpath'));  addpath(here);
cfg  = minfuel_config();
p    = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);  c = p.c;

% --- load the four cases -----------------------------------------------------
F = load(fullfile(cfg.dirs.fronts,'tf_front_results.mat'));
[~,k20] = min(abs([F.results.factor]-1.20));
cases = struct('name',{},'X',{},'U',{},'lamDef',{});
cases(1) = load_case('PASS 1.20x (up)',  F.results(k20));
cases(2) = load_case('PASS 1.85x',       getfield(load(fullfile(cfg.dirs.minfuel,'legacy_ms_f1850.mat')),'out')); %#ok<GFLD>
cases(3) = load_case('FAIL 1.35x (en)',  getfield(load(fullfile(cfg.dirs.minfuel,'minfuel_f1350_en.mat')),'out')); %#ok<GFLD>
cases(4) = load_case('FAIL 1.65x (dn)',  getfield(load(fullfile(cfg.dirs.minfuel,'minfuel_f1650_dn.mat')),'out')); %#ok<GFLD>

fig = figure('Color','w','Position',[60 60 1240 900],'Visible','off');
try, theme(fig,'light'); catch, end
diag = struct('name',{},'burnPct',{},'violNearSwPct',{},'violMidArcPct',{}, ...
              'nArcs',{},'nFlippedArcs',{},'betaDriftCorr',{},'spreadPct',{});

for e = 1:numel(cases)
    X=cases(e).X; U=cases(e).U; lamDef=cases(e).lamDef; s=U(4,:); N=size(X,2)-1;
    mMid=0.5*(X(7,1:end-1)+X(7,2:end));
    lamV=lamDef(4:6,:); nlamV=sqrt(sum(lamV.^2,1)); lamM=lamDef(7,:);
    burnI=(s(1:end-1)>0.5)&(s(2:end)>0.5); coastI=(s(1:end-1)<0.5)&(s(2:end)<0.5);
    swI=find(diff(double(s>0.5))~=0); swI=min(swI,N); swI(swI<1)=[];
    W=nlamV.*c./mMid - lamM;
    beta=sum(W(swI))/sum(W(swI).^2);
    Sfun=1-beta*W;
    burnPct=100*mean(Sfun(burnI)<0);
    spreadPct=100*std(1./W(swI))/abs(mean(1./W(swI)));

    % H1 vs H2: violation locality
    vio=find((burnI & Sfun>=0) | (coastI & Sfun<=0));
    if isempty(swI), dsw=inf(size(vio)); else
        dsw=arrayfun(@(q) min(abs(q-swI)), vio); end
    nearPct=100*mean(dsw<=1);                      % within 1 interval of a switch
    % per-arc violation fraction (arcs delimited by switches)
    arcId=cumsum([1 double(diff(double(s(1:N)>0.5))~=0)]);   % 1 x N interval arc ids
    nArcs=max(arcId); flip=0; arcFrac=zeros(1,nArcs);
    for a=1:nArcs
        ivs=find(arcId==a); ivs=ivs(burnI(ivs)|coastI(ivs));
        if isempty(ivs), continue; end
        arcFrac(a)=mean((burnI(ivs)&Sfun(ivs)>=0)|(coastI(ivs)&Sfun(ivs)<=0));
        flip=flip+(arcFrac(a)>0.5);
    end
    % H3: drift of the implied scale along the trajectory
    invW=1./W(swI);
    if numel(swI)>2, drift=corr(swI(:), invW(:)); else, drift=NaN; end

    diag(e)=struct('name',cases(e).name,'burnPct',burnPct,'violNearSwPct',nearPct, ...
        'violMidArcPct',100-nearPct,'nArcs',nArcs,'nFlippedArcs',flip, ...
        'betaDriftCorr',drift,'spreadPct',spreadPct);
    fprintf('%-16s burn%%=%5.1f | viol: %4.1f%% near-switch, %4.1f%% mid-arc | flipped arcs %d/%d | 1/W drift corr %+.2f | beta-spread %.1f%%\n', ...
        cases(e).name, burnPct, nearPct, 100-nearPct, flip, nArcs, drift, spreadPct);

    % panels: left = sign trace + violations; right = implied scale per switch
    tt=X(8,1:N)*p.tStar/86400;
    subplot(numel(cases),2,2*e-1); hold on; grid on; box on;
    plot(tt, sign(Sfun).*min(abs(Sfun),1), '-', 'Color',[0.4 0.4 0.7]);
    sb=nan(1,N); sb(burnI)=-1.05; plot(tt,sb,'.','Color',[0.85 0.4 0.1],'MarkerSize',4);
    if ~isempty(vio), plot(tt(vio), zeros(size(vio)), '.', 'Color',[0.8 0 0], 'MarkerSize',6); end
    ylim([-1.2 1.2]); ylabel('sgn(S)\cdotmin(|S|,1)');
    title(sprintf('%s: S recovered from duals (red = sign violations, orange = burn arcs)', cases(e).name), 'FontSize',9);
    if e==numel(cases), xlabel('t (days)'); end
    subplot(numel(cases),2,2*e); hold on; grid on; box on;
    plot(1:numel(swI), invW/abs(mean(invW)), 'o-', 'Color',[0.2 0.5 0.2], 'MarkerSize',4);
    yline(1,'--'); ylabel('(1/W_k)/mean'); title('implied scale per switch','FontSize',9);
    if e==numel(cases), xlabel('switch index'); end
end

outP=fullfile(cfg.dirs.plots,'beta_checker_diag.png');
exportgraphics(fig,outP,'Resolution',150); close(fig);
fprintf('WROTE %s\n', outP);
end

% ---------------------------------------------------------------------------
function q = load_case(name, r)
% Normalize a stored solution (front row or `out` struct) into a case record.
q = struct('name',name,'X',r.X,'U',r.U,'lamDef',r.lamDef);
end
