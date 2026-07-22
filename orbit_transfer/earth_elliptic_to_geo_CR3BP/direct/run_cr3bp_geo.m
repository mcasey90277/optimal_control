% RUN_CR3BP_GEO  Front door: elliptic->GEO low-thrust transfer WITH lunar gravity.
%
% One self-contained entry script (house run_gergaud pattern): set the
% parameters in section 1, run, and get a solved transfer + saved data
% products + plots + (optionally) a movie, all under a user-chosen run name.
% The chain is the certified Phase-1 pipeline: two-pass seed -> two-body
% energy solve -> mu-continuation (lunar gain 0->1) -> Bertrand-Epenoy
% eps-homotopy down to epsMin.
%
%   epsMin = 1  -> minimum ENERGY solution (smooth throttle; no sharpening)
%   epsMin = 0  -> minimum FUEL solution (bang-bang; the campaign objective)
%   0<epsMin<1  -> eps-optimal (partially sharpened) solution
%
% Every stage is gated on Solve_Succeeded + defect < 1e-6 and checkpointed
% under the run name (re-invoking resumes; set rerun=true for a cold start).
%
% OUTPUTS (files, under results/):
%   <runName>.mat            - data products: states X (incl. t row), controls
%                              U, discrete costates lamDef, sigma/L mesh, dL,
%                              time vector, throttle, m_f, dV, fp, provenance
%   <runName>_traj.png       - 2D top-down + 3D trajectory (red burn/blue coast)
%   <runName>_throttle.png   - throttle + mass vs time
%   <runName>_movie.mp4/.gif - (movieMode true) house-style transfer movie
%
% REFERENCES:
%   [1] doc/cr3bp_geo_phase1_note.tex (the OCP + pipeline this script runs).
%   [2] ../../earth_elliptic_to_geo/direct/frontdoor/run_gergaud.m (pattern).

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------
thrustN   = 10;          % max thrust [N] (certified t_f anchors exist for
                         %   10/5/2.5/1/0.5/0.2/0.1; other values need tfTargetTU)
phi0      = 0;           % Moon phase at t=0 [rad]: 0 = Moon on +x (the perigee
                         %   direction; spacecraft departs apogee on -x side)
gain      = 1;           % lunar mass scale: 1 = full Moon (CR3BP), 0 = two-body
epsMin    = 0;           % homotopy endpoint: 1 = min-ENERGY, 0 = min-FUEL
ctf       = 1.5;         % t_f = ctf * tfMin(thrustN) (two-body certified anchor)
tfTargetTU= [];          % [] = auto from certified table; or explicit t_f [TU]
x0Elems   = [];          % [] = HMG GTO start [P0;ex0;ey0;hx0;hy0] (P0=11625 km,
                         %   e=0.75 along +x, i=7 deg); or user [5x1] MEE
xfElems   = [];          % [] = GEO target [1;0;0;0;0]; or user [5x1] MEE
runName   = 'cr3bp_T10N_phi0_fuel';   % basename for ALL artifacts of this run
movieMode = true;        % true -> render <runName>_movie.mp4/.gif (adds ~2 min)
rerun     = false;       % true -> ignore checkpoints, solve cold
maxIter   = 1500;        % IPOPT cap per solve

%% ------------------------------------------------------------------------
%% 2. SOLVE  (seed -> 2-body energy -> gain walk -> eps sharpen)
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
resDir = fullfile(here, 'results');  if ~exist(resDir,'dir'), mkdir(resDir); end
ckDir  = fullfile(resDir, 'frontdoor');  if ~exist(ckDir,'dir'), mkdir(ckDir); end

par  = kepler_lt_params(thrustN, 1500, 2000);
pert = lunar_params(par, phi0, gain);
if isempty(xfElems), xfElems = [1;0;0;0;0]; end
if isempty(tfTargetTU)
    cert = table3_certified(thrustN);        % errors informatively off-table
    tfTargetTU = ctf * cert.tfmin;
end
% Seed knobs from the certified per-rung recipe registry (seedThr=0.999 was
% the front door's original sin: a near-full-throttle seed spirals in fewer
% revs and lands a WORSE basin -- 1375.31/17sw vs certified 1377.15/19sw at
% 10 N). Off-table thrusts fall back to the campaign defaults (0.4 / 25).
seedThr = 0.4;  nodesPerRev = 25;
try
    rec = table3_recipes(thrustN);
    seedThr = rec.fuel.seedThr;  nodesPerRev = rec.fuel.npr;
catch
end
fp = struct('thrustN',thrustN,'m0kg',par.m0kg,'ispS',par.ispS, ...
    'tfTarget',tfTargetTU,'muM',pert.muM,'DM',pert.DM,'nM',pert.nM, ...
    'phi0',phi0,'gain',gain,'epsMin',epsMin,'xf',xfElems(:).', ...
    'seedThr',seedThr,'nodesPerRev',nodesPerRev);
fprintf('=== RUN_CR3BP_GEO [%s]: T=%g N, gain=%g, phi0=%g, epsMin=%g, tf=%.4f TU (%.2f d) ===\n', ...
    runName, thrustN, gain, phi0, epsMin, tfTargetTU, tfTargetTU*par.TU_s/86400);

% --- stage A: two-pass seed (probe revs, then full density) ------------------
seedOpts = struct('thr', seedThr, 'betaMode', 'tangential', 'N', 50, 'stopP', xfElems(1));
if ~isempty(x0Elems), seedOpts.initElems = x0Elems; end
ckSeed = fullfile(ckDir, [runName '_seed.mat']);
if ~rerun && isfile(ckSeed)
    S = load(ckSeed);  check_fp_local(S, fp, ckSeed);
    sigma = S.sigma;  X0 = S.X0;  U0 = S.U0;  dL0 = S.dL0;
else
    [~,~,~,~, infoP] = mee_seed(par, seedOpts);
    seedOpts.N = max(60, round(nodesPerRev * infoP.nRev));
    [sigma, X0, U0, dL0] = mee_seed(par, seedOpts);
    save(ckSeed, 'sigma','X0','U0','dL0','fp');
end

% --- stage B: two-body energy solve (eps=1, no Moon) -------------------------
solveOpts = struct('par',par,'mode','fixedtf','eps',1,'tfTarget',tfTargetTU, ...
    'x0',X0(:,1),'xf',xfElems,'maxIter',maxIter,'warmTight',false);
ckE = fullfile(ckDir, [runName '_energy2b.mat']);
if ~rerun && isfile(ckE)
    S = load(ckE);  check_fp_local(S, fp, ckE);  o = S.o;
else
    o = casadi_lt_mee(sigma, X0, U0, dL0, solveOpts);
    assert_gate(o, 'two-body energy');
    save(ckE, 'o','fp');
end
Xk = o.X;  Uk = o.U;  dLk = o.dL;

% --- stage C: mu-continuation (lunar gain walk on the energy objective) ------
if gain > 0
    gsched = [0.25 0.5 0.75 1.0] * gain;
    solveOpts.warmTight = true;
    for gk = gsched
        parG = par;  parG.pert = pert;  parG.pert.gain = gk;
        solveOpts.par = parG;
        ckG = fullfile(ckDir, sprintf('%s_gain%04d.mat', runName, round(1000*gk)));
        if ~rerun && isfile(ckG)
            S = load(ckG);  check_fp_local(S, fp, ckG);  o = S.o;
        else
            o = casadi_lt_mee(sigma, Xk, Uk, dLk, solveOpts);
            assert_gate(o, sprintf('gain=%.3f', gk));
            save(ckG, 'o','fp');
        end
        Xk = o.X;  Uk = o.U;  dLk = o.dL;
        fprintf('  gain=%.3f OK  defect=%.2e  m_f=%.4f kg\n', gk, o.maxDefect, par.m0kg*o.mf);
    end
end

% --- stage D: eps-homotopy down to epsMin at the final physics ---------------
best = o;
if epsMin < 1
    parG = par;  if gain > 0, parG.pert = pert; end
    schedFull = [0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0];
    sched = [schedFull(schedFull > epsMin), epsMin];
    if rerun   % homotopy_mee has its OWN per-step resume caches under this tag;
               % a cold rerun must clear them too, else stale steps from a prior
               % configuration are warn-and-trusted (schema-older fp path)
        delete(fullfile(ckDir, [runName '_eps_step*.mat']));
    end
    ho = struct('par',parG,'x0',Xk(:,1),'xf',xfElems,'tfTarget',tfTargetTU, ...
        'maxIter',maxIter,'sched',sched,'resDir',ckDir,'tag',[runName '_eps'],'fp',fp);
    [best, tbl] = homotopy_mee(sigma, Xk, Uk, dLk, ho); %#ok<ASGLU>
    assert(best.certified && abs(best.epsReached - epsMin) < 1e-12, ...
        'run_cr3bp_geo:uncertified', 'homotopy stalled at eps=%.4g (wanted %.4g)', ...
        best.epsReached, epsMin);
end
assert_gate(best, 'final');
X = best.X;  U = best.U;  dL = best.dL;

%% ------------------------------------------------------------------------
%% 3. DATA PRODUCTS  (saved under the user-defined run name)
%% ------------------------------------------------------------------------
tTU   = X(7,:);  tDays = tTU * par.TU_s/86400;
Lmesh = pi + sigma(:).' * dL;
prov  = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
    'script','run_cr3bp_geo','pipeline','seed->energy2b->gainwalk->epsHomotopy', ...
    'ipoptStatus',best.ipoptStatus,'maxDefect',best.maxDefect, ...
    'epsReached', epsMin*(epsMin<1) + 1*(epsMin>=1));
products = struct('X',X,'U',U,'lamDef',best.lamDef,'sigma',sigma,'L',Lmesh, ...
    'dL',dL,'t_TU',tTU,'t_days',tDays,'throttle',U(4,:), ...
    'm_f_kg',par.m0kg*best.mf,'dV_kms',best.dV_kms,'switches',best.switches, ...
    'fp',fp,'provenance',prov); %#ok<NASGU>
outFile = fullfile(resDir, [runName '.mat']);
save(outFile, '-struct', 'products');
fprintf('DATA PRODUCTS: %s\n  m_f=%.4f kg  dV=%.4f km/s  switches=%d  defect=%.2e\n', ...
    outFile, par.m0kg*best.mf, best.dV_kms, best.switches, best.maxDefect);

%% ------------------------------------------------------------------------
%% 4. PLOTS  (2D top-down + 3D burn/coast; throttle + mass vs time)
%% ------------------------------------------------------------------------
cart = mee_res_to_cart_res(X, U, dL, sigma, thrustN, ctf, 1, 8);
r = cart.fuel.X(1:3,:);  sTh = cart.fuel.U(4,:);  tD = cart.fuel.X(8,:)*par.TU_s/86400;
burn = sTh > 0.05;
fig = figure('Visible','off','Position',[60 60 1280 560]);
thG = linspace(0,2*pi,361);
for sp = 1:2
    ax = subplot(1,2,sp);  hold(ax,'on');
    plot3(ax, cos(thG), sin(thG), 0*thG, '-', 'Color',[.2 .6 .2]);
    seg = @(msk,col) plot3(ax, r(1,msk), r(2,msk), r(3,msk), '.', ...
                           'Color',col, 'MarkerSize',4);
    seg(burn,  [0.85 0.2 0.15]);   % red  = thrusting
    seg(~burn, [0.2 0.35 0.8]);    % blue = coasting
    plot3(ax, 0,0,0,'o','MarkerFaceColor',[.2 .4 .9],'MarkerSize',9);
    quiver3(ax, 0,0,0, 1.3*cos(phi0), 1.3*sin(phi0), 0, 'k','LineWidth',1.1);
    text(ax, 1.35*cos(phi0), 1.35*sin(phi0), 0, 'to Moon (t=0)');
    axis(ax,'equal'); grid(ax,'on');
    xlabel(ax,'x (ND)'); ylabel(ax,'y (ND)');
    if sp==1, view(ax,2); title(ax,'top-down'); else, view(ax,3); zlabel(ax,'z'); title(ax,'3D'); end
end
sgtitle(fig, sprintf('%s:  T=%g N, gain=%g, \\phi_0=%.2f, \\epsilon_{min}=%g  (red burn / blue coast)', ...
    strrep(runName,'_','\_'), thrustN, gain, phi0, epsMin));
exportgraphics(fig, fullfile(resDir,[runName '_traj.png']), 'Resolution', 150);
fig2 = figure('Visible','off','Position',[60 60 1100 500]);
subplot(2,1,1); stairs(tDays, U(4,:), 'k','LineWidth',1); ylim([-0.05 1.05]); grid on
ylabel('throttle \delta'); title(sprintf('thrust profile (%d switches)', best.switches));
subplot(2,1,2); plot(tDays, par.m0kg*X(6,:), 'r','LineWidth',1.2); grid on
xlabel('time (days)'); ylabel('mass (kg)');
exportgraphics(fig2, fullfile(resDir,[runName '_throttle.png']), 'Resolution', 150);
close(fig); close(fig2);
fprintf('PLOTS: %s_traj.png, %s_throttle.png\n', runName, runName);

%% ------------------------------------------------------------------------
%% 5. MOVIE  (optional)
%% ------------------------------------------------------------------------
if movieMode
    addpath(fullfile(here, '..','..','earth_elliptic_to_geo','direct','viz'));
    lbl = 'min-fuel';  if epsMin >= 1, lbl = 'min-energy'; elseif epsMin > 0, lbl = sprintf('\\epsilon=%g-optimal', epsMin); end
    if gain > 0, lbl = [lbl ' (CR3BP)']; else, lbl = [lbl ' (2-body)']; end
    cart.cfg.label = lbl;
    transfer_movie(cart, fullfile(resDir, [runName '_movie']));
end
fprintf('RUN_CR3BP_GEO done [%s]\n', runName);

%% ------------------------------------------------------------------------
function assert_gate(o, stage)
% ASSERT_GATE  Certification gate: full convergence + tight defect.
assert(strcmp(o.ipoptStatus,'Solve_Succeeded') && o.maxDefect < 1e-6, ...
    'run_cr3bp_geo:gate', '%s stage not certified (status=%s, defect=%.2e)', ...
    stage, o.ipoptStatus, o.maxDefect);
end

function check_fp_local(S, fp, file)
% CHECK_FP_LOCAL  Fail-loud fingerprint guard on front-door checkpoints.
if ~isfield(S,'fp'), error('run_cr3bp_geo:noFp', '%s has no fingerprint', file); end
fn = fieldnames(fp);
for kf = 1:numel(fn)
    f = fn{kf};
    if isfield(S.fp,f) && ~isequal(S.fp.(f), fp.(f))
        error('run_cr3bp_geo:fpMismatch', ...
            'checkpoint %s: field ''%s'' differs from current parameters -- change runName or set rerun=true', file, f);
    end
end
end
