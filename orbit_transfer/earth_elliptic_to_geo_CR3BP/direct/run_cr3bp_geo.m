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
%                              U, defectDuals (discrete defect-constraint duals
%                              from opti.dual -- adjoint-proportional up to
%                              collocation scaling; see
%                              E2B/verify/hamiltonian_along_traj.m), sigma/L
%                              mesh, dL, time vector, throttle, m_f, dV, fp,
%                              provenance
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
thrustN   = 5;          % max thrust [N] (certified t_f anchors exist for
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
%runName   = 'cr3bp_T10N_phi0_fuel';   % basename for ALL artifacts of this run
runName   = 'cr3bp_T5N_phi0_fuel';   % basename for ALL artifacts of this run
movieMode = true;        % true -> render <runName>_movie.mp4/.gif (adds ~2 min)
rerun     = false;       % true -> ignore checkpoints, solve cold
maxIter   = 1500;        % IPOPT cap per solve
ipoptExtra= struct();    % opt-in IPOPT overrides (e.g. mumps_mem_percent) merged
                         %   over solver defaults; struct() = unchanged
liftDL    = false;       % true -> per-node lifted DeltaL (block-banded KKT; the
                         %   solver's designed remedy for the arrowhead-column
                         %   MUMPS factorization failure at large N; numerically
                         %   identical formulation)
% Programmatic override hook (ladder driver): if the caller's workspace
% defines LADDER_OVERRIDES (struct), its fields replace the defaults above.
% Absent -> byte-identical interactive behavior.
if exist('LADDER_OVERRIDES','var') && isstruct(LADDER_OVERRIDES)
    ovf = fieldnames(LADDER_OVERRIDES);
    for ovk = 1:numel(ovf), eval([ovf{ovk} ' = LADDER_OVERRIDES.(ovf{ovk});']); end %#ok<EVLDOT>
    clear ovf ovk
end

%% ------------------------------------------------------------------------
%% 2. SOLVE  (seed -> 2-body energy -> gain walk -> eps sharpen)
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;

% --- A7: parameter validation (fail loud before any solve work) ----------
assert(isscalar(thrustN) && isnumeric(thrustN) && isreal(thrustN) && ...
    isfinite(thrustN) && thrustN > 0, 'run_cr3bp_geo:badParam', ...
    'thrustN must be a finite positive real scalar (got %s)', mat2str(thrustN));
assert(isscalar(gain) && isnumeric(gain) && isreal(gain) && isfinite(gain) && ...
    gain >= 0 && gain <= 1, 'run_cr3bp_geo:badParam', ...
    'gain must be a real scalar in [0,1] (got %s)', mat2str(gain));
assert(isscalar(epsMin) && isnumeric(epsMin) && isreal(epsMin) && ...
    isfinite(epsMin) && epsMin >= 0 && epsMin <= 1, 'run_cr3bp_geo:badParam', ...
    'epsMin must be a real scalar in [0,1] (got %s)', mat2str(epsMin));
assert(isscalar(phi0) && isnumeric(phi0) && isreal(phi0) && isfinite(phi0), ...
    'run_cr3bp_geo:badParam', 'phi0 must be a real finite scalar (got %s)', mat2str(phi0));
assert(isscalar(ctf) && isnumeric(ctf) && isreal(ctf) && isfinite(ctf) && ctf > 0, ...
    'run_cr3bp_geo:badParam', 'ctf must be a finite positive real scalar (got %s)', mat2str(ctf));
assert(isempty(x0Elems) || (isnumeric(x0Elems) && numel(x0Elems)==5), ...
    'run_cr3bp_geo:badParam', 'x0Elems must be [] or a 5-element numeric vector (1x5 or 5x1)');
assert(isempty(xfElems) || (isnumeric(xfElems) && numel(xfElems)==5), ...
    'run_cr3bp_geo:badParam', 'xfElems must be [] or a 5-element numeric vector (1x5 or 5x1)');
assert(ischar(runName) && ~isempty(runName), ...
    'run_cr3bp_geo:badParam', 'runName must be a nonempty char');
if ~isempty(x0Elems), x0Elems = reshape(x0Elems(:), 5, 1); end   % 1x5 ok -> 5x1

x0IsDefault = isempty(x0Elems);   % A8: "default endpoints" recorded BEFORE
xfIsDefault = isempty(xfElems);   %     the [] -> literal-default overwrite below

resDir = fullfile(here, 'results');  if ~exist(resDir,'dir'), mkdir(resDir); end
ckDir  = fullfile(resDir, 'frontdoor');  if ~exist(ckDir,'dir'), mkdir(ckDir); end

par  = kepler_lt_params(thrustN, 1500, 2000);
pert = lunar_params(par, phi0, gain);
if isempty(xfElems), xfElems = [1;0;0;0;0]; end
xfElems = reshape(xfElems(:), 5, 1);   % 1x5 ok -> 5x1
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
catch ME
    if ~strcmp(ME.identifier, 'table3_recipes:unknownThrust'), rethrow(ME); end   % A6
end
x0ElemsFp = 'HMG-default';  if ~x0IsDefault, x0ElemsFp = x0Elems(:).'; end   % A1
fp = struct('thrustN',thrustN,'m0kg',par.m0kg,'ispS',par.ispS, ...
    'tfTarget',tfTargetTU,'muM',pert.muM,'DM',pert.DM,'nM',pert.nM, ...
    'phi0',phi0,'gain',gain,'epsMin',epsMin,'xf',xfElems(:).', ...
    'seedThr',seedThr,'nodesPerRev',nodesPerRev, ...
    'x0Elems',x0ElemsFp,'maxIter',maxIter);
fprintf('=== RUN_CR3BP_GEO [%s]: T=%g N, gain=%g, phi0=%g, epsMin=%g, tf=%.4f TU (%.2f d) ===\n', ...
    runName, thrustN, gain, phi0, epsMin, tfTargetTU, tfTargetTU*par.TU_s/86400);

% --- stage A: two-pass seed (probe revs, then full density) ------------------
seedOpts = struct('thr', seedThr, 'betaMode', 'tangential', 'N', 50, 'stopP', xfElems(1));
if ~isempty(x0Elems), seedOpts.initElems = x0Elems; end
[~,~,~,~, infoP] = mee_seed(par, seedOpts);   % probe always runs (A4 needs it on every call)

% A4: seed-rev sanity vs the certified 2-body table (off-config lunar physics
% may legitimately shift revs, so this WARNS -- it does not assert/block).
certRev = try_certified(thrustN);
if ~isempty(certRev)
    revErr = abs(infoP.nRev - certRev.revs) / certRev.revs;
    if revErr > 0.15
        warning('run_cr3bp_geo:revsOffCertified', ['seed probe nRev=%.3f differs ' ...
            'from the certified 2-body revs=%.3f by %.1f%% (>15%%) at thrustN=%g N ' ...
            '-- check seedThr/config'], infoP.nRev, certRev.revs, 100*revErr, thrustN);
    end
end

ckSeed = fullfile(ckDir, [runName '_seed.mat']);
if ~rerun && isfile(ckSeed)
    S = load(ckSeed);  check_fp_local(S, fp, ckSeed);
    sigma = S.sigma;  X0 = S.X0;  U0 = S.U0;  dL0 = S.dL0;
else
    % B2: on-table thrusts use the recipe's node count EXACTLY (recipe
    % fidelity); the max(60,...) floor is kept only for the off-table
    % fallback path, where nodesPerRev/seedThr are just campaign defaults,
    % not a certified recipe, so a thin low-rev mesh still needs a floor.
    if ~isempty(certRev)
        seedOpts.N = round(nodesPerRev * infoP.nRev);
    else
        seedOpts.N = max(60, round(nodesPerRev * infoP.nRev));
    end
    [sigma, X0, U0, dL0] = mee_seed(par, seedOpts);
    save(ckSeed, 'sigma','X0','U0','dL0','fp');
end
fp.N = size(X0,2) - 1;                        % A1: seed mesh N, known only now
S = load(ckSeed);  S.fp = fp;  save(ckSeed, '-struct', 'S');   % re-save with fp.N

% --- stage B: two-body energy solve (eps=1, no Moon) -------------------------
solveOpts = struct('par',par,'mode','fixedtf','eps',1,'tfTarget',tfTargetTU, ...
    'x0',X0(:,1),'xf',xfElems,'maxIter',maxIter,'warmTight',false, ...
    'ipoptExtra',ipoptExtra,'liftDL',liftDL);
ckE = fullfile(ckDir, [runName '_energy2b.mat']);
if ~rerun && isfile(ckE)
    S = load(ckE);  check_fp_local(S, fp, ckE);  o = S.o;
    % B3: guard against a partial cache delete recreating a different-N seed
    % (fp-mismatch alone would not catch a seed regenerated with the SAME
    % fp.N value but a differently-shaped o.X, e.g. a stale ckE surviving a
    % seed-only rerun).
    assert(size(S.o.X,2) == numel(sigma), 'run_cr3bp_geo:meshMismatch', ...
        'checkpoint %s: cached mesh size (%d nodes) does not match the current seed mesh (%d nodes) -- set rerun=true', ...
        ckE, size(S.o.X,2), numel(sigma));
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
            assert(size(S.o.X,2) == numel(sigma), 'run_cr3bp_geo:meshMismatch', ...   % B3
                'checkpoint %s: cached mesh size (%d nodes) does not match the current seed mesh (%d nodes) -- set rerun=true', ...
                ckG, size(S.o.X,2), numel(sigma));
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
        'maxIter',maxIter,'sched',sched,'resDir',ckDir,'tag',[runName '_eps'],'fp',fp, ...
        'ipoptExtra',ipoptExtra,'liftDL',liftDL, ...
        'fpStrict',true);   % A2: E3B front door opts in to fail-closed caches
    [best, tbl] = homotopy_mee(sigma, Xk, Uk, dLk, ho); %#ok<ASGLU>
    assert(best.certified && abs(best.epsReached - epsMin) < 1e-12, ...
        'run_cr3bp_geo:uncertified', 'homotopy stalled at eps=%.4g (wanted %.4g)', ...
        best.epsReached, epsMin);
end
% A5: a boundSaturation warning quietly means a box bound, not the physics,
% is the real blocker (defect/termErr alone can both read machine precision
% anyway) -- check right after the final solve.
[wmsg, wid] = lastwarn();
if strcmp(wid, 'casadi_lt_mee:boundSaturation')
    error('run_cr3bp_geo:boundSaturation', ['casadi_lt_mee:boundSaturation fired ' ...
        'during the final solve: %s -- STOP, do not treat this solution as ' ...
        'certified'], wmsg);
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
mf_kg = par.m0kg*best.mf;

% A8: basin-reference warning -- catches the "silent basin failure" class
% (e.g. the 1.84 kg regression) at the canonical certified configuration only.
certRef = try_certified(thrustN);
if ~isempty(certRef) && gain==1 && epsMin==0 && ctf==1.5 && x0IsDefault && xfIsDefault
    if mf_kg < certRef.m_f_kg - 0.1
        warning('run_cr3bp_geo:wrongBasin', ['final mass m_f=%.4f kg is more than ' ...
            '0.1 kg below the certified 2-body reference m_f=%.4f kg at thrustN=%g N ' ...
            '-- possible silent basin failure'], mf_kg, certRef.m_f_kg, thrustN);
    end
end

products = struct('X',X,'U',U,'defectDuals',best.lamDef,'sigma',sigma,'L',Lmesh, ...
    'dL',dL,'t_TU',tTU,'t_days',tDays,'throttle',U(4,:), ...
    'm_f_kg',mf_kg,'dV_kms',best.dV_kms,'switches',best.switches, ...
    'fp',fp,'provenance',prov); %#ok<NASGU>
outFile = fullfile(resDir, [runName '.mat']);
save(outFile, '-struct', 'products');
fprintf('DATA PRODUCTS: %s\n  m_f=%.4f kg  dV=%.4f km/s  switches=%d  defect=%.2e\n', ...
    outFile, mf_kg, best.dV_kms, best.switches, best.maxDefect);

%% ------------------------------------------------------------------------
%% 4. PLOTS  (2D top-down + 3D burn/coast; throttle + mass vs time)
%% ------------------------------------------------------------------------
% A9: ctfEff is the ACTUAL tf/tfmin ratio realized (NaN off-table -- the raw
% user ctf is not a meaningful ratio there since tfTargetTU may have been
% given explicitly rather than derived from ctf*cert.tfmin).
certCtf = try_certified(thrustN);
if ~isempty(certCtf), ctfEff = tfTargetTU / certCtf.tfmin; else, ctfEff = NaN; end
cart = mee_res_to_cart_res(X, U, dL, sigma, thrustN, ctfEff, 1, 8);
r = cart.fuel.X(1:3,:);  sTh = cart.fuel.U(4,:);   % B4: unused tD extraction removed (plots use tDays, from the MEE state's own time row)
burn = sTh > 0.05;

% A11: lunar-approach envelope check (informational -- the lt_mee_rhs d3
% guard is only VALIDATED for spacecraft-Moon separation >= ~8 LU; this is
% NOT a collision model). Only meaningful when the lunar term is active.
if gain > 0
    angM     = pert.nM*cart.fuel.X(8,:) + pert.phi0;
    rMoon    = [pert.DM*cos(angM); pert.DM*sin(angM); zeros(size(angM))];
    minSepLU = min(vecnorm(cart.fuel.X(1:3,:) - rMoon, 2, 1));
    if minSepLU < 2
        warning('run_cr3bp_geo:moonApproach', ['spacecraft-Moon separation drops ' ...
            'to %.3f LU (<2 LU) along this trajectory -- below the lt_mee_rhs d3 ' ...
            'guard''s validated envelope (>=~8 LU); this result is not certified ' ...
            'in that close-approach regime'], minSepLU);
    end
end
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
    if isnan(ctfEff), lbl = [lbl ' (custom t_f)']; end   % A9
    cart.cfg.label = lbl;
    transfer_movie(cart, fullfile(resDir, [runName '_movie']));
end
fprintf('RUN_CR3BP_GEO done [%s]\n', runName);

%% ------------------------------------------------------------------------
function assert_gate(o, stage)
% ASSERT_GATE  Certification gate: full convergence + tight defect/unit-norm/
% terminal-error (A5: maxDefect alone can read machine precision while the
% control-cone unit norm or the terminal-target miss is not actually tight).
assert(strcmp(o.ipoptStatus,'Solve_Succeeded') && o.maxDefect < 1e-6 && ...
    o.maxUnit < 1e-8 && o.termErr < 1e-8, 'run_cr3bp_geo:gate', ...
    '%s stage not certified (status=%s, defect=%.2e, maxUnit=%.2e, termErr=%.2e)', ...
    stage, o.ipoptStatus, o.maxDefect, o.maxUnit, o.termErr);
end

function cert = try_certified(thrustN)
% TRY_CERTIFIED  table3_certified(thrustN), or [] if thrustN is off-table
% (narrow catch on the registry's OWN unknown-thrust error id -- any other
% error rethrows, per A6's narrow-catch discipline).
try
    cert = table3_certified(thrustN);
catch ME
    if ~strcmp(ME.identifier, 'table3_certified:unknownThrust'), rethrow(ME); end
    cert = [];
end
end

function check_fp_local(S, fp, file)
% CHECK_FP_LOCAL  Fail-loud fingerprint guard on front-door checkpoints.
% B1 (2nd-reviewer wave): these checkpoints are all campaign-fresh (this
% front door has no legacy user base to warn-and-trust) -- a field MISSING
% from the cached fingerprint is now a hard error too, not just a value
% mismatch.
if ~isfield(S,'fp'), error('run_cr3bp_geo:noFp', '%s has no fingerprint', file); end
fn = fieldnames(fp);
for kf = 1:numel(fn)
    f = fn{kf};
    if ~isfield(S.fp,f)
        error('run_cr3bp_geo:fpSchemaOlder', ...
            'checkpoint %s: field ''%s'' is absent from the cached fingerprint -- change runName or set rerun=true', file, f);
    end
    if ~isequal(S.fp.(f), fp.(f))
        error('run_cr3bp_geo:fpMismatch', ...
            'checkpoint %s: field ''%s'' differs from current parameters -- change runName or set rerun=true', file, f);
    end
end
end
