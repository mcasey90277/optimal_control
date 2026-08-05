function S = sweep_phasing_direct(opts)
% SWEEP_PHASING_DIRECT  Phasing map by DIRECT solves; costates via the covector map.
%
% THE PIVOT, and why. Seven variants of costate-continuation marching
% (indirect/sweep_phasing.m) could not cross a 1/6-period grid edge from the
% anchor: initial-state perturbations amplify ~1e3x through the flow, which is
% intrinsic to single shooting -- no predictor, Jacobian, or budget above it
% cured the fragility. The DIRECT transcription does not have that failure
% mode: its constraints are local, a trajectory-space warm start degrades
% benignly, and this campaign's certified results (reference reproduced to 7
% digits; the 3.817 discovery) were all produced this way. So the map is built
% by direct solves, one per phasing pair, and the costate catalog Darin wants
% falls out of the VALIDATED covector mapping (defect multipliers = costates to
% 6e-4 deg, scale 0.99999 -- see certify/costate_compare.m). pumpkyn's tfMin
% can then be seeded per point as a zero-distance polish where an indirect twin
% is wanted; that is optional and off the critical path.
%
% PER GRID POINT (fD, fA):
%   1. endpoints from the propagated DRO / tulip at those phase fractions;
%   2. direct solve: Hermite-Simpson + Sundman, N intervals, warm-started from
%      the nearest already-solved neighbour (BFS from the anchor), resampled to
%      uniform time so the solver's Sundman re-distribution applies cleanly;
%   3. gate: global single-shot accuracy (position km / velocity m/s), NLP
%      defect, throttle floor. (The full per-interval residual is reserved for
%      the map's interesting points afterwards -- it costs more than the solve.)
%   4. record: t_f, mass, dV, periselene, ARRIVAL-POINT altitude (the ceiling
%      any altitude policy has at that phasing -- see FINDINGS), and the
%      sign-resolved initial costates lambda0 from the defect multipliers:
%      THE CATALOG ENTRY.
%
% INPUTS:
%   opts - struct (optional):
%          .nD, .nA   grid sizes                        [default 6, 6]
%          .N         collocation intervals             [default 800]
%          .thrustN / .ispS / .m0kg                     [0.07 / 900 / 150]
%          .globTolKm / .globTolMs   gate tolerances    [1.0 km / 1.0 m/s]
%          .maxIter   IPOPT cap per solve               [default 3000]
%          .verbose                                     [default true]
%
% OUTPUTS:
%   S - struct: .sD .sA (phase grids), and per-point [nD x nA]:
%       .TF .MF .DV .PERIS .ARRALT .GLOBKM .GLOBMS .DEFECT .PASS .WALL
%       .LAM0 [8 x nD x nA] initial costates + tf (the catalog; NaN if the
%        point failed), .anchor, .meta
%       Progress log + checkpoints in results/ alongside the final .mat.
%
% REFERENCES:
%   [1] indirect/sweep_phasing.m -- the marching approach and its verdict.
%   [2] certify/costate_compare.m -- the covector mapping this harvests.
%   [3] DRO_tulip/FINDINGS.md -- the arrival-altitude ceiling channel.

if nargin < 1, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
nD = d('nD',6);  nA = d('nA',6);  N = d('N',800);
thrustN = d('thrustN',0.07);  ispS = d('ispS',900);  m0 = d('m0kg',150);
globTolKm = d('globTolKm',1.0);  globTolMs = d('globTolMs',1.0);
maxIter = d('maxIter',3000);  verbose = d('verbose',true);
floorKm = d('floorKm',500);          % closes the through-the-Moon hole; known
sA0off  = d('sA0',[]);               % not to distort any legitimate basin
sD0off  = d('sD0',0);
% independent = true: solve every node from the CRUDE INTERNAL seed, no
% chaining. Motivated by the transect finding that the solution family BREAKS
% at the tulip's fast segments -- continuation in any representation (costate
% or trajectory, any step size down to 1/160 period) fails to cross them, so
% chained sweeps cannot circle the torus. Independent cold solves also remove
% the warm-start confound from difficulty measurements, which is exactly what
% the arrival-difficulty hypothesis needs.
indep   = d('independent',false);
% ORBIT PARAMETERS AS OPTIONS -- future sweeps vary these (the period axis is
% Darin's third dimension), and every data product must record which orbits it
% belongs to. Defaults reproduce the original DRO(tau0=1) -> tulip(7,-1) pair.
% REFERENCE-GUIDED RE-SOLVE (Mike's shortcut): when a previous run's results
% are supplied, attempt ONLY the cells that run proved green -- failed attempts
% produce no reusable product, so re-running them buys nothing but audit
% completeness. refPass = logical [nD x nA] from the prior map.
refMat  = d('refMat', '');
% w0Cells: [k x 2] (iD,iA) list restricting WAVE 0 to known cold-openers (the
% prior run's [w0]-PASS cells). Cold attempts elsewhere are known to fail and
% only green via the wave-1 flood, so re-paying them buys nothing.
w0Cells = d('w0Cells', []);
tauDRO  = d('tauDRO', 1.0);
NpT     = d('NpTulip', 7);
tauT    = d('tauTulip', 5*2*pi/6);
pmT     = d('pmTulip', -1);

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'));  addpath(fullfile(here,'certify'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
rd = fullfile(here,'results');  if ~isfolder(rd), mkdir(rd); end
logF  = fullfile(rd, sprintf('dsweep_%dx%d_progress.txt', nD, nA));
ckptF = fullfile(rd, sprintf('dsweep_%dx%d_ckpt.mat', nD, nA));
ckptCellsF = fullfile(rd, sprintf('dsweep_%dx%d_cells_ckpt.mat', nD, nA));
log_ = @(varargin) local_filelog(logF, verbose, varargin{:});

muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
rMoonKm = 1737.4;
g0   = 9.80665*tStar^2/(1000*lStar);
Tmax = (thrustN/m0)*tStar^2/(lStar*1000);
c    = (ispS/tStar)*g0;
tW = tic;

%% orbits and anchor (same construction as the marching harness)
[~, rvD0] = pumpkynPie.cr3bp.getDRO(tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, NpT, pmT);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, muStar);
dep = @(f) interp1(tD, rvD, mod(f,1)*tD(end), 'spline');
arr = @(f) interp1(tT, rvT, mod(f,1)*tT(end), 'spline');
dvTheta = pumpkyn.util.bsxAng(rvT(:,4:6), rvD(1,4:6), 2);
[~, idxA] = max(dvTheta);
fA0 = tT(idxA)/tT(end);
if isempty(sA0off), sA0off = fA0; end
sD = mod(sD0off + (0:nD-1)/nD, 1);
sA = mod(sA0off + (0:nA-1)/nA, 1);

% anchor warm start: the indirect reference trajectory, resampled uniformly
lam = [13.9579470518969; 6.26766054139842; 7.02087277890196; -0.969930036390389; ...
       -5.04449989501292; -3.11628049459521; 5.50150618633038; 4.01524259262941];
sol = pumpkyn.cr3bp.tfMin(rvD(1,1:6), rvT(idxA,1:6), lam, Tmax, c, muStar);
[tauR, rvR] = pumpkyn.cr3bp.tfMinProp(sol(8), [rvD(1,1:6), 1, sol(1:7)'], Tmax, c, muStar);
sN = linspace(0, tauR(end), N+1);
seed0.X  = interp1(tauR, rvR(:,1:7), sN, 'spline').';
LV = interp1(tauR, rvR(:,11:13), sN, 'spline');
seed0.U  = [(-LV./max(vecnorm(LV,2,2),eps)).'; ones(1,N+1)];
seed0.tf = sol(8);

log_('DIRECT SWEEP %dx%d, N=%d, sundman HS  (anchor at sD=%.3f sA=%.3f)\n', ...
    nD, nA, N, sD(1), sA(1));

%% the wave engine
% WAVE 0 (Mike's triage): every grid point gets one CHEAP cold attempt -- the
% crude chord seed under a TIGHT iteration budget, so easy points convert in
% ~their usual few hundred iterations and grinders fail fast instead of
% burning an hour. The passes are the patch openers.
% WAVE 1: multi-source BFS -- every pass seeds its neighbours trajectory-wise
% (the measured-strong regime inside patches; the creases stop the flood, and
% where it stops IS the fold-line map). Up to maxTriesW attempts per point,
% one from each solved neighbour.
% Plots refresh at every checkpoint so the torus picture is always current.
TF = nan(nD,nA); MF = TF; DV = TF; PERIS = TF; ARRALT = TF;
GLOBKM = TF; GLOBMS = TF; DEFECT = TF; WALL = TF; ARRSPD = TF; PASS = false(nD,nA);
LAM0 = nan(8,nD,nA);  LAMT = nan(nD,nA);  CELLS = cell(nD,nA);
seeds = cell(nD,nA);
TRIES = zeros(nD,nA);
nbr = [1 0; -1 0; 0 1; 0 -1];
nAttempt = 0;
refPass = true(nD,nA);
if ~isempty(refMat)
    RM = load(refMat);
    if isfield(RM,'S'), refPass = RM.S.PASS; else, refPass = RM.PASS; end
    assert(isequal(size(refPass), [nD nA]), ...
        'refMat PASS is %dx%d, grid is %dx%d -- incompatible reference', ...
        size(refPass,1), size(refPass,2), nD, nA);
    log_('reference-guided: attempting only the %d cells the reference run proved green\n', nnz(refPass));
end
w0Mask = true(nD,nA);
if ~isempty(w0Cells)
    w0Mask = false(nD,nA);
    for kw = 1:size(w0Cells,1), w0Mask(w0Cells(kw,1),w0Cells(kw,2)) = true; end
    log_('wave 0 restricted to %d known cold-opener cells\n', size(w0Cells,1));
end
mode = lower(d('mode','waves'));
if indep, mode = 'independent'; end
wave0Iter = d('wave0Iter',500);
maxTriesW = d('maxTriesW',4);

    function okOut = try_node(iD, iA, seedIn, iterCap, tag)
        % one gated solve attempt at grid point (iD,iA); records everything.
        okOut = false;
        rv0 = dep(sD(iD));  rvf = arr(sA(iA));
        t0 = tic;
        try
            if isempty(seedIn)
                X0 = [];  U0 = [];  tf0 = [];
            else
                X0 = seedIn.X;  U0 = seedIn.U;  tf0 = seedIn.tf;
            end
            cpuCap = 600;
            o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tmax, c, muStar, N, ...
                    X0, U0, tf0, struct('maxIter',iterCap,'scheme','hermite-simpson', ...
                    'sundman',true,'returnModel',true,'minAltKm',floorKm, ...
                    'maxCpuSec',cpuCap));
            % NEVER FLY GARBAGE. The continuous verification integrates the
            % returned trajectory, and for an unconverged iterate those
            % integrations crawl near the singularity UNBOUNDED (measured: one
            % junk point held the sweep 38 minutes after its solver cap). If
            % the NLP itself did not converge cleanly, the point is already
            % fail on primal evidence -- skip verification outright.
            if o.success && o.maxDefect < 1e-9
                C = certify_dro_mintime(o, struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...
                        Tmax, c, struct('tfRef',[],'verbose',false,'posTolKm',inf));
            else
                C = struct('globKm', Inf, 'gates', struct('id','G1bv','value',Inf));
            end
            w = toc(t0);
            gOK = o.success && C.globKm < globTolKm && ...
                  C.gates(strcmp({C.gates.id},'G1bv')).value < globTolMs && o.maxDefect < 1e-9;
            % keep the BEST result per point: a pass always beats a fail; among
            % passes, smaller t_f wins (basin bookkeeping happens downstream)
            better = gOK && (~PASS(iD,iA) || o.tf < TF(iD,iA));
            if isnan(TF(iD,iA)) || better
                r2 = vecnorm(o.X(1:3,:) - [1-muStar;0;0], 2, 1);
                TF(iD,iA) = o.tf;  MF(iD,iA) = o.mf;
                DV(iD,iA) = c*lStar/tStar*log(1/o.mf);
                PERIS(iD,iA)  = min(r2)*lStar - rMoonKm;
                ARRALT(iD,iA) = norm(rvf(1:3) - [1-muStar 0 0])*lStar - rMoonKm;
                ARRSPD(iD,iA) = norm(rvf(4:6))*lStar/tStar;
                GLOBKM(iD,iA) = C.globKm;
                GLOBMS(iD,iA) = C.gates(strcmp({C.gates.id},'G1bv')).value;
                DEFECT(iD,iA) = o.maxDefect;  WALL(iD,iA) = w;
                PASS(iD,iA) = gOK;
                if ~isempty(o.lamDef)
                    lamD = o.lamDef(1:7,:);              % physical costates only
                    % sign resolution: majority vote of primer-vs-control
                    % alignment over the first nodes (single-node comparison
                    % is vulnerable to one noisy station -- review finding)
                    nv = min(10, size(lamD,2));  vote = 0;
                    for kv = 1:nv
                        aS = o.Um(1:3,kv)/max(norm(o.Um(1:3,kv)),eps);
                        pv = -lamD(4:6,kv)/max(norm(lamD(4:6,kv)),eps);
                        vote = vote + sign(pv.'*aS);
                    end
                    sg = 1;  if vote < 0, sg = -1; end
                    LAM0(:,iD,iA) = [sg*lamD(:,1); o.tf];
                    if size(o.lamDef,1) >= 8             % Sundman: lambda_t consistency
                        LAMT(iD,iA) = sg*median(o.lamDef(8,:));
                    end
                    % SELF-CONTAINED cell record: someone holding only this
                    % struct can reconstruct and re-verify the solution --
                    % states, controls, costates, the phases it connects, the
                    % actual endpoint states used, and the physics constants.
                    CELLS{iD,iA} = struct('X',o.X,'U',o.U,'Um',o.Um, ...
                        'tNodes',o.tNodes,'tf',o.tf,'lamDef',sg*o.lamDef, ...
                        'sD',sD(iD),'sA',sA(iA), ...
                        'rv0',rv0(1:6),'rvf',rvf(1:6), ...
                        'Tmax',Tmax,'c',c,'muStar',muStar, ...
                        'floorKm',floorKm,'N',N);
                end
                if gOK
                    tu = linspace(0, o.tf, N+1);
                    sd2.X = interp1(o.tNodes, o.X.', tu, 'spline').';
                    sd2.U = interp1(o.tNodes, o.U.', tu, 'spline').';
                    nr2 = max(vecnorm(sd2.U(1:3,:),2,1), eps);
                    sd2.U(1:3,:) = sd2.U(1:3,:)./nr2;
                    sd2.U(4,:) = min(max(sd2.U(4,:),0),1);
                    sd2.tf = o.tf;
                    seeds{iD,iA} = sd2;
                end
            end
            log_('  [%s] (%2d,%2d) tf=%.5f dV=%.4f peris=%6.0f vArr=%.3f glob=%.2fkm %s (%.0fs)\n', ...
                tag, iD, iA, o.tf, c*lStar/tStar*log(1/o.mf), ...
                min(vecnorm(o.X(1:3,:)-[1-muStar;0;0],2,1))*lStar-rMoonKm, ...
                norm(rvf(4:6))*lStar/tStar, C.globKm, local_pf(gOK), w);
            okOut = gOK;
        catch ME
            log_('  [%s] (%2d,%2d) ERROR after %.0fs: %s\n', tag, iD, iA, toc(t0), ME.message);
        end
        TRIES(iD,iA) = TRIES(iD,iA) + 1;
        nAttempt = nAttempt + 1;
        if true   % save after EVERY square (Mike's directive 2026-08-04):
                  % worst-case loss on any crash = the square in flight
            save(ckptF,'TF','MF','DV','PERIS','ARRALT','ARRSPD','GLOBKM','GLOBMS', ...
                 'DEFECT','PASS','WALL','LAM0','LAMT','sD','sA','TRIES');
            % THE CELLS ARE THE IRREPLACEABLE PRODUCT -- checkpoint them too.
            % End-of-run-only saving is how 11 hours of trajectories and duals
            % were lost once already; never again. (Every 9 attempts: the file
            % rewrite grows to ~60 MB late in a run, a few seconds each.)
            save(ckptCellsF, 'CELLS', '-v7.3');
            try
                plot_phase_torus(ckptF, fullfile(rd, sprintf('phase_torus_%dx%d', nD, nA)));
                close all
            catch, end
        end
    end

switch mode
case 'waves'
    % Wave 0's job is to OPEN each patch, not to solve every point: knock on a
    % coarse sublattice with the FULL budget (measured: cold greens need well
    % over 500 iterations -- a tight budget over all points finds none), and
    % let wave 1 flood the full grid from whatever opens.
    stride = d('wave0Stride',2);
    log_('WAVE 0: cold chord seed, full budget %d, sublattice stride %d (%d points)\n', ...
        maxIter, stride, numel(1:stride:nD)*numel(1:stride:nA));
    for iD = 1:stride:nD
        for iA = 1:stride:nA
            if refPass(iD,iA) && w0Mask(iD,iA)
                try_node(iD, iA, [], maxIter, 'w0');
            end
        end
    end
    log_('WAVE 0 done: %d/%d pass. WAVE 1: multi-source chaining, full budget %d\n', ...
        nnz(PASS), nD*nA, maxIter);
    queue = {};
    for iD = 1:nD
        for iA = 1:nA
            if PASS(iD,iA), queue{end+1} = [iD iA]; end %#ok<AGROW>
        end
    end
    while ~isempty(queue)
        q = queue{1};  queue(1) = [];
        for kk = 1:4
            jD = mod(q(1)-1+nbr(kk,1), nD)+1;  jA = mod(q(2)-1+nbr(kk,2), nA)+1;
            if refPass(jD,jA) && ~PASS(jD,jA) && TRIES(jD,jA) < maxTriesW
                if try_node(jD, jA, seeds{q(1),q(2)}, maxIter, 'w1')
                    queue{end+1} = [jD jA]; %#ok<AGROW>
                end
            end
        end
    end
case 'independent'
    for iD = 1:nD
        for iA = 1:nA
            try_node(iD, iA, [], maxIter, 'ind');
        end
    end
case 'chain'
    % legacy anchor-seeded flood (single-hop; the sub-stepped march is retired
    % -- it cannot cross creases and costs too much failing to)
    try_node(1, 1, seed0, maxIter, 'ch');
    queue = {[1 1]};
    while ~isempty(queue)
        q = queue{1};  queue(1) = [];
        for kk = 1:4
            jD = mod(q(1)-1+nbr(kk,1), nD)+1;  jA = mod(q(2)-1+nbr(kk,2), nA)+1;
            if ~PASS(jD,jA) && TRIES(jD,jA) < maxTriesW && ~isempty(seeds{q(1),q(2)})
                if try_node(jD, jA, seeds{q(1),q(2)}, maxIter, 'ch')
                    queue{end+1} = [jD jA]; %#ok<AGROW>
                end
            end
        end
    end
otherwise
    error('sweep_phasing_direct:mode','unknown mode %s', mode);
end

S = struct('sD',sD,'sA',sA,'TF',TF,'MF',MF,'DV',DV,'PERIS',PERIS,'ARRALT',ARRALT, ...
    'ARRSPD',ARRSPD, ...
    'GLOBKM',GLOBKM,'GLOBMS',GLOBMS,'DEFECT',DEFECT,'PASS',PASS,'WALL',WALL, ...
    'LAM0',LAM0, 'LAMT',LAMT, 'anchor',struct('fA0',fA0,'tfRef',sol(8)), ...
    'meta',struct('N',N,'thrustN',thrustN,'ispS',ispS,'m0kg',m0, ...
                  'orbit',struct('tauDRO',tauDRO,'NpTulip',NpT, ...
                       'tauTulip',tauT,'pmTulip',pmT,'muStar',muStar, ...
                       'lStar',lStar,'tStar',tStar,'rMoonKm',rMoonKm, ...
                       'periodDRO',tD(end),'periodTulip',tT(end)), ...
                  'floorKm',floorKm,'globTolKm',globTolKm,'globTolMs',globTolMs, ...
                  'wallMin',toc(tW)/60,'opts',opts));
save(fullfile(rd, sprintf('dsweep_%dx%d.mat', nD, nA)), 'S', '-v7.3');
cellsMeta = S.meta;  cellsGrids = struct('sD',sD,'sA',sA); %#ok<NASGU>
save(fullfile(rd, sprintf('dsweep_%dx%d_cells.mat', nD, nA)), ...
     'CELLS','cellsMeta','cellsGrids','-v7.3');
nOK = nnz(PASS);
log_('DONE: %d/%d solved+gated, %d attempted, %.1f min total\n', ...
    nOK, nD*nA, nnz(~isnan(TF)), toc(tW)/60);
if nOK > 0
    fly = TF;  fly(~PASS) = NaN;
    [tfBest, ib] = min(fly(:));  [bD,bA] = ind2sub([nD nA], ib);
    log_('fastest GATED point: tf=%.6f at (sD,sA)=(%.3f,%.3f), peris %.0f km, arrAlt %.0f km\n', ...
        tfBest, sD(bD), sA(bA), PERIS(bD,bA), ARRALT(bD,bA));
end
end

% ---------------------------------------------------------------------------
function s = local_pf(p)
% LOCAL_PF  PASS/fail marker. INPUTS: p (logical)  OUTPUTS: s [char]
if p, s = 'PASS'; else, s = 'fail'; end
end

% ---------------------------------------------------------------------------
function local_filelog(f, echo, varargin)
% LOCAL_FILELOG  Append-and-flush progress logging (MATLAB -batch buffers
% stdout). INPUTS: f; echo; varargin
fid = fopen(f, 'a');
if fid > 0, fprintf(fid, varargin{:}); fclose(fid); end
if echo, fprintf(varargin{:}); end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
