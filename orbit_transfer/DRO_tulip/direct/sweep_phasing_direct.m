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

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'));  addpath(fullfile(here,'certify'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
rd = fullfile(here,'results');  if ~isfolder(rd), mkdir(rd); end
logF  = fullfile(rd, sprintf('dsweep_%dx%d_progress.txt', nD, nA));
ckptF = fullfile(rd, sprintf('dsweep_%dx%d_ckpt.mat', nD, nA));
log_ = @(varargin) local_filelog(logF, verbose, varargin{:});

muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
rMoonKm = 1737.4;
g0   = 9.80665*tStar^2/(1000*lStar);
Tmax = (thrustN/m0)*tStar^2/(lStar*1000);
c    = (ispS/tStar)*g0;
tW = tic;

%% orbits and anchor (same construction as the marching harness)
[~, rvD0] = pumpkynPie.cr3bp.getDRO(1.0);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, 1.0, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(1.0, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(5*2*pi/6, 7, -1);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, 5*2*pi/6, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(5*2*pi/6, rvT0, muStar);
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

%% BFS over the grid, trajectory-space warm starts
TF = nan(nD,nA); MF = TF; DV = TF; PERIS = TF; ARRALT = TF;
GLOBKM = TF; GLOBMS = TF; DEFECT = TF; WALL = TF; ARRSPD = TF; PASS = false(nD,nA);
LAM0 = nan(8,nD,nA);
seeds = cell(nD,nA);
state = zeros(nD,nA);
queue = {struct('iD',1,'iA',1,'seed',seed0,'pfD',sD(1),'pfA',sA(1))};
nbr = [1 0; -1 0; 0 1; 0 -1];
nDone = 0;

while ~isempty(queue)
    q = queue{1};  queue(1) = [];
    iD = q.iD;  iA = q.iA;
    if state(iD,iA) ~= 0, continue, end
    state(iD,iA) = 1;                        % one attempt per point (direct is robust;
    rv0 = dep(sD(iD));  rvf = arr(sA(iA));   %  failures recorded, not retried endlessly)
    t0 = tic;
    try
        % SUB-STEPPED PHASE MARCH from the parent phasing to the target. One
        % 1/18 arrival step moves the tulip target ~160,000 km -- measured to
        % defeat any warm start outright (the solver wanders to 60-80 ND
        % monsters). So the edge is walked adaptively: try the full step; on a
        % gate-fail, halve the phase step and chain gate-checked intermediate
        % solves as seed carriers. Same medicine as the costate march, but in
        % trajectory space, where steps demonstrably chain.
        if indep
            o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tmax, c, muStar, N, ...
                    [], [], [], struct('maxIter',maxIter,'scheme','hermite-simpson', ...
                    'sundman',true,'returnModel',true,'minAltKm',floorKm));
            C = certify_dro_mintime(o, struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...
                    Tmax, c, struct('tfRef',[],'verbose',false,'posTolKm',inf));
        else
            [o, C] = local_phase_march(q.pfD, q.pfA, q.seed, sD(iD), sA(iA), ...
                dep, arr, Tmax, c, muStar, N, maxIter, floorKm, ...
                globTolKm, globTolMs, lStar, tStar);
        end
        w = toc(t0);
        r2 = vecnorm(o.X(1:3,:) - [1-muStar;0;0], 2, 1);
        TF(iD,iA) = o.tf;  MF(iD,iA) = o.mf;
        DV(iD,iA) = c*lStar/tStar*log(1/o.mf);
        PERIS(iD,iA)  = min(r2)*lStar - rMoonKm;
        ARRALT(iD,iA) = norm(rvf(1:3) - [1-muStar 0 0])*lStar - rMoonKm;
        ARRSPD(iD,iA) = norm(rvf(4:6))*lStar/tStar;      % km/s: the difficulty hypothesis variable
        GLOBKM(iD,iA) = C.globKm;  GLOBMS(iD,iA) = C.gates(strcmp({C.gates.id},'G1bv')).value;
        DEFECT(iD,iA) = o.maxDefect;  WALL(iD,iA) = w;
        PASS(iD,iA) = o.success && C.globKm < globTolKm && GLOBMS(iD,iA) < globTolMs ...
                      && o.maxDefect < 1e-9;
        % THE CATALOG ENTRY: sign-resolved initial costates from the duals.
        % Sign convention resolved against the primal control (dual-free).
        if ~isempty(o.lamDef)
            lamD = o.lamDef;
            aS = o.Um(1:3,1)/max(norm(o.Um(1:3,1)),eps);
            sgn = 1;  if (-lamD(4:6,1)/max(norm(lamD(4:6,1)),eps)).'*aS < 0, sgn = -1; end
            LAM0(:,iD,iA) = [sgn*lamD(:,1); o.tf];
        end
        % uniform-time resample of this solution = the neighbours' warm start
        tu = linspace(0, o.tf, N+1);
        sd.X = interp1(o.tNodes, o.X.', tu, 'spline').';
        sd.U = interp1(o.tNodes, o.U.', tu, 'spline').';
        nrm = max(vecnorm(sd.U(1:3,:),2,1), eps);
        sd.U(1:3,:) = sd.U(1:3,:)./nrm;  sd.U(4,:) = min(max(sd.U(4,:),0),1);
        sd.tf = o.tf;
        seeds{iD,iA} = sd;
        % SEED HYGIENE (pilot finding: junk begets junk). A gated-fail node's
        % own trajectory must not warm-start its neighbours; hand the parent's
        % seed onward instead.
        if PASS(iD,iA), sdPush = sd; else, sdPush = q.seed; end
        nDone = nDone + 1;
        log_('  (%d,%d) tf=%.5f dV=%.4f peris=%6.0f arrAlt=%6.0f vArr=%.3f glob=%.3fkm %s [%d/%d, %.0fs]\n', ...
            iD, iA, o.tf, DV(iD,iA), PERIS(iD,iA), ARRALT(iD,iA), ARRSPD(iD,iA), C.globKm, ...
            local_pf(PASS(iD,iA)), nDone, nD*nA, w);
        for kk = 1:4
            jD = mod(iD-1+nbr(kk,1), nD)+1;  jA = mod(iA-1+nbr(kk,2), nA)+1;
            if state(jD,jA) == 0
                if PASS(iD,iA), pfD = sD(iD); pfA = sA(iA); else, pfD = q.pfD; pfA = q.pfA; end
                queue{end+1} = struct('iD',jD,'iA',jA,'seed',sdPush,'pfD',pfD,'pfA',pfA); %#ok<AGROW>
            end
        end
        if mod(nDone,3) == 0
            save(ckptF,'TF','MF','DV','PERIS','ARRALT','ARRSPD','GLOBKM','GLOBMS', ...
                 'DEFECT','PASS','WALL','LAM0','sD','sA');
        end
    catch ME
        log_('  (%d,%d) FAILED after %.0fs: %s\n', iD, iA, toc(t0), ME.message);
        for kk = 1:4                          % still expand so the map routes around
            jD = mod(iD-1+nbr(kk,1), nD)+1;  jA = mod(iA-1+nbr(kk,2), nA)+1;
            if state(jD,jA) == 0
                queue{end+1} = struct('iD',jD,'iA',jA,'seed',q.seed, ...
                    'pfD',q.pfD,'pfA',q.pfA); %#ok<AGROW>
            end
        end
    end
end

S = struct('sD',sD,'sA',sA,'TF',TF,'MF',MF,'DV',DV,'PERIS',PERIS,'ARRALT',ARRALT, ...
    'ARRSPD',ARRSPD, ...
    'GLOBKM',GLOBKM,'GLOBMS',GLOBMS,'DEFECT',DEFECT,'PASS',PASS,'WALL',WALL, ...
    'LAM0',LAM0, 'anchor',struct('fA0',fA0,'tfRef',sol(8)), ...
    'meta',struct('N',N,'thrustN',thrustN,'ispS',ispS,'m0kg',m0, ...
                  'wallMin',toc(tW)/60,'opts',opts));
save(fullfile(rd, sprintf('dsweep_%dx%d.mat', nD, nA)), 'S', '-v7.3');
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
function [o, C] = local_phase_march(pfD, pfA, seed, fDt, fAt, dep, arr, ...
    Tmax, c, muStar, N, maxIter, floorKm, globTolKm, globTolMs, lStar, tStar)
% LOCAL_PHASE_MARCH  Adaptive phase-space continuation between two phasings,
% every sub-step a full gate-checked direct solve, seeds chained.
% INPUTS:  parent phases + seed; target phases; handles and solver params
% OUTPUTS: o, C for the TARGET phasing (the final accepted solve; if the march
%          dies early, the last attempt at wherever it died -- caller gates it)
dD = mod(fDt - pfD + 0.5, 1) - 0.5;
dA = mod(fAt - pfA + 0.5, 1) - 0.5;
L  = max(abs(dD), abs(dA));  if L == 0, L = eps; end
hMin = 1/160;  budget = 14;  tEdge = tic;
pos = 0;  h = L;  nUsed = 0;
while pos < L
    pos2 = min(pos + h, L);
    fD = pfD + dD*(pos2/L);  fA = pfA + dA*(pos2/L);
    rv0 = dep(fD);  rvf = arr(fA);
    o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tmax, c, muStar, N, ...
            seed.X, seed.U, seed.tf, ...
            struct('maxIter',maxIter,'scheme','hermite-simpson', ...
                   'sundman',true,'returnModel',true,'minAltKm',floorKm));
    C = certify_dro_mintime(o, struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...
            Tmax, c, struct('tfRef',[],'verbose',false,'posTolKm',inf));
    nUsed = nUsed + 1;
    okS = o.success && C.globKm < globTolKm && ...
          C.gates(strcmp({C.gates.id},'G1bv')).value < globTolMs && o.maxDefect < 1e-9;
    if okS
        pos = pos2;
        tu = linspace(0, o.tf, N+1);
        seed.X = interp1(o.tNodes, o.X.', tu, 'spline').';
        seed.U = interp1(o.tNodes, o.U.', tu, 'spline').';
        nrm = max(vecnorm(seed.U(1:3,:),2,1), eps);
        seed.U(1:3,:) = seed.U(1:3,:)./nrm;  seed.U(4,:) = min(max(seed.U(4,:),0),1);
        seed.tf = o.tf;
        h = min(h*2, L - pos + eps);
    else
        h = h/2;
        if h < hMin || nUsed >= budget || toc(tEdge) > 1800, return, end
    end
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
