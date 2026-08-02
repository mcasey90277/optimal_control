function S = sweep_phasing(opts)
% SWEEP_PHASING  Map the DRO->tulip min-time family over departure/arrival phase.
%
% Darin's program (call of 2026-08-03): map the indirect method's costates over
% essentially all initial and final orbit phasings, so the extremal FAMILY is
% seen rather than sampled -- single point solves land in arbitrary basins, as
% the direct campaign's cold-start experiment measured. Each converged point
% seeds its neighbours (continuation), which is the costate-catalog walk
% mechanism generalized from abstracts/data/bht1500_continuation.m to two
% dimensions.
%
% WHAT IT DOES
%   1. Builds the departure DRO and arrival tulip once, each propagated over one
%      period. Phase = fraction of period, so the grid is [0,1) x [0,1) with
%      wraparound.
%   2. Breadth-first CONTINUATION from the demo's phasing pair, whose converged
%      costates are known: each grid point is solved by pumpkyn.cr3bp.tfMin
%      seeded from the converged costates of the neighbour it was reached from.
%      Failed points are retried from other converged neighbours.
%   3. VERIFIES every point independently -- tfMin returns no convergence flag,
%      so each solution is re-propagated and the terminal state error measured.
%      Unverified points are recorded as failures, never silently kept.
%   4. Records, per point: t_f, the eight costates (THE CATALOG), terminal
%      residual, PERISELENE, revolutions about the Moon, final mass, dV.
%   5. Optional FOLD PASS: where t_f jumps between neighbours, re-solve from
%      every converged neighbour separately. Distinct verified answers mean two
%      sheets of the solution surface cross there; the second sheet is stored
%      rather than discarded.
%
% TWO CHANNELS, NOT ONE, and the reason is measured: the indirect solver is
% UNCONSTRAINED, and the min-time family improves with deeper lunar flybys (the
% direct campaign's floor experiment found a floor-riding extremal UNDER the
% reference t_f). Some phasings will converge to extremals with periselene at
% hundreds of km -- or below the surface, which the point-mass model permits. A
% dV map alone would present those as "the global minimum". The periselene map
% and the flyability mask are therefore first-class outputs.
%
% INPUTS:
%   opts - struct (optional):
%          .nD, .nA     grid sizes, departure x arrival     [default 12, 12]
%          .tauDRO      DRO period parameter                [default 1.0]
%          .NpTulip     tulip petal count                   [default 7]
%          .tauTulip    tulip period parameter              [default 5*2*pi/6]
%          .pmTulip     tulip hemisphere                    [default -1]
%          .thrustN     thrust, N                           [default 0.07]
%          .ispS        specific impulse, s                 [default 900]
%          .m0kg        wet mass, kg                        [default 150]
%          .termTol     verification: terminal 6-state error [default 1e-7]
%          .maxTries    seed attempts per point before giving up [default 4]
%          .foldTol     t_f jump flagged as a fold, ND      [default 0.05]
%          .foldPass    run the multi-neighbour fold pass   [default true]
%          .perisFloorKm flyability mask level, km          [default 500]
%          .saveMat     write results/<name>.mat            [default true]
%          .plots       render the four maps                [default true]
%          .verbose                                         [default true]
%
% OUTPUTS:
%   S - struct:
%       .sD, .sA     phase grids [1xnD], [1xnA] (fractions of period)
%       .TF          final times [nD x nA], NaN where unverified
%       .LAM         costates    [8 x nD x nA] -- the catalog
%       .PERIS       periselene altitude, km [nD x nA]
%       .REVS        revolutions about the Moon [nD x nA]
%       .DV          delta-V, km/s [nD x nA]
%       .RESID       terminal verification residual [nD x nA]
%       .TRIES       seed attempts used [nD x nA]
%       .FOLD        fold flag [nD x nA logical]
%       .TF2, .LAM2  second sheet where the fold pass found one (NaN elsewhere)
%       .anchor      the demo phasing pair and costates the walk started from
%       .meta        constants, options, timing
%
% REFERENCES:
%   [1] abstracts/data/bht1500_continuation.m -- the 1-D walk this generalizes.
%   [2] DRO_tulip/FINDINGS.md -- basin multiplicity + floor results motivating
%       the fold pass and the periselene channel.
%   [3] pumpkyn.cr3bp.tfMin / tfMinProp (the solver and its propagator).

if nargin < 1, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
nD       = d('nD', 12);        nA       = d('nA', 12);
tauDRO   = d('tauDRO', 1.0);   NpT      = d('NpTulip', 7);
tauT     = d('tauTulip', 5*2*pi/6);  pmT = d('pmTulip', -1);
thrustN  = d('thrustN', 0.07); ispS     = d('ispS', 900);
m0       = d('m0kg', 150);
termTol  = d('termTol', 1e-7); maxTries = d('maxTries', 4);
foldTol  = d('foldTol', 0.05); foldPass = d('foldPass', true);
perisFloorKm = d('perisFloorKm', 500);
saveMat  = d('saveMat', true); doPlots  = d('plots', true);
verbose  = d('verbose', true);

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;
rMoonKm= 1737.4;
g0   = 9.80665*tStar^2/(1000*lStar);
Tmax = (thrustN/m0)*tStar^2/(lStar*1000);
c    = (ispS/tStar)*g0;
tWall = tic;

%% the two orbits, each over one period
[~, rvD0] = pumpkynPie.cr3bp.getDRO(tauDRO);
rvD0      = pumpkyn.cr3bp.cont_np(rvD0, tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, NpT, pmT);
rvT0      = pumpkyn.cr3bp.cont_np(rvT0, tauT, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, muStar);
depState = @(f) interp1(tD, rvD, mod(f,1)*tD(end), 'spline');
arrState = @(f) interp1(tT, rvT, mod(f,1)*tT(end), 'spline');

%% the anchor: the demo's phasing pair, whose costates are known
% Departure = the DRO propagation start (fraction 0). Arrival = the demo's
% max-velocity-angle pick, converted to a fraction of the tulip period.
dvTheta   = pumpkyn.util.bsxAng(rvT(:,4:6), rvD(1,4:6), 2);
[~, idxA] = max(dvTheta);
fA0       = tT(idxA)/tT(end);
lamAnchor = [13.9579470518969; 6.26766054139842; 7.02087277890196; ...
             -0.969930036390389; -5.04449989501292; -3.11628049459521; ...
             5.50150618633038; 4.01524259262941];

% THE GRID IS ALIGNED TO THE ANCHOR. The first version used a plain (0:n-1)/n
% grid and snapped the anchor to its nearest node -- which sat on the far side
% of a continuation barrier (measured near arrival phase ~0.028), so the very
% first solve failed and the whole walk starved. Offsetting the grid so the
% anchor IS a node makes the first solve trivially convergent, and the BFS then
% routes around barriers on the torus instead of dying on one.
sD = mod(0   + (0:nD-1)/nD, 1);
sA = mod(fA0 + (0:nA-1)/nA, 1);
iD0 = 1;  iA0 = 1;

if verbose
    fprintf('\n===== PHASING SWEEP: DRO(tau0=%.3g) -> tulip(Np=%d), %dx%d grid =====\n', ...
        tauDRO, NpT, nD, nA);
    fprintf('  anchor: (sD, sA) = (%.3f, %.3f), demo costates; thrust %.3g N\n', ...
        sD(iD0), sA(iA0), thrustN);
end

%% breadth-first continuation over the grid
TF    = nan(nD,nA);   LAM   = nan(8,nD,nA);  PERIS = nan(nD,nA);
REVS  = nan(nD,nA);   DV    = nan(nD,nA);    RESID = nan(nD,nA);
TRIES = zeros(nD,nA); state = zeros(nD,nA);  % 0 open, 1 converged, -1 dead
% Each edge of the walk is a CONTINUATION in phase space, not a single solve.
% Measured necessity: the tulip moves 1.23 ND (~477,000 km) across one 6-grid
% arrival step, and tfMin cannot jump that from the neighbour's costates in one
% call (the first version of this harness tried, and the anchor itself failed).
% local_continue bisects the phase step until each sub-step converges, exactly
% the adaptive walk bht1500_continuation.m uses for thrust.
queue = {struct('iD',iD0,'iA',iA0,'pfD',0,'pfA',fA0,'lam',lamAnchor)};
nbr   = [1 0; -1 0; 0 1; 0 -1];

while ~isempty(queue)
    q = queue{1};  queue(1) = [];
    iD = q.iD;  iA = q.iA;
    if state(iD,iA) ~= 0, continue, end
    TRIES(iD,iA) = TRIES(iD,iA) + 1;

    [ok, sol, res, peris, revs, mf] = local_continue( ...
        q.pfD, q.pfA, q.lam, sD(iD), sA(iA), depState, arrState, ...
        Tmax, c, muStar, termTol, lStar, rMoonKm, 6);

    if ok
        state(iD,iA) = 1;
        TF(iD,iA)    = sol(8);   LAM(:,iD,iA) = sol;
        PERIS(iD,iA) = peris;    REVS(iD,iA)  = revs;
        DV(iD,iA)    = c*lStar/tStar*log(1/mf);
        RESID(iD,iA) = res;
        for kk = 1:4                          % expand, with wraparound
            jD = mod(iD-1+nbr(kk,1), nD)+1;  jA = mod(iA-1+nbr(kk,2), nA)+1;
            if state(jD,jA) == 0
                queue{end+1} = struct('iD',jD,'iA',jA, ...
                    'pfD',sD(iD),'pfA',sA(iA),'lam',sol); %#ok<AGROW>
            end
        end
        if verbose
            fprintf('  (%2d,%2d) tf=%.4f peris=%7.0f km  [%d/%d]\n', ...
                iD, iA, sol(8), peris, nnz(state==1), nD*nA);
        end
    elseif TRIES(iD,iA) >= maxTries
        state(iD,iA) = -1;                    % dead; stop re-queueing it
    end
    % a failed-but-alive point is NOT re-queued here -- it re-enters when
    % another neighbour converges and pushes it again with a fresh seed
end

%% fold detection, and optionally the second sheet
FOLD = false(nD,nA);
for iD = 1:nD
    for iA = 1:nA
        if isnan(TF(iD,iA)), continue, end
        for kk = 1:4
            jD = mod(iD-1+nbr(kk,1), nD)+1;  jA = mod(iA-1+nbr(kk,2), nA)+1;
            if ~isnan(TF(jD,jA)) && abs(TF(iD,iA)-TF(jD,jA)) > foldTol
                FOLD(iD,iA) = true;
            end
        end
    end
end

TF2 = nan(nD,nA);  LAM2 = nan(8,nD,nA);
if foldPass && any(FOLD(:))
    if verbose, fprintf('  fold pass: %d flagged points\n', nnz(FOLD)); end
    for iD = 1:nD
        for iA = 1:nA
            if ~FOLD(iD,iA), continue, end
            for kk = 1:4
                jD = mod(iD-1+nbr(kk,1), nD)+1;  jA = mod(iA-1+nbr(kk,2), nA)+1;
                if isnan(TF(jD,jA)), continue, end
                [ok, sol] = local_continue(sD(jD), sA(jA), LAM(:,jD,jA), ...
                    sD(iD), sA(iA), depState, arrState, ...
                    Tmax, c, muStar, termTol, lStar, rMoonKm, 4);
                if ~ok, continue, end
                if sol(8) < TF(iD,iA) - 1e-6
                    % a faster verified extremal: promote it, demote the old
                    TF2(iD,iA) = TF(iD,iA);   LAM2(:,iD,iA) = LAM(:,iD,iA);
                    TF(iD,iA)  = sol(8);      LAM(:,iD,iA)  = sol;
                elseif abs(sol(8) - TF(iD,iA)) > foldTol && isnan(TF2(iD,iA))
                    TF2(iD,iA) = sol(8);      LAM2(:,iD,iA) = sol;  % second sheet
                end
            end
        end
    end
end

%% assemble, save, plot
S = struct('sD',sD, 'sA',sA, 'TF',TF, 'LAM',LAM, 'PERIS',PERIS, 'REVS',REVS, ...
    'DV',DV, 'RESID',RESID, 'TRIES',TRIES, 'FOLD',FOLD, 'TF2',TF2, 'LAM2',LAM2, ...
    'anchor', struct('iD',iD0,'iA',iA0,'fA0',fA0,'lam',lamAnchor), ...
    'meta', struct('muStar',muStar,'lStar',lStar,'tStar',tStar,'Tmax',Tmax, ...
        'c',c,'thrustN',thrustN,'ispS',ispS,'m0kg',m0,'tauDRO',tauDRO, ...
        'NpTulip',NpT,'tauTulip',tauT,'pmTulip',pmT,'opts',opts, ...
        'wallSec',toc(tWall)));

here = fileparts(mfilename('fullpath'));
if saveMat
    rd = fullfile(here,'results');  if ~isfolder(rd), mkdir(rd); end
    fn = fullfile(rd, sprintf('sweep_tau%g_%dx%d.mat', tauDRO, nD, nA));
    save(fn, 'S', '-v7.3');
    if verbose, fprintf('  results -> %s\n', fn); end
end

nOK = nnz(~isnan(TF));
if verbose
    fprintf(['  %d/%d verified (%d dead, %d folds, %d second-sheet), ' ...
        '%.1f min wall\n'], nOK, nD*nA, nnz(state==-1), nnz(FOLD), ...
        nnz(~isnan(TF2)), toc(tWall)/60);
    if nOK > 0
        [tfBest, ib] = min(TF(:));  [bD,bA] = ind2sub([nD nA], ib);
        fprintf('  fastest verified: tf = %.6f at (sD,sA) = (%.3f, %.3f), peris %.0f km', ...
            tfBest, sD(bD), sA(bA), PERIS(bD,bA));
        if PERIS(bD,bA) < perisFloorKm
            fprintf('  <-- BELOW the %d km flyability floor\n', perisFloorKm);
        else
            fprintf('\n');
        end
        flyable = TF;  flyable(PERIS < perisFloorKm) = NaN;
        if any(~isnan(flyable(:)))
            [tfFly, ifb] = min(flyable(:));  [fD,fA] = ind2sub([nD nA], ifb);
            fprintf('  fastest FLYABLE (peris >= %d km): tf = %.6f at (%.3f, %.3f), peris %.0f km\n', ...
                perisFloorKm, tfFly, sD(fD), sA(fA), PERIS(fD,fA));
        end
    end
end

if doPlots && nOK > 0
    fig = figure('Color','w','Position',[50 50 1500 850]);
    local_map(1, sD, sA, TF.',   'final time t_f [ND]');
    hold on;  spy2(FOLD.', sD, sA);
    local_map(2, sD, sA, DV.',   '\DeltaV [km/s]');
    local_map(3, sD, sA, PERIS.','periselene altitude [km]');
    hold on;  contour(sD, sA, PERIS.', [perisFloorKm perisFloorKm], 'r-', 'LineWidth',2);
    fly = TF;  fly(PERIS < perisFloorKm) = NaN;
    local_map(4, sD, sA, fly.',  sprintf('t_f, FLYABLE only (peris \\geq %d km)', perisFloorKm));
    sgtitle(sprintf(['DRO(\\tau_0=%.3g) \\rightarrow tulip: min-time family over ' ...
        'phasing  |  %d/%d verified, %d folds'], tauDRO, nOK, nD*nA, nnz(FOLD)), ...
        'FontWeight','bold');
    if saveMat
        exportgraphics(fig, fullfile(here,'results', ...
            sprintf('sweep_tau%g_%dx%d.png', tauDRO, nD, nA)), 'Resolution',140);
    end
end
end

% ---------------------------------------------------------------------------
function [ok, sol, res, perisKm, revs, mf] = local_continue( ...
    fDa, fAa, lamA, fDb, fAb, depState, arrState, ...
    Tmax, c, muStar, termTol, lStar, rMoonKm, maxDepth)
% LOCAL_CONTINUE  Walk the costates from phase pair (fDa,fAa) to (fDb,fAb).
%
% Try the full step first. If the solve fails verification, split the phase
% step in half (shortest way around the circle in each coordinate) and recurse:
% converge to the midpoint, then from the midpoint to the target. maxDepth = 5
% permits up to 32 sub-steps per grid edge, each an independent solve+verify.
% This is the 2-D generalization of bht1500_continuation.m's thrust walk.
%
% INPUTS:  fDa,fAa - converged phases; lamA [8x1] their costates
%          fDb,fAb - target phases; depState/arrState - phase -> state handles
%          Tmax; c; muStar; termTol; lStar; rMoonKm; maxDepth [scalar]
% OUTPUTS: as local_solve_verify, for the TARGET point
res = inf;  perisKm = NaN;  revs = NaN;  mf = NaN;
rvf = arrState(fAb);  rv0 = depState(fDb);
[ok, sol, res, perisKm, revs, mf] = local_solve_verify( ...
    rv0(1:6), rvf(1:6), lamA, Tmax, c, muStar, termTol, lStar, rMoonKm);
if ok || maxDepth <= 0, return, end
% bisect along the shortest circular path in each phase
fDm = fDa + (mod(fDb - fDa + 0.5, 1) - 0.5)/2;
fAm = fAa + (mod(fAb - fAa + 0.5, 1) - 0.5)/2;
[okm, solm] = local_continue(fDa, fAa, lamA, fDm, fAm, depState, arrState, ...
    Tmax, c, muStar, termTol, lStar, rMoonKm, maxDepth-1);
if ~okm, ok = false; sol = nan(8,1); return, end
[ok, sol, res, perisKm, revs, mf] = local_continue(fDm, fAm, solm, fDb, fAb, ...
    depState, arrState, Tmax, c, muStar, termTol, lStar, rMoonKm, maxDepth-1);
end

% ---------------------------------------------------------------------------
function [ok, sol, res, perisKm, revs, mf] = local_solve_verify( ...
    rv0, rvf, lamSeed, Tmax, c, muStar, termTol, lStar, rMoonKm)
% LOCAL_SOLVE_VERIFY  One tfMin solve plus the independent check it lacks.
% tfMin returns no convergence flag, so every solution is re-propagated and the
% terminal state compared against the target. Anything else is not a result.
% INPUTS:  rv0, rvf [1x6]; lamSeed [8x1]; Tmax; c; muStar; termTol; lStar; rMoonKm
% OUTPUTS: ok [logical]; sol [8x1]; res (terminal 6-state error norm);
%          perisKm; revs; mf (final mass fraction)
sol = nan(8,1);  res = inf;  perisKm = NaN;  revs = NaN;  mf = NaN;
try
    w = warning('off','all');  cleanW = onCleanup(@() warning(w)); %#ok<NASGU>
    sol = pumpkyn.cr3bp.tfMin(rv0, rvf, lamSeed, Tmax, c, muStar);
    if ~all(isfinite(sol)) || sol(8) <= 0, ok = false; return, end
    [~, rv] = pumpkyn.cr3bp.tfMinProp(sol(8), [rv0, 1, sol(1:7)'], Tmax, c, muStar);
    res = norm(rv(end,1:6) - rvf);
    r2  = vecnorm(rv(:,1:3) - [1-muStar 0 0], 2, 2);
    perisKm = min(r2)*lStar - rMoonKm;
    th   = unwrap(atan2(rv(:,2), rv(:,1)-(1-muStar)));
    revs = (max(th)-min(th))/(2*pi);
    mf   = rv(end,7);
    ok   = res < termTol;
catch
    ok = false;
end
end

% ---------------------------------------------------------------------------
function local_map(k, sD, sA, Z, ttl)
% LOCAL_MAP  One pcolor panel of the sweep. INPUTS: k (subplot); sD; sA; Z; ttl
subplot(2,2,k);
pcolor(sD, sA, Z);  shading flat;  colorbar;
xlabel('departure phase s_D');  ylabel('arrival phase s_A');
title(ttl);
end

% ---------------------------------------------------------------------------
function spy2(M, sD, sA)
% SPY2  Overlay fold markers on the current map. INPUTS: M [nAxnD logical]; sD; sA
[ia, id] = find(M);
if ~isempty(ia), plot(sD(id), sA(ia), 'k.', 'MarkerSize', 10); end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
