function R = run_minenergy_pilot(cfg)
% RUN_MINENERGY_PILOT  Fixed-tf MIN-ENERGY pilot on flagship 12x12 DRO->tulip cells.
%
% The first non-min-time run of the costate pipeline: for a few flagship
% torus cells, at t_f = gamma * t_f^min(cell),
%
%   (1) DIRECT  casadi_mintime_dro(objective='energy', tfFix) warm-started
%               from the stored min-time cell (Sundman HS, N = 800, 500 km
%               floor) -> X, U, Um, lamDef, J = Int s^2 dt;
%   (2) HARVEST harvest_ms_seed (Hager covector mapping, rows 1:7);
%   (3) REFINE  ms_minenergy over the K-ladder 12 -> 24 -> 48 (fixed-tf
%               ms_bvp), z7 = lambda_0;
%   (4) GATES   G1b flown control on the direct solution (< 100 km / 10 m/s);
%               single-shooting acceptance ss_bvp_accept (|dz| < 1e-6 at
%               the single-shooting residual floor 1e-6 -- the role pumpkyn
%               tfMin plays for min-time); Hamiltonian first integral along
%               the indirect flight (ABSOLUTE drift < 1e-6; measured
%               3e-9..4e-8 at ode113 RelTol 1e-10 over 18-30 d);
%               direct-vs-indirect cost agreement; throttle interior
%               somewhere (else it is min-time in disguise).
%
% Every (cell, gamma) record is saved on completion to
% direct/results/minenergy_pilot.mat (resumable: existing records are
% skipped unless cfg.redo). Figures: direct/results/minenergy_<iD>_<iA>_g<gam>.png.
%
% INPUTS:
%   cfg - (optional) struct overriding the ADJUSTABLE PARAMETERS below:
%         .cells [n x 2] (iD, iA) torus indices; .gammas [1 x m] t_f
%         multipliers applied to EVERY cell; .Kladder; .maxIter; .redo;
%         .regate (reuse a stored direct solution, redo harvest/ms/gates);
%         .stage 'run' (default) | 'report' (print the table only)
%
% OUTPUTS:
%   R - struct array of records (one per (cell, gamma)); also saved
%
% REFERENCES:
%   [1] process/COSTATE_LIBRARY_PIPELINE.md (the three-step pipeline).
%   [2] costate_common/cr3bp_minenergy_pmp.m (the PMP field + conventions).

%% ADJUSTABLE PARAMETERS ---------------------------------------------------
P.cells   = [2 5; 6 8; 1 2];     % golden harvest cell; a short-tf cell; a long-tf cell
P.gammas  = 1.2;                 % t_f / t_f^min; ladder e.g. [1.1 1.2 1.4]
P.Kladder = [12 24 48];          % ms segmentation escalation
P.maxIter = 1500;                % IPOPT budget for the direct energy solve
P.maxCpuSec = 1200;              % IPOPT wall cap per solve
P.tolR    = 1e-10;               % ms residual tolerance
P.tolDz   = 1e-6;                % acceptance move tolerance
P.globTolKm = 100;  P.globTolMs = 10;   % flown-arrival tolerances (map tolerance)
P.redo    = false;
P.regate  = false;
P.stage   = 'run';
if nargin >= 1 && ~isempty(cfg)
    fn = fieldnames(cfg);
    for k = 1:numel(fn), P.(fn{k}) = cfg.(fn{k}); end
end

%% Paths -------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'direct', 'lib'));
addpath(fullfile(here, 'indirect'));
addpath(fullfile(fileparts(here), 'costate_common'));
addpath(fullfile(fileparts(fileparts(here)), 'oclib'));
if isempty(which('pumpkyn.cr3bp.tfMinEoM'))
    run(fullfile(fileparts(here), 'cr3bp_common', 'setup_cr3bp_common.m'));
end
if isempty(which('casadi.Opti')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end
resDir  = fullfile(here, 'direct', 'results');
outMat  = fullfile(resDir, 'minenergy_pilot.mat');
S = load(fullfile(resDir, 'dsweep_12x12_cells.mat'));
lStar = S.cellsMeta.orbit.lStar;  tStar = S.cellsMeta.orbit.tStar;
kmPerNd = lStar;  msPerNd = 1e3*lStar/tStar;

R = struct([]);
if exist(outMat, 'file'), L = load(outMat); R = L.R; end
if strcmp(P.stage, 'report'), report(R); return, end

%% Loop --------------------------------------------------------------------
for ic = 1:size(P.cells, 1)
    iD = P.cells(ic,1);  iA = P.cells(ic,2);
    cc = S.CELLS{iD, iA};
    if isempty(cc), fprintf('cell (%d,%d) empty in the torus -- skip\n', iD, iA); continue, end
    for gam = P.gammas
        tag = sprintf('(%d,%d) gamma=%.2f', iD, iA, gam);
        hit = [];
        if ~isempty(R), hit = find([R.iD]==iD & [R.iA]==iA & abs([R.gam]-gam)<1e-12, 1); end
        if ~P.redo && ~P.regate && ~isempty(hit)
            fprintf('%s already done -- skip (cfg.redo / cfg.regate)\n', tag);  continue
        end
        fprintf('\n=== %s : tf_min = %.4f ND -> tf = %.4f ND (%.2f d) ===\n', ...
                tag, cc.tf, gam*cc.tf, gam*cc.tf*tStar/86400);
        rec = struct('iD',iD,'iA',iA,'gam',gam,'tfMin',cc.tf,'tf',gam*cc.tf, ...
                     'Tmax',cc.Tmax,'c',cc.c,'muStar',cc.muStar,'rv0',cc.rv0,'rvf',cc.rvf);

        % --- (1) direct energy solve, warm from the min-time cell ----------
        nN = cc.N + 1;
        tU = linspace(0, cc.tf, nN);
        X0 = interp1(cc.tNodes, cc.X', tU, 'pchip')';
        U0 = interp1(cc.tNodes, cc.U', tU, 'pchip')';
        U0(1:3,:) = U0(1:3,:) ./ max(vecnorm(U0(1:3,:), 2, 1), eps);
        U0(4,:)   = min(max(U0(4,:), 0), 1);
        if P.regate && ~isempty(hit) && isfield(R(hit), 'X') && ~isempty(R(hit).X)
            % reuse the stored direct solution (the expensive step)
            r0 = R(hit);
            o = struct('X',r0.X,'U',r0.U,'Um',r0.Um,'tNodes',r0.tNodes,'lamDef',r0.lamDef, ...
                'tf',r0.tf,'s',linspace(0,1,size(r0.X,2)),'success',r0.direct.success, ...
                'ipoptStatus',r0.direct.status,'J',r0.direct.J,'mf',r0.direct.mf, ...
                'thrMin',r0.direct.thrMin,'maxDefect',r0.direct.maxDefect, ...
                'altMinKm',r0.direct.altMinKm);
            rec.direct = r0.direct;
            fprintf('  direct (stored): %s  J = %.6f  mf = %.6f\n', o.ipoptStatus, o.J, o.mf);
        else
            tic
            o = casadi_mintime_dro(cc.rv0(:), cc.rvf(:), cc.Tmax, cc.c, cc.muStar, cc.N, ...
                    X0, U0, gam*cc.tf, struct('scheme','hermite-simpson','sundman',true, ...
                    'minAltKm',cc.floorKm,'returnModel',true,'objective','energy', ...
                    'tfFix',gam*cc.tf,'maxIter',P.maxIter,'maxCpuSec',P.maxCpuSec));
            rec.direct = struct('success',o.success,'status',o.ipoptStatus,'wall',toc, ...
                'J',o.J,'mf',o.mf,'thrMin',o.thrMin,'thrMax',max(o.U(4,:)), ...
                'maxDefect',o.maxDefect,'altMinKm',o.altMinKm,'tf',o.tf);
            fprintf('  direct: %s  J = %.6f  mf = %.6f  thr in [%.3f, %.3f]  alt %.0f km  %.0f s\n', ...
                    o.ipoptStatus, o.J, o.mf, o.thrMin, rec.direct.thrMax, o.altMinKm, rec.direct.wall);
        end
        rec.X = o.X;  rec.U = o.U;  rec.Um = o.Um;  rec.tNodes = o.tNodes;  rec.lamDef = o.lamDef;
        if ~o.success
            rec.verdict = 'DIRECT FAIL';  R = append(R, rec);  save(outMat, 'R');  continue
        end
        % --- G1b: fly the direct control ---------------------------------
        [erNd, evNd] = flown_control_error(o, cc.muStar, cc.Tmax, cc.c);
        rec.direct.flownKm = erNd*kmPerNd;  rec.direct.flownMs = evNd*msPerNd;
        fprintf('  G1b flown: %.2f km / %.3f m/s\n', rec.direct.flownKm, rec.direct.flownMs);

        % --- (2) harvest --------------------------------------------------
        o7 = o;  o7.lamDef = o.lamDef(1:7,:);       % no lambda_t=+1 check: fixed tf
        % --- (3) refine over the K-ladder ---------------------------------
        msBest = [];  Kused = NaN;
        for K = P.Kladder
            [seed, hd] = harvest_ms_seed(o7, K);
            [z, mi] = ms_minenergy(cc.rv0(:), cc.rvf(:), gam*cc.tf, seed, cc.Tmax, cc.c, ...
                        cc.muStar, struct('tolR',P.tolR,'accept',true,'tolDz',P.tolDz));
            fprintf('  ms K=%2d: conv=%d ||R||=%.1e iters=%d %.1fs  Hdrift=%.1e  accept: dz=%.1e normR0=%.1e normR=%.1e %d\n', ...
                    K, mi.converged, mi.normR, mi.iters, mi.wall, mi.Hdrift, ...
                    mi.accept.dz, mi.accept.normR0, mi.accept.normR, mi.accept.accepted);
            if isempty(msBest) || mi.normR < msBest.normR
                msBest = mi;  msBest.z = z;  msBest.K = K;  msBest.hd = hd;  Kused = K;
            end
            if mi.converged, break, end
        end
        rec.harvest = struct('sign',msBest.hd.sign,'voteMargin',msBest.hd.voteMargin);
        rec.ms = struct('z',msBest.z,'K',Kused,'converged',msBest.converged,'normR',msBest.normR, ...
            'iters',msBest.iters,'wall',msBest.wall,'Hjunction',msBest.H,'Hdrift',msBest.Hdrift, ...
            'sJunction',msBest.s,'acceptDz',msBest.accept.dz,'acceptNormR0',msBest.accept.normR0, ...
            'acceptNormR',msBest.accept.normR,'accepted',msBest.accept.accepted, ...
            'acceptIters',msBest.accept.iters);
        % lambda_0: harvested seed vs converged
        [seed1, ~] = harvest_ms_seed(o7, Kused);
        lamSeed = seed1.Y(8:14,1);
        rec.ms.lam0SeedRelErr = max(abs(lamSeed - msBest.z)) / max(abs(msBest.z));

        % --- indirect flight from the converged lambda_0 ------------------
        y1 = [cc.rv0(:); 1; msBest.z];
        [yh, ~, T, Y] = cr3bp_minenergy_prop(gam*cc.tf, y1, false, cc.Tmax, cc.c, cc.muStar);
        nT = numel(T);  sT = zeros(nT,1);  HT = zeros(nT,1);
        for k = 1:nT
            [~, ~, ax] = cr3bp_minenergy_pmp(Y(k,:)', cc.Tmax, cc.c, cc.muStar);
            sT(k) = ax.s;  HT(k) = ax.H;
        end
        rec.ind = struct('arrKm', sqrt(sum((yh(1:3)-cc.rvf(1:3)').^2))*kmPerNd, ...
            'arrMs', sqrt(sum((yh(4:6)-cc.rvf(4:6)').^2))*msPerNd, ...
            'J', trapz(T, sT.^2), 'mf', yh(7), 'lamMf', yh(14), ...
            'HdriftAbs', max(abs(HT - HT(1))), ...
            'HdriftRel', max(abs(HT - HT(1)))/max(abs(HT(1)),1e-12), 'H', HT(1), ...
            'sMin', min(sT), 'sMax', max(sT), 'T', T, 's', sT);
        rec.costRelDiff = abs(rec.ind.J - o.J)/o.J;
        fprintf('  indirect flight: arrival %.3f km / %.4f m/s, J = %.6f (direct %.6f, rel %.1e), mf = %.6f, H = %.6f (drift %.1e), s in [%.3f, %.3f]\n', ...
                rec.ind.arrKm, rec.ind.arrMs, rec.ind.J, o.J, rec.costRelDiff, rec.ind.mf, ...
                rec.ind.H, rec.ind.HdriftAbs, rec.ind.sMin, rec.ind.sMax);

        % --- verdict ------------------------------------------------------
        gates = [rec.direct.flownKm < P.globTolKm && rec.direct.flownMs < P.globTolMs, ...
                 rec.ms.converged, rec.ms.accepted, rec.ind.HdriftAbs < 1e-6, ...
                 rec.ind.arrKm < P.globTolKm && rec.ind.arrMs < P.globTolMs, ...
                 rec.costRelDiff < 1e-4, rec.ind.sMax < 0.999 || rec.ind.sMin < 0.999];
        rec.gates = gates;
        if all(gates), rec.verdict = 'PASS'; else, rec.verdict = 'FAIL'; end
        fprintf('  gates [flown ms accept H arrive cost interior] = %s -> %s\n', ...
                sprintf('%d', gates), rec.verdict);
        R = append(R, rec);  save(outMat, 'R');
        make_figure(rec, o, resDir, tStar);
    end
end
report(R);
end

% ---------------------------------------------------------------------------
function R = append(R, rec)
% APPEND  Add a record, replacing an existing (iD,iA,gam) entry.
% INPUTS: R struct array; rec struct. OUTPUTS: R.
if isempty(R), R = rec; return, end
hit = find([R.iD]==rec.iD & [R.iA]==rec.iA & abs([R.gam]-rec.gam)<1e-12, 1);
% field sets may differ (failed direct solves carry fewer); merge by name
if isempty(hit), hit = numel(R) + 1; end
fn = fieldnames(rec);
for k = 1:numel(fn), R(hit).(fn{k}) = rec.(fn{k}); end
end

% ---------------------------------------------------------------------------
function report(R)
% REPORT  One line per record.
if isempty(R), fprintf('no records\n'); return, end
fprintf('\n%-6s %-5s %-8s %-8s %-9s %-9s %-7s %-6s %-8s %-8s %-8s %-6s\n', ...
    'cell', 'gam', 'tf[d]', 'J', 'mf_dir', 'mf_ind', 'thrMin', 'K', '||R||', 'acc dz', 'Hdrift', 'verdict');
for k = 1:numel(R)
    r = R(k);
    if isfield(r, 'ms') && ~isempty(r.ms)
        if ~isfield(r.ind, 'HdriftAbs'), r.ind.HdriftAbs = NaN; end   % pre-regate record
        fprintf('(%2d,%2d) %-5.2f %-8.3f %-8.5f %-9.6f %-9.6f %-7.3f %-6d %-8.1e %-8.1e %-8.1e %s\n', ...
            r.iD, r.iA, r.gam, r.tf*382981.289129055/86400, r.direct.J, r.direct.mf, r.ind.mf, ...
            r.direct.thrMin, r.ms.K, r.ms.normR, r.ms.acceptDz, r.ind.HdriftAbs, r.verdict);
    else
        fprintf('(%2d,%2d) %-5.2f %s\n', r.iD, r.iA, r.gam, r.verdict);
    end
end
end

% ---------------------------------------------------------------------------
function make_figure(rec, o, resDir, tStar)
% MAKE_FIGURE  Throttle (direct nodes vs indirect flight) + costates.
f = figure('Visible','off','Position',[100 100 1000 420]);
tD = o.tNodes*tStar/86400;  tI = rec.ind.T*tStar/86400;
subplot(1,2,1)
plot(tD, o.U(4,:), '.', 'MarkerSize', 6); hold on
plot(tI, rec.ind.s, '-', 'LineWidth', 1.2);
ylim([0 1.05]); grid on; xlabel('t [days]'); ylabel('throttle s');
title(sprintf('cell (%d,%d), \\gamma = %.2f: J_{dir} = %.4f, J_{ind} = %.4f', ...
      rec.iD, rec.iA, rec.gam, rec.direct.J, rec.ind.J));
legend('direct nodes', 'indirect flight', 'Location', 'best');
subplot(1,2,2)
lam = rec.lamDef(1:7,:);  tM = 0.5*(o.tNodes(1:end-1)+o.tNodes(2:end))*tStar/86400;
plot(tM, rec.harvest.sign*lam(1:3,:), '-'); hold on
plot(tM, rec.harvest.sign*lam(4:6,:), '--'); plot(tM, rec.harvest.sign*lam(7,:), ':', 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('costates (defect duals, sign-resolved)');
title(sprintf('\\lambda_0 (ms) = [%.3g %.3g %.3g | %.3g %.3g %.3g | %.3g]', rec.ms.z));
print(f, fullfile(resDir, sprintf('minenergy_%d_%d_g%.2f.png', rec.iD, rec.iA, rec.gam)), '-dpng', '-r110');
close(f);
end
