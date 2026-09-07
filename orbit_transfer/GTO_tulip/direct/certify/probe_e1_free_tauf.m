function res = probe_e1_free_tauf(opts)
% PROBE_E1_FREE_TAUF  Does freeing tau_f move the certified min-fuel optimum?
%
% Motivation (GPT-6 Astra code review of the direct core chain, finding E1,
% 2026-09-06, doc/reviews/direct_core_chain_gpt6astra_review_2026-09-06.md):
% casadi_minfuel_sundman holds the regularized length tau_f FIXED (= tauf0 from
% the seed) AND pins t(tau_f) = t_f. Every physical trajectory has its own
% regularized length tau_f = int_0^tf dt/kappa(x), so fixing it adds the
% isoperimetric constraint int dt/kappa = tauf0 to the fixed-time problem. The
% campaign's P0 run (TODO.md, 2026-07-26) already showed the constraint binds:
% freeing tau_f via the cScale slack at the flagship t_f moved cScale to 1.0051.
% What is not known is the size of the effect on m_f at eps = 0.
%
% This probe takes stored fixed-tau_f solutions, re-solves each with the
% campaign engine at eps = 0 (like-for-like baseline in THIS process/mesh),
% then re-solves with the free-time solver casadi_energy_freetf (single-primary
% clock, moonZone <= 0; tau_f free through cScale; t(tau_f) = t_f PINNED) and
% reports the shift in final mass.
%
% PRE-REGISTERED READING of |dm_f| = |m_f(free) - m_f(fixed re-solve)|:
%   < 1e-7          E1 is a formal defect with no practical consequence here
%   1e-7 .. 1e-6    ambiguous; report both numbers, no headline change
%   > 1e-6 (15 mg)  a REAL shift: the certified rows are extremals of the
%                   restricted problem and published m_f values need free-tau_f
%                   re-solves (compare against the basin spreads, 1e-3 .. 2e-3)
%
% INPUTS:
%   opts - struct (all optional):
%     .files    - cellstr of fixed-tau_f solution .mats (fields out.X [8xN+1],
%                 out.U [4xN+1], sigma, tauf0, rv0, rvf)    [flagship + basin24]
%     .outDir   - results folder                            [results/e1_freetauf]
%     .logFile  - append-mode log                           [outDir/probe.log]
%     .maxIter  - IPOPT cap per solve                       [1500]
%     .sched    - fallback eps schedule if the direct eps=0 free re-solve does
%                 not converge tight                        [cfg.schedNeighbor]
%
% OUTPUTS:
%   res - struct array, one per file: .file .mf_stored .mf_fixed .mf_free
%         .cScale .dmf .dprop_g .sw_fixed .sw_free .defect_fixed .defect_free
%         .status_fixed .status_free .route ('direct'|'schedule'|'FAILED')
%         Full solver outputs saved per file to outDir/probe_e1_<base>.mat.
%
% REFERENCES:
%   [1] doc/reviews/direct_core_chain_gpt6astra_review_2026-09-06.md (E1, E4).
%   [2] ../TODO.md "Free-span reformulation", P0 (cScale = 1.0051 at eps = 1).
%   [3] GTO_ELFO/direct/elfo/casadi_energy_freetf.m (Betts cScale slack state).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts, f, d);
here   = fileparts(mfilename('fullpath'));
libDir = fullfile(here, '..', 'lib');
files  = gd('files', {fullfile(libDir, 'sundman_minfuel_certified.mat'), ...
                      fullfile(libDir, 'sundman_minfuel_basin24_f1150.mat')});
outDir = gd('outDir', fullfile(here, '..', 'results', 'e1_freetauf'));
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = gd('logFile', fullfile(outDir, 'probe.log'));
maxIter = gd('maxIter', 1500);
cfg     = minfuel_config();
sched   = gd('sched', cfg.schedNeighbor);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
lg('=== PROBE E1 free tau_f  %s ===', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
lg('thrust %.4g N  m0 %g kg  Isp %g s  pSund %g  maxIter %d', cfg.thrustN, cfg.m0kg, cfg.ispS, cfg.pSund, maxIter);

res = struct([]);
for k = 1:numel(files)
    f = files{k};  [~, base] = fileparts(f);
    lg('--- [%d/%d] %s', k, numel(files), base);
    R = load(f);
    if isfield(R, 'out'), X = R.out.X(1:8,:);  U = R.out.U; else, X = R.X(1:8,:); U = R.U; end
    sigma = R.sigma(:);  tauf0 = R.tauf0;  rv0 = R.rv0;  rvf = R.rvf;
    tf = X(8, end);
    mfStored = X(7, end);
    lg('    stored: N=%d  tf=%.6f ND  tauf0=%.6f  m_f=%.8f', numel(sigma)-1, tf, tauf0, mfStored);

    r = struct('file', f, 'mf_stored', mfStored, 'mf_fixed', NaN, 'mf_free', NaN, ...
               'cScale', NaN, 'dmf', NaN, 'dprop_g', NaN, 'sw_fixed', NaN, 'sw_free', NaN, ...
               'defect_fixed', NaN, 'defect_free', NaN, 'status_fixed', '', 'status_free', '', ...
               'route', 'FAILED', 'tf', tf, 'tauf0', tauf0);

    % --- (1) like-for-like baseline: fixed-tau_f engine, eps = 0, tight -----
    tA = tic;
    oF = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                                X, U, tauf0, cfg.pSund, maxIter, 0, true);
    r.mf_fixed = oF.mf;  r.sw_fixed = oF.switches;  r.defect_fixed = oF.maxDefect;
    r.status_fixed = oF.ipoptStatus;
    lg('    FIXED re-solve: %s  defect=%.2e  m_f=%.8f  sw=%d  (%.0f s)', ...
       oF.ipoptStatus, oF.maxDefect, oF.mf, oF.switches, toc(tA));
    okF = strcmp(oF.ipoptStatus, 'Solve_Succeeded') && oF.maxDefect < 1e-6;
    if okF, X = oF.X;  U = oF.U; end
    save(fullfile(outDir, ['probe_e1_' base '.mat']), 'r', 'oF', 'sigma', 'tauf0', 'rv0', 'rvf', 'tf');

    % --- (2) free tau_f: cScale slack, t_f pinned, eps = 0 tight from (1) ----
    o = struct('epsilon', 0, 'tfTarget', tf, 'moonZone', -1, 'pSund', cfg.pSund, ...
               'maxIter', maxIter, 'warmTight', true, 'c0', 1);
    tB = tic;
    oE = casadi_energy_freetf(sigma, rv0, rvf, p.Tmax, p.c, p.muStar, X, U, tauf0, o);
    lg('    FREE direct eps=0: %s  defect=%.2e  m_f=%.8f  cScale=%.6f  sw=%d  (%.0f s)', ...
       oE.ipoptStatus, oE.maxDefect, oE.mf, oE.cScale, oE.switches, toc(tB));
    okE = oE.success && oE.maxDefect < 1e-6;
    route = 'direct';
    if ~okE
        % fallback: short eps schedule, first step loose (a genuine move), rest tight
        lg('    direct eps=0 not tight -> schedule fallback %s', mat2str(sched));
        Xk = X;  Uk = U;  okE = false;  route = 'schedule';
        for ke = 1:numel(sched)
            o.epsilon = sched(ke);  o.warmTight = ke > 1;
            oS = casadi_energy_freetf(sigma, rv0, rvf, p.Tmax, p.c, p.muStar, Xk, Uk, tauf0, o);
            okS = oS.success && oS.maxDefect < 1e-6;
            lg('      eps=%.4g: %s defect=%.2e m_f=%.8f cScale=%.6f sw=%d', ...
               sched(ke), oS.ipoptStatus, oS.maxDefect, oS.mf, oS.cScale, oS.switches);
            if okS, Xk = oS.X;  Uk = oS.U;  oE = oS;  okE = (sched(ke) == 0); end
            if ~okS, break; end
        end
    end
    if okE
        r.mf_free = oE.mf;  r.cScale = oE.cScale;  r.sw_free = oE.switches;
        r.defect_free = oE.maxDefect;  r.status_free = oE.ipoptStatus;  r.route = route;
        r.dmf = oE.mf - r.mf_fixed;  r.dprop_g = -1000 * cfg.m0kg * r.dmf;
        lg('    RESULT: m_f fixed %.8f -> free %.8f   dm_f = %+.3e (%+.2f g propellant)   cScale = %.6f', ...
           r.mf_fixed, r.mf_free, r.dmf, r.dprop_g, r.cScale);
        if abs(r.dmf) < 1e-7,     lg('    READING: |dm_f| < 1e-7 -> E1 formal only at this t_f');
        elseif abs(r.dmf) < 1e-6, lg('    READING: 1e-7 <= |dm_f| < 1e-6 -> ambiguous');
        else,                     lg('    READING: |dm_f| >= 1e-6 -> REAL shift; certified rows are restricted-problem extremals');
        end
    else
        r.status_free = oE.ipoptStatus;  r.route = 'FAILED';
        lg('    FREE re-solve FAILED to converge tight (status %s)', oE.ipoptStatus);
    end
    save(fullfile(outDir, ['probe_e1_' base '.mat']), 'r', 'oF', 'oE', 'sigma', 'tauf0', 'rv0', 'rvf', 'tf');
    if isempty(res), res = r; else, res(end+1) = r; end %#ok<AGROW>
end

lg('=== SUMMARY ===');
lg('%-40s %12s %12s %12s %10s %8s %8s', 'file', 'mf_fixed', 'mf_free', 'dmf', 'dprop_g', 'cScale', 'route');
for k = 1:numel(res)
    [~, b] = fileparts(res(k).file);
    lg('%-40s %12.8f %12.8f %+12.3e %+10.2f %8.5f %8s', b, res(k).mf_fixed, res(k).mf_free, ...
       res(k).dmf, res(k).dprop_g, res(k).cScale, res(k).route);
end
save(fullfile(outDir, 'probe_e1_summary.mat'), 'res');
lg('PROBE E1 DONE');
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
% GETDEF  s.(f) if present and nonempty, else d.
% INPUTS: s [struct]; f [char]; d [any].   OUTPUTS: v [any].   REFERENCES: none.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% ---------------------------------------------------------------------------
function logmsg(f, s)
% LOGMSG  Append one line to the log file (or stdout if no file).
% INPUTS: f [char path or '']; s [char].   OUTPUTS: none.   REFERENCES: matlab-campaign skill.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
