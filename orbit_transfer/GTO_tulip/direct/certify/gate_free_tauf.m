function res = gate_free_tauf(opts)
% GATE_FREE_TAUF  First-order gate + results-set export for the free-tau_f rows.
%
% The 2026-09-06/07 probe (probe_e1_free_tauf) re-solved every stored row of
% the dV-t_f front with tau_f released through a cScale slack state at the
% same pinned t_f, and found the fixed-tau_f rows restricted by up to
% 9.4e-3 in m_f. Those free rows were converged IPOPT optima with clean raw
% duals, but they had not been through the generic first-order gate, and
% they lived only as probe artifacts. This driver does both, per row:
%
%   1. RE-SOLVE with the campaign engine in its ported free-tau_f mode
%      (casadi_minfuel_sundman, opts.freeTauf + returnModel), warm-started AT
%      the probe's 9-row free solution WITH its multipliers (opts.lamG0), eps
%      = 0, tight. This is also the acceptance test of the port: the engine
%      must sit still at the number casadi_energy_freetf produced. (The first
%      pass, 2026-09-07 08:21, ran without the multipliers: 10/17 rows sat
%      still anyway, 7 wandered off within 800 iterations -- the zero-dual
%      restart, see the engine header; those 7 were re-run with lamG0.)
%   2. GUARD (verify_common/certified_guard): converged, machine-tight, m_f
%      not worse than the probe's.
%   3. GATE (verify_common/foc_check, manifest 'tulip_free', nx = 9): full
%      KKT stationarity, bang-bang sign law, signed primer minimum condition,
%      Hager-mapped free-mass transversality, singular-arc / regular-switching
%      diagnostics, IPOPT inertia verdict; standard report + sidecar via
%      foc_report into certify/results/.
%   4. EXPORT a results-set row in the campaign's minfuel_at_tf layout to
%      results/minfuel_freetauf/minfuel_f####_free*.mat: 8-row out.X/out.U/
%      out.lamDef (the slack is a constant, exposed as out.cScale, with the
%      full 9-row primal/costates in out.X9/out.lamDef9), out.factor .tf
%      .tf_days .dV .prop_kg .certified .meta, and top-level sigma, tauf0
%      (= the EFFECTIVE regularized length cScale*tauf0_seed -- a free-tau_f
%      solution IS a fixed-tau_f KKT point at that length, so either engine
%      mode re-solves it in place), rv0, rvf, factor, fp. aggregate_front
%      reads this directory beside results/minfuel/.
%
% INPUTS:
%   opts - struct (all optional):
%     .files   - cellstr of probe_e1_<base>.mat files  [all in results/e1_freetauf]
%     .outDir  - results-set folder                    [cfg.dirs.minfuelFree]
%     .logFile - append-mode log                       [outDir/gate.log]
%     .maxIter - IPOPT cap for the sit-still re-solve  [800]
%     .resume  - true: skip rows whose results-set file exists   [false]
%     .save    - write the results-set rows            [true]
%
% OUTPUTS:
%   res - struct array, one per row: .base .factor .mf_probe .mf .cScale
%         .dmf_vs_probe .switches .maxDefect .ipoptStatus .kktStatInf
%         .signPct .dirSignedMax .dirSignedPct .lamMassEndMapped .sdotMinRel .lamTimeCoV
%         .focPass .inertia .guard ('OK'|'IMPROVED'|message) .outFile .status
%         ('gated'|'FAILED'). Also written to outDir/gate_summary.{mat,txt}.
%
% REFERENCES:
%   [1] probe_e1_free_tauf.m (the free rows; OPTIMALITY_CERTIFICATION.md
%       LEAD-5 FRONT RESULT for the numbers).
%   [2] verify_common/foc_check.m, foc_manifest.m ('tulip_free'), foc_report.m,
%       foc_ipopt_inertia.m, certified_guard.m.
%   [3] lib/casadi_minfuel_sundman.m header (freeTauf port, 2026-09-07) and
%       lib/minfuel_at_tf.m (the results-row layout this reproduces).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts, f, d);
here = fileparts(mfilename('fullpath'));
addpath(here);  addpath(fullfile(here, '..', 'lib'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
vcDir = fullfile(here, '..', '..', '..', 'verify_common');
addpath(vcDir);  setup_verify_common();
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

probeDir = fullfile(here, '..', 'results', 'e1_freetauf');
files   = gd('files', default_files(probeDir));
outDir  = gd('outDir', cfg.dirs.minfuelFree);
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = gd('logFile', fullfile(outDir, 'gate.log'));
maxIter = gd('maxIter', 800);
resume  = gd('resume', false);
doSave  = gd('save', true);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
lg('=== GATE free tau_f  %s ===', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
lg('engine casadi_minfuel_sundman freeTauf=true  manifest tulip_free  maxIter %d  rows %d', maxIter, numel(files));

res = struct([]);
for k = 1:numel(files)
    f = files{k};  [~, base] = fileparts(f);
    base = regexprep(base, '^probe_e1_', '');
    [factor, outFile] = row_name(base, cfg, outDir);
    r = struct('base', base, 'factor', factor, 'mf_probe', NaN, 'mf', NaN, 'cScale', NaN, ...
        'dmf_vs_probe', NaN, 'switches', NaN, 'maxDefect', NaN, 'ipoptStatus', '', ...
        'kktStatInf', NaN, 'signPct', NaN, 'dirSignedMax', NaN, 'dirSignedPct', NaN, 'lamMassEndMapped', NaN, ...
        'sdotMinRel', NaN, 'lamTimeCoV', NaN, 'focPass', false, 'inertia', '', ...
        'guard', '', 'outFile', outFile, 'status', 'FAILED');
    if resume && isfile(outFile)
        lg('--- [%d/%d] %s  (resume: %s exists, skipping)', k, numel(files), base, outFile);
        R = load(outFile, 'out');
        r.mf = R.out.mf;  r.cScale = R.out.cScale;  r.switches = R.out.switches;
        r.maxDefect = R.out.maxDefect;  r.status = 'gated';
        if isfield(R.out.meta, 'foc')
            fo = R.out.meta.foc;
            r.kktStatInf = fo.kktStatInf;  r.signPct = fo.signPct;  r.dirSignedMax = fo.dirSignedMax;
            r.dirSignedPct = fo.dirSignedPct;
            r.lamMassEndMapped = fo.lamMassEndMapped;  r.sdotMinRel = fo.sdotMinRel;
            r.lamTimeCoV = fo.lamTimeCoV;  r.focPass = fo.pass;  r.inertia = fo.inertia;
        end
        res = append_row(res, r);
        continue
    end
    lg('--- [%d/%d] %s  (factor %.3f)', k, numel(files), base, factor);
    P = load(f);
    if ~isfield(P, 'oE') || ~isfield(P.oE, 'X') || size(P.oE.X,1) ~= 9
        lg('    no 9-row free solution in the probe file -- skipped');
        res = append_row(res, r);  continue
    end
    sigma = P.sigma(:);  tauf0 = P.tauf0;  rv0 = P.rv0;  rvf = P.rvf;  tf = P.tf;
    r.mf_probe = P.oE.mf;
    lg('    probe: m_f %.8f cScale %.6f sw %d tf %.6f tauf0(seed) %.6f', ...
       P.oE.mf, P.oE.cScale, P.oE.switches, tf, tauf0);

    % --- 1. sit-still re-solve in the ported free-tau_f mode, model attached
    tA = tic;
    try
        % dual warm start from the probe's own multipliers (same NLP, same
        % constraint order as casadi_energy_freetf): without it IPOPT starts
        % from zero multipliers and 7/17 rows wandered off within 800 iterations
        lamG0 = [];  if isfield(P.oE, 'lamAll'), lamG0 = P.oE.lamAll; end
        out = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
              P.oE.X, P.oE.U, tauf0, cfg.pSund, maxIter, 0, true, ...
              struct('freeTauf', true, 'returnModel', true, 'lamG0', lamG0));
    catch err
        lg('    ENGINE ERROR: %s', err.message);
        res = append_row(res, r);  continue
    end
    r.mf = out.mf;  r.cScale = out.cScale;  r.switches = out.switches;
    r.maxDefect = out.maxDefect;  r.ipoptStatus = out.ipoptStatus;
    r.dmf_vs_probe = out.mf - P.oE.mf;
    lg('    re-solve: %s defect %.2e m_f %.8f (d vs probe %+.2e) cScale %.6f (probe %.6f) sw %d (%.0f s)', ...
       out.ipoptStatus, out.maxDefect, out.mf, r.dmf_vs_probe, out.cScale, P.oE.cScale, out.switches, toc(tA));

    % --- 2. certified-quantity guard (one-sided: not worse than the probe)
    try
        gi = certified_guard( ...
            struct('success', out.success, 'ipoptStatus', out.ipoptStatus, ...
                   'maxDefect', out.maxDefect, 'value', out.mf), ...
            struct('caller', 'gate_free_tauf', 'label', base, 'saved', P.oE.mf, ...
                   'name', 'm_f', 'errName', 'mass', 'better', 'higher', ...
                   'feasTol', 1e-8, 'tol', 1e-6));
        if isfield(gi, 'improved') && gi.improved, r.guard = 'IMPROVED'; else, r.guard = 'OK'; end
    catch err
        r.guard = err.message;
        lg('    GUARD REFUSED: %s', err.message);
        res = append_row(res, r);  continue
    end

    % --- 3. generic first-order gate on the 9-state model
    tB = tic;
    out9 = out;  out9.X = out.X9;
    rep = foc_check(out9, sigma, foc_manifest('tulip_free'), struct('eps', 0));
    rep.ipopt = foc_ipopt_inertia(out.regHistory);
    r.kktStatInf = rep.kktStatInf;  r.signPct = rep.signPct;
    r.dirSignedMax = rep.dirSignedMax;  r.dirSignedPct = rep.dirSignedPct;
    r.lamMassEndMapped = rep.lamMassEndMapped;
    r.sdotMinRel = rep.sdotMinRel;  r.lamTimeCoV = rep.lamTimeCoV;  r.focPass = rep.pass;
    if isfield(rep.ipopt, 'verdict'), r.inertia = char(rep.ipopt.verdict); end
    lg('    foc_check: kktStat %.2e sLag %+d signPct %.2f%% dirTan %.2e dirSigned max %.2e pct %.2f%% lamMassMapped %.2e (fullGL %.2e) lamTimeCoV %.3f sdotMinRel %.2e (phys %.2e) singular %d sw %d PASS=%d inertia=%s (%.0f s)', ...
       rep.kktStatInf, rep.sLag, rep.signPct, rep.dirTanMax, rep.dirSignedMax, rep.dirSignedPct, ...
       rep.lamMassEndMapped, rep.lamMassEndFullGL, rep.lamTimeCoV, rep.sdotMinRel, ...
       rep.sdotMinRelPhys, rep.singularArcNodes, rep.nSwitches, rep.pass, r.inertia, toc(tB));
    tag = ['free_' base];
    foc_report(rep, tag, fullfile(here, 'results'));

    % --- 4. results-set row (minfuel_at_tf layout, 8-row contract + slack)
    focSummary = struct('kktStatInf', rep.kktStatInf, 'signPct', rep.signPct, ...
        'dirTanMax', rep.dirTanMax, 'dirSignedMax', rep.dirSignedMax, 'dirSignedPct', rep.dirSignedPct, ...
        'lamMassEndMapped', rep.lamMassEndMapped, 'lamMassEndFullGL', rep.lamMassEndFullGL, ...
        'lamTimeCoV', rep.lamTimeCoV, 'sdotMinRel', rep.sdotMinRel, ...
        'sdotMinRelPhys', rep.sdotMinRelPhys, 'singularArcNodes', rep.singularArcNodes, ...
        'nSwitches', rep.nSwitches, 'pass', rep.pass, 'inertia', r.inertia, ...
        'sidecar', tag, 'manifest', 'tulip_free');
    row = out;
    row.certified  = true;              % converged tight at eps = 0 (the row's definition)
    row.epsReached = 0;
    row.factor  = factor;  row.tf = tf;  row.tf_days = tf*p.tStar/86400;
    row.dV      = p.c*log(1/out.mf)*p.lStar/p.tStar;
    row.prop_kg = p.m0kg*(1-out.mf);
    row.meta = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
        'githash', git_hash(here), 'seed', ['probe_e1_' base '.mat (free-tau_f re-solve of ' base ')'], ...
        'branch', 'free', 'sched', 0, 'maxIter', maxIter, 'pSund', cfg.pSund, ...
        'tfMin', cfg.tfMin, 'freeTauf', true, 'taufSeed', tauf0, 'cScale', out.cScale, ...
        'taufEffective', out.tauf, 'mfFixedTauf', P.r.mf_fixed, 'foc', focSummary, ...
        'guard', r.guard, 'dualWarmStart', ~isempty(lamG0), ...
        'solver', 'casadi_minfuel_sundman freeTauf=true (CasADi+IPOPT, Sundman trapezoid, cScale slack)');
    if isfield(row, 'model'), row = rmfield(row, 'model'); end   % never persist the Opti object
    if doSave
        % top-level tauf0 is the EFFECTIVE length cScale*tauf0_seed (see header)
        A = struct('out', row, 'sigma', sigma, 'tauf0', out.tauf, 'rv0', rv0, 'rvf', rvf, ...
                   'factor', factor, ...
                   'fp', cr3bp_fingerprint(p, struct('tf', tf, 'factor', factor, 'freeTauf', true)));
        save(outFile, '-struct', 'A');
        lg('    WROTE %s  (tauf0 stored = effective %.6f)', outFile, out.tauf);
    end
    r.status = 'gated';
    res = append_row(res, r);
end

% --- summary --------------------------------------------------------------
save(fullfile(outDir, 'gate_summary.mat'), 'res');
fid = fopen(fullfile(outDir, 'gate_summary.txt'), 'w');
fprintf(fid, '%-34s %6s %11s %11s %9s %8s %4s %8s %8s %7s %7s %8s %8s %7s %5s %-14s %s\n', ...
    'row', 'factor', 'mf_probe', 'mf', 'd_probe', 'cScale', 'sw', 'defect', 'kktStat', 'sign%', ...
    'sgn%', 'lamMass', 'sdotRel', 'ltCoV', 'pass', 'inertia', 'status');
for k = 1:numel(res)
    q = res(k);
    fprintf(fid, '%-34s %6.3f %11.8f %11.8f %+9.1e %8.5f %4d %8.1e %8.1e %6.2f%% %6.2f%% %8.1e %8.1e %7.3f %5d %-14s %s\n', ...
        q.base, q.factor, q.mf_probe, q.mf, q.dmf_vs_probe, q.cScale, q.switches, q.maxDefect, ...
        q.kktStatInf, q.signPct, q.dirSignedPct, q.lamMassEndMapped, q.sdotMinRel, q.lamTimeCoV, ...
        q.focPass, q.inertia, q.status);
end
fclose(fid);
lg('summary -> %s', fullfile(outDir, 'gate_summary.txt'));
lg('GATE_ALL_DONE  gated %d / %d', sum(strcmp({res.status}, 'gated')), numel(res));
end

% ---------------------------------------------------------------------------
function files = default_files(probeDir)
% DEFAULT_FILES  Every per-row probe artifact in the free-tau_f probe folder.
% INPUTS: probeDir [char].  OUTPUTS: files [cellstr].  REFERENCES: none.
d = dir(fullfile(probeDir, 'probe_e1_*.mat'));
d = d(~strcmp({d.name}, 'probe_e1_summary.mat'));      % the probe's own summary is not a row
files = fullfile({d.folder}, {d.name});
if isempty(files), error('gate_free_tauf:noProbes', 'no probe_e1_*.mat in %s', probeDir); end
end

function [factor, outFile] = row_name(base, cfg, outDir)
% ROW_NAME  Factor + results-set filename for a probe row.
% Two rows sit at 1.150 (the published flagship and the basin-24 winner), so
% those carry their origin in the suffix; every other row is minfuel_f####_free.
% INPUTS: base [char]; cfg [struct]; outDir [char].
% OUTPUTS: factor [scalar]; outFile [char].  REFERENCES: none.
tok = regexp(base, '_f(\d{4})', 'tokens', 'once');
if isempty(tok), factor = 1.150; else, factor = str2double(tok{1})/1000; end
name = strrep(cfg.fname('minfuel', factor), '.mat', '_free.mat');
if contains(base, 'certified'), name = strrep(name, '_free.mat', '_free_flagship.mat'); end
if contains(base, 'basin24'),   name = strrep(name, '_free.mat', '_free_basin24.mat');  end
outFile = fullfile(outDir, name);
end

function res = append_row(res, r)
% APPEND_ROW  Grow the result struct array.
% INPUTS: res [struct array]; r [struct].  OUTPUTS: res.  REFERENCES: none.
if isempty(res), res = r; else, res(end+1) = r; end
end

function h = git_hash(here)
% GIT_HASH  Short git hash for provenance (empty on failure -- never blocks).
% INPUTS: here [char].  OUTPUTS: h [char].  REFERENCES: none.
[rc, s] = system(sprintf('cd "%s" && git rev-parse --short HEAD 2>/dev/null', here));
if rc == 0, h = strtrim(s); else, h = ''; end
end

function v = getdef(s, f, d)
% GETDEF  s.(f) if present and nonempty, else d.
% INPUTS: s [struct]; f [char]; d [any].  OUTPUTS: v [any].  REFERENCES: none.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function logmsg(f, s)
% LOGMSG  Append one line to the log file (and echo to stdout).
% INPUTS: f [char]; s [char].  OUTPUTS: none.  REFERENCES: matlab-campaign skill.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); fprintf('%s\n', s); end
end
