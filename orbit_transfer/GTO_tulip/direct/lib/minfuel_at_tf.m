function out = minfuel_at_tf(factor, varargin)
% MINFUEL_AT_TF  Canonical single-t_f min-fuel solve (THE per-t_f driver).
%
% Consolidates the four superseded drivers (solve_tf_minfuel, tf_step, and the
% step logic copy-pasted inside run_tf_front / run_tf_2anchor) into one
% function with an explicit seed source, one schedule policy from
% MINFUEL_CONFIG, and provenance-stamped output files. Recipe per seed type:
%
%   'energy'   (default) seed = energy backbone file at THIS factor:
%              (1) TIGHT re-clean at eps=1 (loose-continued backbone energy has
%                  inconsistent KKT duals; sharpening it directly blows up
%                  inf_du ~1e10), (2) fine energy->fuel sharpen, warmTight,
%                  schedule cfg.schedSharpen (ends at exactly eps=0).
%   'neighbor' seed = an existing bang-bang solution at opts.seedFactor:
%              rescale the time state to this t_f, FIRST step loose (genuine
%              continuation move), remaining steps tight, cfg.schedNeighbor.
%   <path>     seed = explicit .mat with X,U (top level or in `out`), treated
%              like 'neighbor'.
%
% FREE tau_f (default since 2026-09-07). The engine is run with
% opts.freeTauf = true: the regularized length is released through a cScale
% slack state while t(tau_f) = t_f stays pinned, so the NLP IS the fixed-time
% problem. Before this every row inherited its seed's tau_f0 unchanged -- the
% neighbour path rescaled the time state but never tau_f0, and the energy
% backbones carried the original 1.150x length -- which added the
% isoperimetric constraint int dt/kappa = tau_f0(seed) to every solve and cost
% up to 9.4e-3 in m_f across the front (review E1; OPTIMALITY_CERTIFICATION.md
% LEAD-5). The stored top-level tauf0 of a free row is the EFFECTIVE length
% cScale*tau_f0(seed), so a later neighbour seed starts from a consistent
% value, and free rows are written to cfg.dirs.minfuelFree with '_free' in
% the name (results/minfuel/ keeps the fixed-tau_f originals). Pass
% 'freeTauf', false to reproduce the pre-2026-09-07 fixed-tau_f behaviour
% (a warning then names the inherited tau_f0).
%
% INPUTS:
%   factor  - t_f / t_f^min [scalar]
%   options (name-value):
%     'seed'       - 'energy' | 'neighbor' | file path      [default 'energy']
%     'seedFactor' - factor of the neighbor solution        [required for 'neighbor']
%     'sched'      - homotopy epsilon schedule override     [default by seed type]
%     'maxIter'    - IPOPT iteration cap                    [default cfg.maxIter]
%     'thrustN'    - max thrust [N] for an off-nominal (ladder) rung
%                    [default cfg.thrustN = campaign nominal]. Off-nominal
%                    requires an explicit seed file: the energy backbones are
%                    nominal-thrust, and the fingerprint check will refuse a
%                    mismatched seed rather than silently warm-start wrong.
%     'branch'     - branch tag recorded in meta + filename suffix, e.g.
%                    'up','dn','en'                          [default 'en'|'nb']
%     'freeTauf'   - release tau_f through the engine's cScale slack
%                    (see above)                             [default true]
%     'outFile'    - output path override                   [default results/minfuel/]
%     'save'       - write the .mat                          [default true]
%
% OUTPUTS:
%   out - solver struct (X,U,lamDef,switches,edge,maxDefect,primerAlignDeg,...;
%         in free-tau_f mode also .cScale .X9 .lamDef9 .taufSeed, with .tauf
%         the effective length) plus .factor .tf .tf_days .dV .prop_kg
%         .certified (logical: at least one schedule step converged tight)
%         .meta (provenance: date, git hash, seed source, schedule, per-step
%         table, ipopt statuses, freeTauf, cScale, taufSeed). An uncertified
%         result is NOT saved -- a loose iterate must never become a neighbor seed.
%
% REFERENCES:
%   [1] LOW_THRUST_MINFUEL_CAMPAIGN.md ("Down-sweep CRACKED": backbone+sharpen).
%   [2] CODE_CLEANUP_PLAN.md (driver consolidation rationale).
%   [3] casadi_minfuel_sundman.m header (freeTauf port) and
%       certify/probe_e1_free_tauf.m (the measurement behind the default).

here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));  % cr3bp_fingerprint/check_cr3bp_fp
cfg  = minfuel_config();
op   = parse_opts(varargin, cfg);
% Off-nominal thrust (ladder rungs). Default is the campaign's nominal value, so
% omitting 'thrustN' is byte-identical to the pre-2026-07-27 behaviour. Only
% thrustN is overridden -- schedules, dirs and filename rules stay nominal.
% NOTE the energy-backbone seed path is nominal-thrust by construction; an
% off-nominal request must pass an explicit seed file (a ladder rung). The
% fingerprint check below enforces that rather than trusting the caller.
if abs(op.thrustN - cfg.thrustN) > 1e-15
    cfg = minfuel_config(struct('thrustN', op.thrustN));
end
p    = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
tf   = factor * cfg.tfMin;
fp   = cr3bp_fingerprint(p, struct('tf', tf, 'factor', factor));

% --- resolve the seed -------------------------------------------------------
switch op.seedKind
    case 'energy'
        efile = find_energy_file(here, cfg, factor);
        E = load(efile);
        check_cr3bp_fp(E, fp, efile, 'energy-seed');
        sigma=E.sigma; rv0=E.rv0; rvf=E.rvf; tauf0=E.tauf0;
        Xk=E.X; Uk=E.U; firstLoose=false; needClean=true;
        seedDesc = sprintf('energy backbone f=%.3f', factor);
    otherwise   % 'neighbor' or explicit file
        [S, seedDesc, seedFn] = load_seed_file(here, cfg, op);
        check_cr3bp_fp(S, fp, seedFn, 'neighbor-seed');
        sigma=S.sigma; rv0=S.rv0; rvf=S.rvf; tauf0=S.tauf0;
        Xk=S.X(1:8,:); Uk=S.U;                  % 8-row contract (free rows carry the slack separately)
        Xk(8,:) = Xk(8,:) * (tf / Xk(8,end));   % rescale time state to new t_f
        firstLoose=true; needClean=false;       % first sched step is the move
        if ~op.freeTauf
            warning('minfuel_at_tf:inheritedTauf0', ...
                ['fixed-tau_f mode: the neighbour seed''s regularized length tau_f0 = %.6f is ' ...
                 'inherited UNCHANGED at the new t_f -- this is the isoperimetric restriction ' ...
                 'of review E1 (LEAD-5); pass ''freeTauf'', true (the default) to release it'], tauf0);
        end
end
engOpts = struct('freeTauf', op.freeTauf);

fprintf('MINFUEL_AT_TF: factor=%.3f  t_f=%.4f ND (%.2f d)  seed=%s\n', ...
        factor, tf, tf*p.tStar/86400, seedDesc);

% --- (1) tight re-clean of an energy seed (same t_f -> no wedge) ------------
stat = {};
if needClean
    oT = casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar, ...
                                Xk,Uk,tauf0,cfg.pSund,op.maxIter,1,true,engOpts);
    fprintf('  re-clean energy: ok=%d defect=%.2g\n', oT.success, oT.maxDefect);
    stat{end+1} = sprintf('reclean:%s', oT.ipoptStatus);
    if strcmp(oT.ipoptStatus,'Solve_Succeeded') && oT.maxDefect < 1e-6, Xk=oT.X; Uk=oT.U; end   % triage C1: full convergence only
end

% --- (2) homotopy sharpen ---------------------------------------------------
best = [];  o = [];  bestEps = NaN;  tbl = zeros(numel(op.sched), 4);
for ke = 1:numel(op.sched)
    e = op.sched(ke);
    tight = ~(firstLoose && ke==1);
    o = casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar, ...
                               Xk,Uk,tauf0,cfg.pSund,op.maxIter,e,tight,engOpts);
    ok = strcmp(o.ipoptStatus,'Solve_Succeeded') && o.maxDefect < 1e-6;   % triage C1
    tbl(ke,:) = [e, o.maxDefect, o.switches, 100*o.edge];
    stat{end+1} = sprintf('eps=%.4g:%s', e, o.ipoptStatus); %#ok<AGROW>
    if op.freeTauf
        fprintf('  eps=%.4g: ok=%d defect=%.2g sw=%d edge=%.1f%% cScale=%.6f\n', ...
                e, ok, o.maxDefect, o.switches, 100*o.edge, o.cScale);
    else
        fprintf('  eps=%.4g: ok=%d defect=%.2g sw=%d edge=%.1f%%\n', ...
                e, ok, o.maxDefect, o.switches, 100*o.edge);
    end
    if ok, Xk=o.X; Uk=o.U; best=o; bestEps=e; end
end
anyClean = ~isempty(best);
if ~anyClean
    warning('minfuel_at_tf:noCleanStep', ...
        ['no schedule step converged tight at factor %.3f; returning the last ' ...
         'UNCERTIFIED attempt (will NOT be saved)'], factor);
    best = o;  bestEps = NaN;
end
% 2026-07-21 triage C2: certification requires the REQUESTED endpoint eps
% (last schedule entry), not merely some clean step.
certified = anyClean && abs(bestEps - op.sched(end)) < 1e-12;
if anyClean && ~certified
    warning('minfuel_at_tf:endpointNotReached', ...
        ['factor %.3f stalled at eps=%.4g (requested %.4g): clean INTERMEDIATE ' ...
         'solution, NOT certified (will NOT be saved)'], factor, bestEps, op.sched(end));
end

% --- package with provenance ------------------------------------------------
out = best;
out.certified = certified;
out.epsReached = bestEps;
out.factor  = factor;  out.tf = tf;  out.tf_days = tf*p.tStar/86400;
out.dV      = p.c*log(1/best.mf)*p.lStar/p.tStar;
out.prop_kg = p.m0kg*(1-best.mf);
cScale = 1;  if op.freeTauf && isfield(best, 'cScale'), cScale = best.cScale; end
out.meta = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
    'githash', git_hash(here), 'seed', seedDesc, 'branch', op.branch, ...
    'sched', op.sched, 'maxIter', op.maxIter, 'pSund', cfg.pSund, ...
    'tfMin', cfg.tfMin, 'stepTable', tbl, 'ipoptStatuses', {stat}, ...
    'freeTauf', op.freeTauf, 'taufSeed', tauf0, 'cScale', cScale, ...
    'taufEffective', best.tauf, ...
    'solver', 'casadi_minfuel_sundman (CasADi+IPOPT, Sundman trapezoid)');
if op.freeTauf, out.meta.solver = [out.meta.solver ' freeTauf=true (cScale slack, t_f pinned)']; end

fprintf('MINFUEL_AT_TF done: f=%.3f dV=%.4f km/s sw=%d edge=%.1f%% defect=%.2g primer=%.3f\n', ...
        factor, out.dV, best.switches, 100*best.edge, best.maxDefect, best.primerAlignDeg);

if op.save && ~certified
    warning('minfuel_at_tf:skipSaveUncertified', ...
        ['factor %.3f did not converge tight; NOT writing an output file (a loose ' ...
         'iterate would poison neighbor-seed lookups). Inspect the returned struct instead.'], factor);
elseif op.save
    if isempty(op.outFile)
        base = cfg.fname('minfuel', factor);
        if op.freeTauf
            if ~exist(cfg.dirs.minfuelFree,'dir'), mkdir(cfg.dirs.minfuelFree); end
            op.outFile = fullfile(cfg.dirs.minfuelFree, strrep(base, '.mat', ['_free_' op.branch '.mat']));
        else
            if ~exist(cfg.dirs.minfuel,'dir'), mkdir(cfg.dirs.minfuel); end
            op.outFile = fullfile(cfg.dirs.minfuel, strrep(base, '.mat', ['_' op.branch '.mat']));
        end
    end
    % top-level tauf0 = the EFFECTIVE regularized length (cScale*tauf0 in free
    % mode, tauf0 itself in fixed mode), so a neighbour seed never inherits a
    % length the solution did not actually have.
    tauf0 = best.tauf;
    save(op.outFile, 'out', 'sigma', 'tauf0', 'rv0', 'rvf', 'factor', 'fp');
    fprintf('  WROTE %s  (tauf0 stored = %.6f%s)\n', op.outFile, tauf0, ternary(op.freeTauf, ' effective', ''));
end
end

% ---------------------------------------------------------------------------
function op = parse_opts(args, cfg)
% Name-value option parsing with seed-dependent defaults.
op = struct('seedKind','energy','seedFactor',NaN,'sched',[],'maxIter',cfg.maxIter, ...
            'thrustN',cfg.thrustN, 'freeTauf',true, ...
            'branch','','outFile','','save',true,'seedFile','');
for k = 1:2:numel(args)
    switch lower(args{k})
        case 'seed'
            v = args{k+1};
            if any(strcmpi(v, {'energy','neighbor'})), op.seedKind = lower(v);
            else, op.seedKind = 'file'; op.seedFile = v; end
        case 'seedfactor', op.seedFactor = args{k+1};
        case 'sched',      op.sched      = args{k+1};
        case 'maxiter',    op.maxIter    = args{k+1};
        case 'thrustn',    op.thrustN    = args{k+1};
        case 'branch',     op.branch     = args{k+1};
        case 'outfile',    op.outFile    = args{k+1};
        case 'save',       op.save       = args{k+1};
        case 'freetauf',   op.freeTauf   = logical(args{k+1});
        otherwise, error('minfuel_at_tf:badOption','unknown option %s', args{k});
    end
end
if isempty(op.sched)
    if strcmp(op.seedKind,'energy'), op.sched = cfg.schedSharpen;
    else,                            op.sched = cfg.schedNeighbor; end
end
if isempty(op.branch)
    if strcmp(op.seedKind,'energy'), op.branch = 'en'; else, op.branch = 'nb'; end
end
if strcmp(op.seedKind,'neighbor') && isnan(op.seedFactor)
    error('minfuel_at_tf:needSeedFactor','seed ''neighbor'' requires ''seedFactor''');
end
end

% ---------------------------------------------------------------------------
function f = find_energy_file(here, cfg, factor) %#ok<INUSL>
% Locate the backbone energy solution in the canonical results layout
% (results/energy/energy_f####.mat; legacy root names migrated 2026-07-09).
f = fullfile(cfg.dirs.energy, cfg.fname('energy', factor));
if ~isfile(f)
    error('minfuel_at_tf:noEnergySeed', ...
          'no energy backbone file for factor %.3f (%s); run orchestrate/backbone_walk.sh first', ...
          factor, f);
end
end

% ---------------------------------------------------------------------------
function [S, desc, fn] = load_seed_file(here, cfg, op)
% Load a bang-bang seed: explicit file, or resolve from seedFactor (new
% minfuel results, then legacy ms_<f>.mat, then the certified anchor).
if strcmp(op.seedKind, 'file')
    fn = op.seedFile;
else
    % free-tau_f rows first (results/minfuel_freetauf: gate_free_tauf rows and
    % free-mode minfuel_at_tf rows), then the fixed-tau_f originals.
    mf = cfg.fname('minfuel', op.seedFactor);
    sfx = {'_free.mat','_free_basin24.mat','_free_flagship.mat','_free_en.mat','_free_nb.mat','_free_dn.mat','_free_up.mat'};
    cand = cellfun(@(x) fullfile(cfg.dirs.minfuelFree, strrep(mf, '.mat', x)), sfx, 'UniformOutput', false);
    cand = [cand, {fullfile(cfg.dirs.minfuel, strrep(mf,'.mat','_en.mat')), ...
            fullfile(cfg.dirs.minfuel, strrep(mf,'.mat','_nb.mat')), ...
            fullfile(cfg.dirs.minfuel, cfg.fname('legacy_ms', op.seedFactor)), ...
            fullfile(here, 'sundman_minfuel_certified.mat')}];
    fn = '';
    for k = 1:numel(cand)
        if isfile(cand{k}), fn = cand{k}; break; end
    end
    if isempty(fn), error('minfuel_at_tf:noSeed','no seed found for factor %.3f', op.seedFactor); end
end
R = load(fn);
if isfield(R,'out'), S.X = R.out.X; S.U = R.out.U; else, S.X = R.X; S.U = R.U; end
S.sigma = R.sigma;  S.tauf0 = R.tauf0;  S.rv0 = R.rv0;  S.rvf = R.rvf;
if isfield(R,'fp'), S.fp = R.fp; end
[~, b, ext] = fileparts(fn);  desc = [b ext];
end

% ---------------------------------------------------------------------------
function v = ternary(c, a, b)
% TERNARY  a if c else b.
% INPUTS: c [logical]; a, b [any].  OUTPUTS: v [any].  REFERENCES: none.
if c, v = a; else, v = b; end
end

% ---------------------------------------------------------------------------
function h = git_hash(here)
% Short git hash for provenance (empty on failure -- never blocks a solve).
[rc, s] = system(sprintf('cd "%s" && git rev-parse --short HEAD 2>/dev/null', here));
if rc==0, h = strtrim(s); else, h = ''; end
end
