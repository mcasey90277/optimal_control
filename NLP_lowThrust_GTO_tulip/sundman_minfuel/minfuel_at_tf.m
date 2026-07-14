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
% INPUTS:
%   factor  - t_f / t_f^min [scalar]
%   options (name-value):
%     'seed'       - 'energy' | 'neighbor' | file path      [default 'energy']
%     'seedFactor' - factor of the neighbor solution        [required for 'neighbor']
%     'sched'      - homotopy epsilon schedule override     [default by seed type]
%     'maxIter'    - IPOPT iteration cap                    [default cfg.maxIter]
%     'branch'     - branch tag recorded in meta + filename suffix, e.g.
%                    'up','dn','en'                          [default 'en'|'nb']
%     'outFile'    - output path override                   [default results/minfuel/]
%     'save'       - write the .mat                          [default true]
%
% OUTPUTS:
%   out - solver struct (X,U,lamDef,switches,edge,maxDefect,primerAlignDeg,...)
%         plus .factor .tf .tf_days .dV .prop_kg .certified (logical: at least
%         one schedule step converged tight) .meta (provenance: date, git hash,
%         seed source, schedule, per-step table, ipopt statuses). An uncertified
%         result is NOT saved -- a loose iterate must never become a neighbor seed.
%
% REFERENCES:
%   [1] LOW_THRUST_MINFUEL_CAMPAIGN.md ("Down-sweep CRACKED": backbone+sharpen).
%   [2] CODE_CLEANUP_PLAN.md (driver consolidation rationale).

here = fileparts(mfilename('fullpath'));  addpath(here);
cfg  = minfuel_config();
op   = parse_opts(varargin, cfg);
p    = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
tf   = factor * cfg.tfMin;

% --- resolve the seed -------------------------------------------------------
switch op.seedKind
    case 'energy'
        E = load(find_energy_file(here, cfg, factor));
        sigma=E.sigma; rv0=E.rv0; rvf=E.rvf; tauf0=E.tauf0;
        Xk=E.X; Uk=E.U; firstLoose=false; needClean=true;
        seedDesc = sprintf('energy backbone f=%.3f', factor);
    otherwise   % 'neighbor' or explicit file
        [S, seedDesc] = load_seed_file(here, cfg, op);
        sigma=S.sigma; rv0=S.rv0; rvf=S.rvf; tauf0=S.tauf0;
        Xk=S.X; Uk=S.U;
        Xk(8,:) = Xk(8,:) * (tf / Xk(8,end));   % rescale time state to new t_f
        firstLoose=true; needClean=false;       % first sched step is the move
end

fprintf('MINFUEL_AT_TF: factor=%.3f  t_f=%.4f ND (%.2f d)  seed=%s\n', ...
        factor, tf, tf*p.tStar/86400, seedDesc);

% --- (1) tight re-clean of an energy seed (same t_f -> no wedge) ------------
stat = {};
if needClean
    oT = casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar, ...
                                Xk,Uk,tauf0,cfg.pSund,op.maxIter,1,true);
    fprintf('  re-clean energy: ok=%d defect=%.2g\n', oT.success, oT.maxDefect);
    stat{end+1} = sprintf('reclean:%s', oT.ipoptStatus);
    if oT.success && oT.maxDefect < 1e-6, Xk=oT.X; Uk=oT.U; end
end

% --- (2) homotopy sharpen ---------------------------------------------------
best = [];  o = [];  tbl = zeros(numel(op.sched), 4);
for ke = 1:numel(op.sched)
    e = op.sched(ke);
    tight = ~(firstLoose && ke==1);
    o = casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar, ...
                               Xk,Uk,tauf0,cfg.pSund,op.maxIter,e,tight);
    ok = o.success && o.maxDefect < 1e-6;
    tbl(ke,:) = [e, o.maxDefect, o.switches, 100*o.edge];
    stat{end+1} = sprintf('eps=%.4g:%s', e, o.ipoptStatus); %#ok<AGROW>
    fprintf('  eps=%.4g: ok=%d defect=%.2g sw=%d edge=%.1f%%\n', ...
            e, ok, o.maxDefect, o.switches, 100*o.edge);
    if ok, Xk=o.X; Uk=o.U; best=o; end
end
certified = ~isempty(best);
if ~certified
    warning('minfuel_at_tf:noCleanStep', ...
        ['no schedule step converged tight at factor %.3f; returning the last ' ...
         'UNCERTIFIED attempt (will NOT be saved)'], factor);
    best = o;
end

% --- package with provenance ------------------------------------------------
out = best;
out.certified = certified;
out.factor  = factor;  out.tf = tf;  out.tf_days = tf*p.tStar/86400;
out.dV      = p.c*log(1/best.mf)*p.lStar/p.tStar;
out.prop_kg = p.m0kg*(1-best.mf);
out.meta = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
    'githash', git_hash(here), 'seed', seedDesc, 'branch', op.branch, ...
    'sched', op.sched, 'maxIter', op.maxIter, 'pSund', cfg.pSund, ...
    'tfMin', cfg.tfMin, 'stepTable', tbl, 'ipoptStatuses', {stat}, ...
    'solver', 'casadi_minfuel_sundman (CasADi+IPOPT, Sundman trapezoid)');

fprintf('MINFUEL_AT_TF done: f=%.3f dV=%.4f km/s sw=%d edge=%.1f%% defect=%.2g primer=%.3f\n', ...
        factor, out.dV, best.switches, 100*best.edge, best.maxDefect, best.primerAlignDeg);

if op.save && ~certified
    warning('minfuel_at_tf:skipSaveUncertified', ...
        ['factor %.3f did not converge tight; NOT writing an output file (a loose ' ...
         'iterate would poison neighbor-seed lookups). Inspect the returned struct instead.'], factor);
elseif op.save
    if isempty(op.outFile)
        if ~exist(cfg.dirs.minfuel,'dir'), mkdir(cfg.dirs.minfuel); end
        base = cfg.fname('minfuel', factor);
        op.outFile = fullfile(cfg.dirs.minfuel, strrep(base, '.mat', ['_' op.branch '.mat']));
    end
    save(op.outFile, 'out', 'sigma', 'tauf0', 'rv0', 'rvf', 'factor');
    fprintf('  WROTE %s\n', op.outFile);
end
end

% ---------------------------------------------------------------------------
function op = parse_opts(args, cfg)
% Name-value option parsing with seed-dependent defaults.
op = struct('seedKind','energy','seedFactor',NaN,'sched',[],'maxIter',cfg.maxIter, ...
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
        case 'branch',     op.branch     = args{k+1};
        case 'outfile',    op.outFile    = args{k+1};
        case 'save',       op.save       = args{k+1};
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
function [S, desc] = load_seed_file(here, cfg, op)
% Load a bang-bang seed: explicit file, or resolve from seedFactor (new
% minfuel results, then legacy ms_<f>.mat, then the certified anchor).
if strcmp(op.seedKind, 'file')
    fn = op.seedFile;
else
    cand = {fullfile(cfg.dirs.minfuel, strrep(cfg.fname('minfuel',op.seedFactor),'.mat','_en.mat')), ...
            fullfile(cfg.dirs.minfuel, strrep(cfg.fname('minfuel',op.seedFactor),'.mat','_nb.mat')), ...
            fullfile(cfg.dirs.minfuel, cfg.fname('legacy_ms', op.seedFactor)), ...
            fullfile(here, 'sundman_minfuel_certified.mat')};
    fn = '';
    for k = 1:numel(cand)
        if isfile(cand{k}), fn = cand{k}; break; end
    end
    if isempty(fn), error('minfuel_at_tf:noSeed','no seed found for factor %.3f', op.seedFactor); end
end
R = load(fn);
if isfield(R,'out'), S.X = R.out.X; S.U = R.out.U; else, S.X = R.X; S.U = R.U; end
S.sigma = R.sigma;  S.tauf0 = R.tauf0;  S.rv0 = R.rv0;  S.rvf = R.rvf;
[~, b, ext] = fileparts(fn);  desc = [b ext];
end

% ---------------------------------------------------------------------------
function h = git_hash(here)
% Short git hash for provenance (empty on failure -- never blocks a solve).
[rc, s] = system(sprintf('cd "%s" && git rev-parse --short HEAD 2>/dev/null', here));
if rc==0, h = strtrim(s); else, h = ''; end
end
