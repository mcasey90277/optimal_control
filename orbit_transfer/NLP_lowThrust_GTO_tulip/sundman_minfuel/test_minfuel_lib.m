function test_minfuel_lib()
% TEST_MINFUEL_LIB  Cheap non-solve guardrail checks for the min-fuel library.
%
% Runs in seconds with NO CasADi/IPOPT solves -- safe to run before/after any
% refactor. Checks config consistency against the certified artifacts,
% filename encode/parse round-trips, stored-solution data integrity (lamDef
% shape vs states), and schedule sanity. Solver-level validation (reproducing
% a banked solve) is Phase 2 of CODE_CLEANUP_PLAN.md, not this file.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per check; errors on first FAIL)

here = fileparts(mfilename('fullpath'));  addpath(here);
cfg  = minfuel_config();
nOK  = 0;

% 1. tfMin matches the certified anchor (certified solved at 1.15x)
C = load(fullfile(here,'sundman_minfuel_certified.mat'));
tfAnchor = C.out.X(8,end);
assert(abs(tfAnchor/1.15 - cfg.tfMin) < 1e-6, ...
    'tfMin (%.10f) disagrees with certified anchor/1.15 (%.10f)', cfg.tfMin, tfAnchor/1.15);
nOK = nOK + 1;  fprintf('PASS %d: cfg.tfMin consistent with certified anchor\n', nOK);

% 2. schedules: strictly decreasing, in [0,1], sharpen ends at exactly 0
for sc = {cfg.schedSharpen, cfg.schedNeighbor}
    s = sc{1};
    assert(all(diff(s) < 0), 'schedule not strictly decreasing');
    assert(all(s >= 0 & s <= 1), 'schedule outside [0,1]');
end
assert(cfg.schedSharpen(end) == 0, 'schedSharpen must end at exactly eps=0');
assert(cfg.schedNeighbor(end) == 0, 'schedNeighbor must end at exactly eps=0');
nOK = nOK + 1;  fprintf('PASS %d: homotopy schedules sane (both end at eps=0)\n', nOK);

% 3. filename encode/parse round-trip at 0.001 granularity (incl. old %.2f
%    collision cases like 1.125 vs 1.13)
for f = [1.0 1.01 1.12 1.125 1.13 1.20 1.45 1.85 2.0]
    nm = cfg.fname('minfuel', f);
    fb = cfg.fparse(nm);
    assert(abs(fb - f) < 5e-4, 'filename round-trip failed: %.4f -> %s -> %.4f', f, nm, fb);
end
assert(~strcmp(cfg.fname('x',1.125), cfg.fname('x',1.13)), '0.001-granularity collision');
nOK = nOK + 1;  fprintf('PASS %d: filename encode/parse round-trips\n', nOK);

% 4. artifact integrity. The certified .mat predates the costate-return
%    feature (its `out` has no lamDef; the PMP numbers came from a later
%    re-solve), so basics are checked there and costate integrity on a
%    lamDef-bearing artifact (legacy down-sweep solution).
N = size(C.out.X,2) - 1;
assert(C.out.maxDefect < 1e-12, 'certified defect %.2g not machine-tight', C.out.maxDefect);
assert(abs(C.out.mf - (1 - 2.2640/15)) < 5e-4, 'certified final mass drifted');
L = load(fullfile(cfg.dirs.minfuel, 'legacy_ms_f1120.mat'));
NL = size(L.out.X,2) - 1;
assert(isequal(size(L.out.lamDef), [8 NL]), 'lamDef shape %s != 8x%d', mat2str(size(L.out.lamDef)), NL);
relTrans = abs(L.out.lamDef(7,end)) / max(abs(L.out.lamDef(7,:)));
assert(relTrans < 1e-3, 'legacy_ms_f1120 relative transversality %.2g >= 1e-3', relTrans);
assert(L.out.primerAlignDeg < 0.2, 'legacy_ms_f1120 primer %.3f deg >= 0.2', L.out.primerAlignDeg);
nOK = nOK + 1;  fprintf('PASS %d: artifact integrity (certified defect %.1e; lamDef 8x%d, relTrans %.1e)\n', ...
                        nOK, C.out.maxDefect, NL, relTrans);

% 5. physics constants: params round-trip and ND time scale
p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
assert(abs(p.tStar - 382981.289129055) < 1e-6, 'tStar drifted');
assert(abs(cfg.tfMin*p.tStar/86400 - 27.8845) < 1e-3, 'tfMin != 27.8845 d');
nOK = nOK + 1;  fprintf('PASS %d: physics constants (tfMin = %.4f d)\n', nOK, cfg.tfMin*p.tStar/86400);

% 6. new-layout result files (if any) carry provenance meta
if exist(cfg.dirs.minfuel,'dir')
    dd = dir(fullfile(cfg.dirs.minfuel,'minfuel_f*.mat'));
    for k = 1:numel(dd)
        R = load(fullfile(dd(k).folder, dd(k).name));
        assert(isfield(R,'out') && isfield(R.out,'meta'), '%s lacks meta', dd(k).name);
        assert(all(isfield(R.out.meta, {'date','githash','seed','sched'})), ...
               '%s meta incomplete', dd(k).name);
    end
    fprintf('     (checked meta on %d new-layout result files)\n', numel(dd));
end
nOK = nOK + 1;  fprintf('PASS %d: result-file provenance\n', nOK);

fprintf('\nALL %d CHECKS PASSED\n', nOK);
end
