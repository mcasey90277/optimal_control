% TEST_FREETAUF_PORT  Contract test for the cScale (free-tau_f) branch of
% casadi_minfuel_sundman.m (ported 2026-09-07 after review finding E1).
%
% Three tiny (maxIter=5) warm re-solves of the certified flagship row:
%   (1) fixed mode, no opts               -- the byte-identity reference
%   (2) fixed mode, opts.freeTauf=false   -- must be IDENTICAL to (1): the
%                                            option parsing alone must not
%                                            perturb the solve
%   (3) free mode, opts.freeTauf=true + returnModel -- 8-row X/lamDef contract
%       kept, slack exposed (cScale, X9 [9 rows], lamDef9, taufSeed, tauf =
%       cScale*taufSeed), creg 'defect' block is 9*N rows, manifest
%       'tulip_free', and the generic gate (foc_check with that manifest)
%       runs mechanically on the 9-state model.
%   (4) a 9-row seed handed to the FIXED engine is sliced to 8 rows, not
%       rejected.
%   (5) opts.lamG0 (dual warm start) is accepted when it matches opti.g and
%       ignored with a warning when it does not.
% Convergence is NOT the subject (maxIter=5): the sit-still / reproduce-the-
% probe checks live in certify/gate_free_tauf.m and the 2026-09-07 record.
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths();
here = fileparts(mfilename('fullpath'));  cd(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
vcDir = fullfile(here, '..', '..', '..', 'verify_common');
addpath(vcDir);  setup_verify_common();

cacheFile = fullfile(here, '..', 'lib', 'sundman_minfuel_certified.mat');
if ~isfile(cacheFile)
    fprintf('SKIPPED -- cache absent\n'); return;
end
S = load(cacheFile);
p = cr3bp_lt_params(0.025, 15, 2100);
pSund = 1.5;   if isfield(S,'pSund') && ~isempty(S.pSund), pSund = S.pSund; end
tf = S.out.X(8,end);  N = numel(S.sigma) - 1;

o1 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    S.out.X, S.out.U, S.tauf0, pSund, 5, 0, true);
o2 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    S.out.X, S.out.U, S.tauf0, pSund, 5, 0, true, struct('freeTauf', false));
assert(isequal(o1.X, o2.X) && isequal(o1.U, o2.U) && isequal(o1.lamAll, o2.lamAll), ...
    'TEST_FREETAUF_PORT: freeTauf=false is not byte-identical to the no-opts call');
assert(~isfield(o1, 'cScale') && ~isfield(o2, 'cScale'), ...
    'TEST_FREETAUF_PORT: fixed mode must not expose cScale');

o3 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    S.out.X, S.out.U, S.tauf0, pSund, 5, 0, true, struct('freeTauf', true, 'returnModel', true));
assert(isequal(size(o3.X), [8 N+1]),  'TEST_FREETAUF_PORT: out.X must keep the 8-row contract');
assert(isequal(size(o3.X9), [9 N+1]), 'TEST_FREETAUF_PORT: out.X9 must be 9 rows');
assert(isequal(size(o3.lamDef), [8 N]) && isequal(size(o3.lamDef9), [9 N]), ...
    'TEST_FREETAUF_PORT: lamDef / lamDef9 shapes');
assert(isscalar(o3.cScale) && abs(o3.taufSeed - S.tauf0) == 0 && ...
       abs(o3.tauf - o3.cScale*S.tauf0) < 1e-12, 'TEST_FREETAUF_PORT: tauf must equal cScale*taufSeed');
assert(all(abs(o3.X9(9,:) - o3.cScale) < 1e-9), 'TEST_FREETAUF_PORT: cScale row must be constant');
assert(strcmp(o3.model.manifest, 'tulip_free'), 'TEST_FREETAUF_PORT: model.manifest');
labels = {o3.model.creg.label};
for lab = {'defect', 'betaNorm', 'thrLo', 'thrHi'}
    assert(any(strcmp(labels, lab{1})), 'TEST_FREETAUF_PORT: creg missing ''%s''', lab{1});
end
defRows = o3.model.creg(strcmp(labels, 'defect')).rows;
assert(numel(defRows) == 9*N, 'TEST_FREETAUF_PORT: defect block has %d rows, expected 9*N = %d', numel(defRows), 9*N);
assert(isfield(o3.boundSat, 'worst') && o3.cScale > 0.2 && o3.cScale < 5, 'TEST_FREETAUF_PORT: cScale outside its box');

o9 = o3;  o9.X = o3.X9;
rep = foc_check(o9, S.sigma, foc_manifest('tulip_free'), struct('eps', 0));
assert(isfinite(rep.kktStatInf) && isequal(size(rep.lam), [9 N+1]), ...
    'TEST_FREETAUF_PORT: foc_check did not run on the 9-state model');

o4 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    o3.X9, o3.U, S.tauf0, pSund, 5, 0, true);
assert(isequal(size(o4.X), [8 N+1]) && ~isfield(o4, 'cScale'), ...
    'TEST_FREETAUF_PORT: a 9-row seed must be sliced by the fixed engine');

% (5) dual warm start: the free solve's own multipliers are accepted (same NLP)
%     and a wrong-length vector is ignored with a warning, not an error
o5 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    o3.X9, o3.U, S.tauf0, pSund, 5, 0, true, struct('freeTauf', true, 'lamG0', o3.lamAll));
assert(isequal(size(o5.X), [8 N+1]) && numel(o5.lamAll) == numel(o3.lamAll), ...
    'TEST_FREETAUF_PORT: lamG0 re-solve shape');
% (the engine's later boundSaturation warning overwrites lastwarn, so capture
% the printed stream instead)
txt = evalc(['o6 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ' ...
             'S.out.X, S.out.U, S.tauf0, pSund, 5, 0, true, struct(''lamG0'', o3.lamAll));']);   % 9-state duals into the 8-state NLP
assert(contains(txt, 'opts.lamG0 has') && isequal(size(o6.X), [8 N+1]), ...
    'TEST_FREETAUF_PORT: mismatched lamG0 must warn and be ignored');

fprintf(['TEST_FREETAUF_PORT: PASS (fixed mode byte-identical; free mode keeps the 8-row ' ...
         'contract, exposes cScale=%.6f, 9*N defect rows, manifest tulip_free, foc_check runs)\n'], o3.cScale);
