% TEST_FOC_ELFO_SMOKE  Byte-identity regression for the Task 9 returnModel hook
% on BOTH ELFO solvers: casadi_energy_freetf.m and casadi_mintime_freetf.m.
%
% Same pattern as GTO_tulip/direct/sundman_minfuel/test_foc_tulip_smoke.m
% (Task 8): two tiny (maxIter=5) warm re-solves of a certified row -- one
% WITHOUT opts.returnModel, one WITH -- must produce IDENTICAL X/U (the
% registry bookkeeping only reads size(opti.g,1); it must never perturb the
% solve). Also checks the expected creg labels are present when the flag is
% on: 'defect','betaNorm','thrLo','thrHi' for the energy/fuel solver (all
% four throttle-inclusive groups), 'defect','betaNorm' ONLY for the mintime
% solver (hard all-burn s==1, no throttle rows).
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();

%% ---- (1) casadi_energy_freetf.m, at the certified energy seed -------------
energyCacheFile = fullfile(here, 'results', 'energy_elfo_freetf.mat');
if ~isfile(energyCacheFile)
    fprintf('SKIPPED -- cache absent\n'); return;
end
S = load(energyCacheFile);
oBase = struct('moonZone', S.moonZone, 'pSund', S.pSund, 'qSund', S.qSund, ...
    'muGain', S.muGain, 'epsilon', 1, 'tfTarget', S.X(8,end), 'maxIter', 5, ...
    'warmTight', true);

% Tmax/cEx/muStar: the campaign's own values (nominal 25 mN rung, matches
% this seed's provenance) -- from cr3bp_lt_params, not zeros.
p = cr3bp_lt_params(0.025, 15, 2100);
o1 = casadi_energy_freetf(S.sigma, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, S.X, S.U, S.tauf0, oBase);
o2opts = oBase;  o2opts.returnModel = true;
o2 = casadi_energy_freetf(S.sigma, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, S.X, S.U, S.tauf0, o2opts);

assert(isequal(o1.X, o2.X), 'TEST_FOC_ELFO_SMOKE: energy X differs with returnModel on -- not byte-identical');
assert(isequal(o1.U, o2.U), 'TEST_FOC_ELFO_SMOKE: energy U differs with returnModel on -- not byte-identical');
assert(~isfield(o1, 'model'), 'TEST_FOC_ELFO_SMOKE: energy out.model present without returnModel');
assert(isfield(o2, 'model') && isfield(o2.model, 'creg'), ...
    'TEST_FOC_ELFO_SMOKE: energy out.model.creg missing with returnModel=true');

labelsE = {o2.model.creg.label};
expectedE = {'defect', 'betaNorm', 'thrLo', 'thrHi'};
for k = 1:numel(expectedE)
    assert(any(strcmp(labelsE, expectedE{k})), ...
        'TEST_FOC_ELFO_SMOKE: energy creg missing label ''%s''', expectedE{k});
end
fprintf('TEST_FOC_ELFO_SMOKE [energy]: PASS (X/U byte-identical; creg labels %s present)\n', ...
    strjoin(expectedE, ','));

%% ---- (2) casadi_mintime_freetf.m, at the certified min-time anchor --------
mintimeCacheFile = fullfile(here, 'results', 'mintime_elfo.mat');
if ~isfile(mintimeCacheFile)
    fprintf('SKIPPED -- cache absent\n'); return;
end
M = load(mintimeCacheFile);
mBase = struct('moonZone', M.moonZone, 'pSund', M.pSund, 'qSund', M.qSund, ...
    'maxIter', 5, 'warmTight', true);

m1 = casadi_mintime_freetf(M.sigma, M.rv0, M.rvf, p.Tmax, p.c, p.muStar, M.X, M.U, M.tauf0, mBase);
m2opts = mBase;  m2opts.returnModel = true;
m2 = casadi_mintime_freetf(M.sigma, M.rv0, M.rvf, p.Tmax, p.c, p.muStar, M.X, M.U, M.tauf0, m2opts);

assert(isequal(m1.X, m2.X), 'TEST_FOC_ELFO_SMOKE: mintime X differs with returnModel on -- not byte-identical');
assert(isequal(m1.U, m2.U), 'TEST_FOC_ELFO_SMOKE: mintime U differs with returnModel on -- not byte-identical');
assert(~isfield(m1, 'model'), 'TEST_FOC_ELFO_SMOKE: mintime out.model present without returnModel');
assert(isfield(m2, 'model') && isfield(m2.model, 'creg'), ...
    'TEST_FOC_ELFO_SMOKE: mintime out.model.creg missing with returnModel=true');

labelsM = {m2.model.creg.label};
expectedM = {'defect', 'betaNorm'};
for k = 1:numel(expectedM)
    assert(any(strcmp(labelsM, expectedM{k})), ...
        'TEST_FOC_ELFO_SMOKE: mintime creg missing label ''%s''', expectedM{k});
end
assert(~any(strcmp(labelsM, 'thrLo')) && ~any(strcmp(labelsM, 'thrHi')), ...
    'TEST_FOC_ELFO_SMOKE: mintime creg unexpectedly carries a throttle-bound label (all-burn solver)');
fprintf('TEST_FOC_ELFO_SMOKE [mintime]: PASS (X/U byte-identical; creg labels %s present, no throttle rows)\n', ...
    strjoin(expectedM, ','));

fprintf('TEST_FOC_ELFO_SMOKE: ALL PASS\n');
