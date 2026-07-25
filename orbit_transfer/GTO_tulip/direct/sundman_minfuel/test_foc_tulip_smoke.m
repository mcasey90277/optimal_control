% TEST_FOC_TULIP_SMOKE  Byte-identity regression for the Task 8 returnModel
% hook on casadi_minfuel_sundman.m: two tiny (maxIter=5) warm re-solves of
% the certified flagship row -- one WITHOUT opts.returnModel, one WITH --
% must produce IDENTICAL X/U (the registry bookkeeping only reads
% size(opti.g,1); it must never perturb the solve). Also checks the expected
% creg labels ('defect','betaNorm','thrLo','thrHi') are present when the flag
% is on.
here = fileparts(mfilename('fullpath'));  cd(here);
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));

S = load(fullfile(here, 'sundman_minfuel_certified.mat'));
p = cr3bp_lt_params(0.025, 15, 2100);
pSund = 1.5;   if isfield(S,'pSund') && ~isempty(S.pSund), pSund = S.pSund; end
tf = S.out.X(8,end);

o1 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    S.out.X, S.out.U, S.tauf0, pSund, 5, 0, true);                        % no opts arg at all
o2 = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    S.out.X, S.out.U, S.tauf0, pSund, 5, 0, true, struct('returnModel', true));

assert(isequal(o1.X, o2.X), 'TEST_FOC_TULIP_SMOKE: X differs with returnModel on -- not byte-identical');
assert(isequal(o1.U, o2.U), 'TEST_FOC_TULIP_SMOKE: U differs with returnModel on -- not byte-identical');
assert(~isfield(o1, 'model'), 'TEST_FOC_TULIP_SMOKE: out.model present without returnModel');
assert(isfield(o2, 'model') && isfield(o2.model, 'creg'), ...
    'TEST_FOC_TULIP_SMOKE: out.model.creg missing with returnModel=true');

labels = {o2.model.creg.label};
expected = {'defect', 'betaNorm', 'thrLo', 'thrHi'};
for k = 1:numel(expected)
    assert(any(strcmp(labels, expected{k})), ...
        'TEST_FOC_TULIP_SMOKE: creg missing label ''%s''', expected{k});
end

fprintf('TEST_FOC_TULIP_SMOKE: PASS (X/U byte-identical with/without returnModel; creg labels %s present)\n', ...
    strjoin(expected, ','));
