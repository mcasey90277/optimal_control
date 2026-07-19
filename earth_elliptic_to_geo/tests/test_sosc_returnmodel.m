% TEST_SOSC_RETURNMODEL  (a) returnModel=false leaves the solve output
% unchanged (numerics invariant); (b) returnModel=true exposes opti + a
% registry whose row ranges partition opti.g exactly once.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

saved = sosc_load_row(fullfile(module_root(),'results','MEE_M2_10N.mat'));
par = kepler_lt_params(saved.thrustN,saved.m0kg,saved.ispS);
base = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',saved.maxIter,'warmTight',true,'printLevel',0);

oOff = casadi_lt_mee(saved.sigma,saved.X,saved.U,saved.dL, base);
optsOn = base; optsOn.returnModel = true;
oOn  = casadi_lt_mee(saved.sigma,saved.X,saved.U,saved.dL, optsOn);

% (a) numerics invariant
assert(isequal(oOff.X,oOn.X) && isequal(oOff.U,oOn.U) && oOff.dL==oOn.dL, ...
    'returnModel changed the numeric solution');
assert(~isfield(oOff,'model'), 'model leaked with flag off');

% (b) registry partitions opti.g
m = size(oOn.model.opti.g,1);
covered = [];
for i = 1:numel(oOn.model.creg), covered = [covered, oOn.model.creg(i).rows]; end %#ok<AGROW>
assert(isequal(sort(covered(:)'),1:m), 'creg rows must partition 1..m exactly once');
labels = {oOn.model.creg.label};
assert(any(strcmp(labels,'defect')) && any(strcmp(labels,'betaNorm')) && ...
       any(strcmp(labels,'termBC')), 'expected core labels present');
fprintf('test_sosc_returnmodel PASSED (m=%d constraint rows, %d groups)\n', m, numel(oOn.model.creg));
