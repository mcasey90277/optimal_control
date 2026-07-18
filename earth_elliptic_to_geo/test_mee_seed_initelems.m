% test_mee_seed_initelems.m
here = fileparts(mfilename('fullpath')); cd(here);
par = kepler_lt_params(10,1500,2000);
% default (no initElems) MUST equal the legacy literal, byte-for-byte
[~, Xd] = mee_seed(par, struct('thr',0.4,'betaMode','tangential','N',20,'nRev',1));
legacy = [11625/par.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
assert(isequal(Xd(:,1), legacy), 'default initial node must be the legacy literal');
% custom initElems honored at node 1
ci = [0.30; 0.60; 0.0; 0.10; 0.0; 1; 0];
[~, Xc] = mee_seed(par, struct('thr',0.4,'betaMode','tangential','N',20,'nRev',1,'initElems',ci));
assert(isequal(Xc(:,1), ci), 'custom initElems must set node 1');
fprintf('test_mee_seed_initelems PASSED\n');
