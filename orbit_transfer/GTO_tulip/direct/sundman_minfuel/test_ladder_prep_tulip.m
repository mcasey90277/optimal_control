% TEST_LADDER_PREP_TULIP  boundSat fields + opts back-compat + chain helper.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','..','..','cr3bp_common'));  setup_cr3bp_common();
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
assert(isfile('sundman_minfuel_certified.mat'), 'need the certified cache');
C = load('sundman_minfuel_certified.mat');
% (a) chain helper from the certified solution to a 20 mN config
p20 = cr3bp_lt_params(0.020, cfg.m0kg, cfg.ispS);
tfNew = C.out.X(8,end);                     % same tf, new thrust (pilot pattern)
[sg, X0, U0, tauf0, fp] = chain_rung_seed_tulip(C, tfNew, p20, struct('note','test'));
assert(size(X0,1)==8 && size(U0,1)==4 && numel(sg)==size(X0,2), 'shapes');
assert(abs(X0(8,end)-tfNew) < 1e-9 && abs(X0(8,1)) < 1e-12, 'time row endpoints');
assert(tauf0 > 0 && issorted(sg) && abs(sg(1))<1e-12 && abs(sg(end)-1)<1e-12, 'sigma/tauf0');
assert(fp.thrustN==0.020 && isfield(fp,'chainedFrom'), 'fp thrust + provenance');
assert(max(abs(X0(1:6,1).'  - C.rv0(:).')) < 1e-9, 'rv0 pinned');
assert(max(abs(X0(1:6,end).'- C.rvf(:).')) < 1e-9, 'rvf pinned');
% (b) same-thrust chain refused
ok = false;
try, chain_rung_seed_tulip(C, tfNew, p, struct());
catch err, ok = strcmp(err.identifier,'chain_rung_seed_tulip:sameThrust'); end
assert(ok, 'same-thrust chain must error');
% (c) solver back-compat + boundSat: 3-iter probe, 14-arg call == 15-arg default
o14 = casadi_minfuel_sundman(C.sigma, tfNew, C.rv0, C.rvf, p.Tmax, p.c, p.muStar, ...
        C.out.X, C.out.U, C.tauf0, cfg.pSund, 3, 0, true);
o15 = casadi_minfuel_sundman(C.sigma, tfNew, C.rv0, C.rvf, p.Tmax, p.c, p.muStar, ...
        C.out.X, C.out.U, C.tauf0, cfg.pSund, 3, 0, true, struct());
assert(isfield(o14,'boundSat') && isfield(o14.boundSat,'minSlack') && ...
       isfield(o14.boundSat,'worst') && islogical(o14.boundSat.hit), 'boundSat fields');
assert(max(abs(o14.X(:)-o15.X(:))) < 1e-12, 'empty opts must be byte-compatible');
fprintf('test_ladder_prep_tulip: ALL PASS\n');
