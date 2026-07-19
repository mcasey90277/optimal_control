% test_sosc_verdict_logic  Map (K,AS,IN) -> verdict per DESIGN sec 12.2 (7 ordered
% rules). Driven by the DIRECT reduced-Hessian inertia IN.red plus IN.robust and
% IN.sensStable; AS.nWeak gates only FAIL-vs-INCONCLUSIVE; AS.licq/redConsistent
% are retired (never referenced).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
mk = @(pass,signOK) struct('signOK',signOK,'pass',pass,'stat',0,'primalEq',0, ...
    'primalIneq',0,'dualFeas',0,'comp',0,'sign',1);
okK = mk(true,true);
% helper builders for the reduced-inertia struct + active set (new fields:
% red, robust, sensStable, method; NO redConsistent).
mkIN = @(rp,rn,rz) struct('red',struct('npos',rp,'nneg',rn,'nzero',rz), ...
    'robust',true,'sensStable',true,'method','reduced-eig','rankA',1, ...
    'nnegBand',[rn rn rn rn],'redMinEig',0);
mkAS = @(nWeak,licq) struct('licq',licq,'nWeak',nWeak,'m_active',1, ...
    'weakLabels',{{}},'nEq',1,'nStrong',0);

% (1) ERROR: KKT not pass
v = sosc_decide(mk(false,true), mkAS(0,true), mkIN(1,0,0));
assert(strcmp(v.verdict,'ERROR') && ...
    strcmp(v.status,'certified-feasibility+sosc-inconclusive'), 'ERROR pass=false');
% (1) ERROR: sign not resolved
v = sosc_decide(mk(true,false), mkAS(0,true), mkIN(1,0,0));
assert(strcmp(v.verdict,'ERROR'), 'ERROR signOK=false');

% (2) ~robust -> INCONCLUSIVE (scale)
INbig = mkIN(1,0,0); INbig.robust = false; INbig.method = 'scale-skip';
v = sosc_decide(okK, mkAS(0,true), INbig);
assert(strcmp(v.verdict,'INCONCLUSIVE') && ...
    strcmp(v.status,'certified-feasibility+sosc-inconclusive'), 'robust=false -> INCONCLUSIVE');

% (3) ~sensStable -> INCONCLUSIVE (zt-sensitive near-flat directions)
INsens = mkIN(1,0,0); INsens.sensStable = false; INsens.nnegBand = [2 2 1 0];
v = sosc_decide(okK, mkAS(0,true), INsens);
assert(strcmp(v.verdict,'INCONCLUSIVE') && ...
    strcmp(v.status,'certified-feasibility+sosc-inconclusive'), 'sensStable=false -> INCONCLUSIVE');

% (4) FAIL: red.nneg>0 AND nWeak==0 (and stable) -> feasible-only
v = sosc_decide(okK, mkAS(0,true), mkIN(1,1,0));
assert(strcmp(v.verdict,'FAIL') && strcmp(v.status,'feasible-only'), 'FAIL');

% (5) INCONCLUSIVE: red.nneg>0 but nWeak>0 -> cone strict subset
v = sosc_decide(okK, mkAS(1,true), mkIN(1,1,0));
assert(strcmp(v.verdict,'INCONCLUSIVE') && ...
    strcmp(v.status,'certified-feasibility+sosc-inconclusive'), 'INCONCLUSIVE weak');

% (6) PASS: red.nneg=0, red.nzero=0 -> certified-sosc
v = sosc_decide(okK, mkAS(0,true), mkIN(1,0,0));
assert(strcmp(v.verdict,'PASS') && strcmp(v.status,'certified-sosc'), 'PASS');

% (7) WEAK_MIN: red.nneg=0, red.nzero>0 -> certified-weak-min
v = sosc_decide(okK, mkAS(0,true), mkIN(1,0,3));
assert(strcmp(v.verdict,'WEAK_MIN') && strcmp(v.status,'certified-weak-min'), 'WEAK_MIN');
% WEAK_MIN robust to weakly-active junctions (never gated by them)
v = sosc_decide(okK, mkAS(2,false), mkIN(1,0,3));
assert(strcmp(v.verdict,'WEAK_MIN'), 'WEAK_MIN robust to weak/LICQ');

fprintf('test_sosc_verdict_logic PASSED\n');
