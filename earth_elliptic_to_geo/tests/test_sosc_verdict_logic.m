% test_sosc_verdict_logic  Map (K,AS,IN) -> verdict per DESIGN sec 11.5 (ordered).
% Verdicts driven by the REDUCED-Hessian inertia IN.red (+ IN.redConsistent);
% AS.nWeak gates only FAIL-vs-INCONCLUSIVE; AS.licq is reported, never gates.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
mk = @(pass,signOK) struct('signOK',signOK,'pass',pass,'stat',0,'primalEq',0, ...
    'primalIneq',0,'dualFeas',0,'comp',0,'sign',1);
okK = mk(true,true);
% helper builders for the reduced-inertia struct + active set
mkIN = @(rp,rn,rz,cons) struct('red',struct('npos',rp,'nneg',rn,'nzero',rz), ...
    'redConsistent',cons,'rankA',1,'subspaceOK',rn==0&&rz==0, ...
    'npos',rp+1,'nneg',rn+1,'nzero',rz,'expected',[2 1 0]);
mkAS = @(nWeak,licq) struct('licq',licq,'nWeak',nWeak,'m_active',1, ...
    'weakLabels',{{}},'nEq',1,'nStrong',0);

% PASS: red=(1,0,0), consistent, K ok  -> certified-sosc
v = sosc_decide(okK, mkAS(0,true), mkIN(1,0,0,true));
assert(strcmp(v.verdict,'PASS') && strcmp(v.status,'certified-sosc'), 'PASS');

% WEAK_MIN: red.nneg=0, red.nzero>0, consistent -> certified-weak-min
v = sosc_decide(okK, mkAS(0,true), mkIN(1,0,3,true));
assert(strcmp(v.verdict,'WEAK_MIN') && strcmp(v.status,'certified-weak-min'), 'WEAK_MIN');
% WEAK_MIN robust to weakly-active junctions + LICQ deficiency (never gated by them)
v = sosc_decide(okK, mkAS(2,false), mkIN(1,0,3,true));
assert(strcmp(v.verdict,'WEAK_MIN'), 'WEAK_MIN robust to weak/LICQ');

% FAIL: red.nneg>0 AND nWeak==0 -> feasible-only (proven saddle)
v = sosc_decide(okK, mkAS(0,true), mkIN(1,1,0,true));
assert(strcmp(v.verdict,'FAIL') && strcmp(v.status,'feasible-only'), 'FAIL');

% INCONCLUSIVE (red.nneg>0 but nWeak>0): critical cone strict subset
v = sosc_decide(okK, mkAS(1,true), mkIN(1,1,0,true));
assert(strcmp(v.verdict,'INCONCLUSIVE') && ...
    strcmp(v.status,'certified-feasibility+sosc-inconclusive'), 'INCONCLUSIVE weak');

% INCONCLUSIVE (~redConsistent): rank estimate untrustworthy
v = sosc_decide(okK, mkAS(0,true), mkIN(1,0,0,false));
assert(strcmp(v.verdict,'INCONCLUSIVE'), 'INCONCLUSIVE redConsistent=false');

% ERROR: KKT not pass (also feasibility+sosc-inconclusive status, non-demoting)
v = sosc_decide(mk(false,true), mkAS(0,true), mkIN(1,0,0,true));
assert(strcmp(v.verdict,'ERROR') && ...
    strcmp(v.status,'certified-feasibility+sosc-inconclusive'), 'ERROR pass=false');
% ERROR: sign not resolved
v = sosc_decide(mk(true,false), mkAS(0,true), mkIN(1,0,0,true));
assert(strcmp(v.verdict,'ERROR'), 'ERROR signOK=false');

fprintf('test_sosc_verdict_logic PASSED\n');
