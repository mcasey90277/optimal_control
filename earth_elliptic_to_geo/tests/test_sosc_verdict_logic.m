root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
mk = @(stat,pass) struct('signOK',stat,'pass',pass,'stat',0,'primalEq',0, ...
    'primalIneq',0,'dualFeas',0,'comp',0,'sign',1);
okK = mk(true,true);
% PASS: kkt ok, licq ok, no weak, inertia subspaceOK
AS1 = struct('licq',true,'nWeak',0,'m_active',1,'weakLabels',{{}},'nEq',1,'nStrong',0);
IN1 = struct('subspaceOK',true,'nzero',0,'npos',2,'nneg',1,'expected',[2 1 0]);
assert(strcmp(sosc_decide(okK,AS1,IN1).verdict,'PASS'));
% FAIL: inertia wrong, no weak, licq ok
IN2 = IN1; IN2.subspaceOK=false; IN2.nzero=0;
assert(strcmp(sosc_decide(okK,AS1,IN2).verdict,'FAIL'));
% INCONCLUSIVE: weak present
AS3 = AS1; AS3.nWeak=1; AS3.weakLabels={'thrHi, node 204'};
assert(strcmp(sosc_decide(okK,AS3,IN1).verdict,'INCONCLUSIVE'));
% ERROR: kkt not pass
badK = mk(false,false);
assert(strcmp(sosc_decide(badK,AS1,IN1).verdict,'ERROR'));
fprintf('test_sosc_verdict_logic PASSED\n');
