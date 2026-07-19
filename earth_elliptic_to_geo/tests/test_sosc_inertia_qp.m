% Hand-built KKT inertia cases (no NLP), incl. reduced-Hessian inertia (sec 11.4):
%   PD reduced Hessian -> KKT (n, m_a, 0), red=(1,0,0), subspaceOK=true
%   Indefinite reduced Hessian -> red=(0,1,0) (one negative curvature dir)
%   Rank-deficient A -> KKT nzero>0, but red=(2,0,0) redConsistent (LICQ
%       deficiency is separated out; the reduced Hessian is actually PD)
%   WEAK_MIN -> red=(0,0,1): PSD reduced Hessian with one flat direction
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
tol = sosc_defaults();
% Case PD: H=I2, A=[1 1] -> reduced Hessian 1x1 = 1 > 0
IN = sosc_inertia(sparse(eye(2)), sparse([1 1]), tol);
assert(isequal([IN.npos IN.nneg IN.nzero],[2 1 0]) && IN.subspaceOK, 'PD case');
assert(isequal([IN.red.npos IN.red.nneg IN.red.nzero],[1 0 0]) && IN.redConsistent, ...
    'PD case reduced inertia (1,0,0)');
% Case indefinite: H=diag(1,-3), A=[1 1] -> reduced Hessian = (1-3)/2 = -1 < 0
IN2 = sosc_inertia(sparse(diag([1 -3])), sparse([1 1]), tol);
assert(~IN2.subspaceOK && IN2.nneg==2, 'indefinite case -> FAIL signature');
assert(isequal([IN2.red.npos IN2.red.nneg IN2.red.nzero],[0 1 0]) && IN2.redConsistent, ...
    'indefinite case reduced inertia (0,1,0)');
% Case rank-deficient A: two identical rows
IN3 = sosc_inertia(sparse(eye(3)), sparse([1 0 0; 1 0 0]), tol);
assert(IN3.nzero > 0, 'rank-deficient A -> nonzero nullity');
assert(isequal([IN3.red.npos IN3.red.nneg IN3.red.nzero],[2 0 0]) && IN3.redConsistent, ...
    'rank-deficient A reduced inertia (2,0,0), LICQ deficiency separated out');
% Case WEAK_MIN: H=diag(2,0), A=[1 0] -> r=1, m_a=1, n=2; KKT inertia (1,1,1);
% reduced Hessian PSD with one flat direction -> red=(0,0,1) (bang-bang signature)
IN4 = sosc_inertia(sparse(diag([2 0])), sparse([1 0]), tol);
assert(isequal([IN4.npos IN4.nneg IN4.nzero],[1 1 1]), 'WEAK_MIN KKT inertia (1,1,1)');
assert(IN4.red.nneg==0 && IN4.red.nzero==1 && IN4.redConsistent, ...
    'WEAK_MIN reduced inertia red.nneg==0, red.nzero==1, redConsistent');
% Direct 2x2-block coverage: D with a [[0 1];[1 0]] block (eigs -1,+1) + a 1x1 (+3).
% Bypasses ldl pivot selection so the 2x2 classification arithmetic is exercised.
Dtest = blkdiag([0 1; 1 0], 3);
[np, nn, nz] = count_inertia(Dtest, 1e-9);
assert(isequal([np nn nz],[2 1 0]), ...
    sprintf('2x2-block count wrong: got [%d %d %d], expected [2 1 0]', np, nn, nz));
fprintf('test_sosc_inertia_qp PASSED\n');
