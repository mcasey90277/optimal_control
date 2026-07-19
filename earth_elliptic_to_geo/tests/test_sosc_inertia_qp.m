% Hand-built KKT inertia cases (no NLP):
%   PD reduced Hessian -> inertia (n, m_a, 0), subspaceOK=true
%   Indefinite reduced Hessian -> wrong inertia, subspaceOK=false
%   Rank-deficient A -> nzero>0
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
tol = sosc_defaults();
% Case PD: H=I2, A=[1 1] -> reduced Hessian 1x1 = 1 > 0
IN = sosc_inertia(sparse(eye(2)), sparse([1 1]), tol);
assert(isequal([IN.npos IN.nneg IN.nzero],[2 1 0]) && IN.subspaceOK, 'PD case');
% Case indefinite: H=diag(1,-3), A=[1 1] -> reduced Hessian = (1-3)/2 = -1 < 0
IN2 = sosc_inertia(sparse(diag([1 -3])), sparse([1 1]), tol);
assert(~IN2.subspaceOK && IN2.nneg==2, 'indefinite case -> FAIL signature');
% Case rank-deficient A: two identical rows
IN3 = sosc_inertia(sparse(eye(3)), sparse([1 0 0; 1 0 0]), tol);
assert(IN3.nzero > 0, 'rank-deficient A -> nonzero nullity');
% Direct 2x2-block coverage: D with a [[0 1];[1 0]] block (eigs -1,+1) + a 1x1 (+3).
% Bypasses ldl pivot selection so the 2x2 classification arithmetic is exercised.
Dtest = blkdiag([0 1; 1 0], 3);
[np, nn, nz] = count_inertia(Dtest, 1e-9);
assert(isequal([np nn nz],[2 1 0]), ...
    sprintf('2x2-block count wrong: got [%d %d %d], expected [2 1 0]', np, nn, nz));
fprintf('test_sosc_inertia_qp PASSED\n');
