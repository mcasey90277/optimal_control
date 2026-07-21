% Hand-built reduced-Hessian inertia cases (no NLP) for the FINAL direct method
% (DESIGN sec 12.1): Z=null(full(A)), RH=Z'HZ, ev=eig(RH). Each case is
% hand-verified; sensStable must be true for all four (single-signed spectra).
%   PD           H=I2,        A=[1 1]      -> Z=[1;-1]/sqrt2, RH=1  -> red=(1,0,0)
%   indefinite   H=diag(1,-3),A=[1 1]      -> RH=-1                 -> red.nneg=1
%   WEAK_MIN     H=diag(2,0), A=[1 0]      -> Z=[0;1],       RH=0  -> red=(0,0,1)
%   rank-defic.  H=I3,        A=[1 0 0;1 0 0] -> Z spans e2,e3, RH=I2 -> red=(2,0,0)
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
tol = sosc_defaults();

% Case PD: reduced Hessian 1x1 = 1 > 0 -> PASS-shape
IN = sosc_inertia(sparse(eye(2)), sparse([1 1]), tol);
assert(isequal([IN.red.npos IN.red.nneg IN.red.nzero],[1 0 0]), 'PD red=(1,0,0)');
assert(IN.sensStable==true, 'PD sensStable');
assert(strcmp(IN.method,'reduced-eig') && IN.robust, 'PD method=reduced-eig robust');
assert(IN.rankA==1, 'PD rankA=1');

% Case indefinite: reduced Hessian = -1 < 0 -> FAIL-shape, stable
IN2 = sosc_inertia(sparse(diag([1 -3])), sparse([1 1]), tol);
assert(isequal([IN2.red.npos IN2.red.nneg IN2.red.nzero],[0 1 0]), 'indef red=(0,1,0)');
assert(IN2.sensStable==true, 'indef sensStable');

% Case WEAK_MIN: reduced Hessian = 0 -> one flat direction
IN3 = sosc_inertia(sparse(diag([2 0])), sparse([1 0]), tol);
assert(IN3.red.nneg==0 && IN3.red.nzero==1, 'WEAK_MIN red.nneg==0 red.nzero==1');
assert(IN3.sensStable==true, 'WEAK_MIN sensStable');

% Case rank-deficient A (redundant constraint, handled by null()): red=(2,0,0)
IN4 = sosc_inertia(sparse(eye(3)), sparse([1 0 0; 1 0 0]), tol);
assert(isequal([IN4.red.npos IN4.red.nneg IN4.red.nzero],[2 0 0]), 'rank-def red=(2,0,0)');
assert(IN4.sensStable==true, 'rank-def sensStable');
assert(IN4.rankA==1, 'rank-def rankA=1 (redundant row collapsed by null())');
fprintf('test_sosc_inertia_qp PASSED\n');
