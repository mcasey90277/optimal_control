% TEST_FOC_DUAL_TO_COSTATE  Unit test for the generic step-weighted map.
% Mirrors earth test_verify_pmp_mee Test 1, but at nx=9 to prove genericity.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root);
setup_verify_common;
tol = 1e-12;
sigma1 = [0; 0.3; 1.0];                       % non-uniform: h1=0.3, h2=0.7
Lam1a  = (2:10).';  Lam1b = (5:13).';         % nx=9, distinct rows
lam1 = foc_dual_to_costate([Lam1a, Lam1b], sigma1);
assert(isequal(size(lam1), [9 3]), 'wrong size');
assert(max(abs(lam1(:,1) - Lam1a)) < tol, 'node 1 one-sided');
assert(max(abs(lam1(:,3) - Lam1b)) < tol, 'node 3 one-sided');
expectMid = 0.3*Lam1a + 0.7*Lam1b;            % h-weighted (h1+h2=1)
assert(max(abs(lam1(:,2) - expectMid)) < tol, 'interior step-weighted');
assert(max(abs(lam1(:,2) - 0.5*(Lam1a+Lam1b))) > 1e-3, 'must differ from plain avg');
fprintf('test_foc_dual_to_costate: ALL PASS\n');
