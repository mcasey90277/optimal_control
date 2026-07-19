% TEST_VERIFY_PMP_MEE  No-solve unit tests for the dual->costate map
% (mee_dual_to_costate) and the primer/switching-function extraction
% (mee_primer_switch), on hand-built synthetic data with an exactly known
% answer -- exercises the two building blocks of verify_pmp_mee.m without
% ever calling casadi_lt_mee (no NLP solve).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
tol = 1e-10;

%% Test 1 -- dual -> nodal costate, step-weighted adjacent-interval average,
% on a NON-uniform 2-interval synthetic mesh (catches an unweighted-average
% bug, which would only agree with the weighted formula on a uniform mesh).
sigma1 = [0; 0.3; 1.0];                 % h1=0.3, h2=0.7
Lam1a  = [2;3;4;5;6;7;8];               % interval-1 dual (distinct rows, catches transposition)
Lam1b  = [5;6;7;8;9;10;11];             % interval-2 dual
LamDef1 = [Lam1a, Lam1b];               % [7x2]
lam1 = mee_dual_to_costate(LamDef1, sigma1);
assert(isequal(size(lam1), [7 3]), 'Test1: wrong output size');
assert(max(abs(lam1(:,1) - Lam1a)) < tol, 'Test1: node 1 must equal interval-1 dual (one-sided)');
assert(max(abs(lam1(:,3) - Lam1b)) < tol, 'Test1: node 3 must equal interval-2 dual (one-sided)');
h1 = 0.3; h2 = 0.7;
lamMidExpected = (h1*Lam1a + h2*Lam1b) / (h1+h2);   % = 0.3*Lam1a + 0.7*Lam1b (h1+h2=1)
assert(max(abs(lam1(:,2) - lamMidExpected)) < tol, 'Test1: interior node must be step-weighted average');
% sanity: NOT the plain (unweighted) average, since h1~=h2
plainAvg = 0.5*(Lam1a + Lam1b);
assert(max(abs(lam1(:,2) - plainAvg)) > 1e-3, 'Test1: weighted and plain average must differ on this non-uniform mesh');
fprintf('Test 1 (dual->costate, step-weighted average) PASSED\n');

%% Common synthetic-state setup for Tests 2-4: P=1,ex=0,ey=0,m=1,t=0, and a
% node longitude chosen (via sigma,dL) to land exactly at L=0 (cL=1,sL=0),
% which makes the closed-form Gauss-matrix algebra hand-computable.
mu = 1;  Tm = 0.01;  cEx = 0.05;         % "low thrust": Tm << Ldot0 ~ O(1)
par = struct('mu', mu, 'Tmax', Tm, 'c', cEx);
sigmaNode = -1;  dLnode = pi;             % L = pi + sigma*dL = pi - pi = 0

%% Test 2 -- K_L=0 case (hx=hy=0, so hterm=0 identically): clean primer
% alignment (0 deg) and correct switching-sign in two lam_m regimes.
% B(X) at this state (hand-derived, lt_mee_rhs.m sec. algebra, sqPmu=1,
% Z=1, A1=2, A2=0, Xh=1, hterm=0):
%   B = [0,2,0; 0,2,0; -1,0,0; 0,0,0.5; 0,0,0]   (rows P,ex,ey,hx,hy)
X2 = [1;0;0;0;0;1;0];
lam2 = zeros(7,1);  lam2(1) = 1;          % lam_P=1 only -> p_el = B(1,:)' = [0;2;0]
beta2 = [0;-1;0];  thr2 = 1;
U2 = [beta2; thr2];

% burn regime: lam_m large enough to flip C1 negative -> S<0 predicted
lamB = lam2;  lamB(6) = 200;
[pv, S, info] = mee_primer_switch(X2, U2, lamB, sigmaNode, dLnode, par);
assert(abs(info.KL) < tol, 'Test2: K_L must be exactly 0 at hx=hy=0');
% expected primerVec = Tm*Ldot0*p_el = 0.01*1*[0;2;0] = [0;0.02;0] (K_L=0 kills the e3 term)
assert(max(abs(pv - [0;0.02;0])) < tol, 'Test2: primerVec (burn regime) mismatch');
predDir = -pv / norm(pv);
assert(max(abs(predDir - beta2)) < tol, 'Test2: predicted direction must exactly match beta2 (0 deg misalignment)');
Sexpected = (Tm*dot([0;2;0],beta2) - (Tm/cEx)*200 + 1) * 1;   % C1*Ldot0, lam_t=0
assert(abs(S - Sexpected) < tol, 'Test2: S (burn regime) mismatch');
assert(S < 0, 'Test2: burn regime must give S<0');

% coast regime: lam_m=0 -> C1~1>0 -> S>0 predicted
lamC = lam2;  lamC(6) = 0;
[~, Sc, ~] = mee_primer_switch(X2, U2, lamC, sigmaNode, dLnode, par);
Scexpected = (Tm*dot([0;2;0],beta2) + 1) * 1;
assert(abs(Sc - Scexpected) < tol, 'Test2: S (coast regime) mismatch');
assert(Sc > 0, 'Test2: coast regime must give S>0');
fprintf('Test 2 (K_L=0 primer/switching, hand-checked both signs) PASSED\n');

%% Test 3 -- K_L~=0 case (hy=0.1 -> hterm=-0.1 at L=0): exercises the
% Ldot-on-control coupling term in both primerVec and S.
%   B = [0,2,0; 0,2,0; -1,0,0; 0,0,0.505; 0,0,0],  K_L = Tm*(-0.1) = -0.001
X3 = [1;0;0;0;0.1;1;0];
lam3 = zeros(7,1);  lam3(1) = 1;  lam3(4) = 2;    % lam_P=1, lam_hx=2 -> p_el=[0;2;1.01]
beta3 = [0;0;1];  thr3 = 1;   % pure normal thrust: w=1, exercises K_L*w in Ldot
U3 = [beta3; thr3];
lam3(7) = 50;                                      % nonzero lam_t, exercises the S coupling term
[pv3, S3, info3] = mee_primer_switch(X3, U3, lam3, sigmaNode, dLnode, par);
assert(abs(info3.KL - (-0.001)) < tol, 'Test3: K_L mismatch');
assert(abs(info3.Ldot0 - 1) < tol, 'Test3: Ldot0 mismatch');
LdotExp = 1 + (-0.001)*1;                          % Ldot0 + K_L*w, w=1
assert(abs(info3.Ldot - LdotExp) < tol, 'Test3: Ldot (actual) mismatch');
pelExp = [0; 2; 1.01];
assert(max(abs(info3.pel - pelExp)) < tol, 'Test3: p_el mismatch');
Gexp = Tm*1.01 + 0 + 50*1 + 1;                      % lam(1:5)'dXdt(1:5) + lam_m*mdot(=0) + lam_t*1 + thr
assert(abs(info3.G - Gexp) < tol, 'Test3: G mismatch');
pvExp = (Tm/1)*LdotExp*pelExp - Gexp*(-0.001)*[0;0;1];
assert(max(abs(pv3 - pvExp)) < tol, 'Test3: primerVec (K_L~=0) mismatch');
C1exp = Tm*dot(pelExp, beta3) - (Tm/cEx)*0 + 1;      % lam_m=0 here
Sexp3 = C1exp*1 - 50*(-0.001)*1;                     % C1*Ldot0 - lam_t*K_L*beta3
assert(abs(S3 - Sexp3) < tol, 'Test3: S (K_L~=0) mismatch');
fprintf('Test 3 (K_L~=0 coupling term, hand-checked) PASSED\n');

%% Test 4 -- B(X) numeric extraction cross-check: 5 independent unit-lam_el
% probes at the Test-3 state (K_L~=0) must reproduce each ROW of the
% hand-derived B matrix exactly (pel = B'*lam_el, so lam_el=e_i extracts
% row i of B). This independently cross-checks the closed-form algebra in
% mee_primer_switch.m's header against its own numerical extraction path.
Bexpected = [0,2,0; 0,2,0; -1,0,0; 0,0,0.505; 0,0,0];
X4  = repmat(X3, 1, 5);
U4  = repmat([0;1;0;1], 1, 5);              % arbitrary fixed valid control, unused by pel
lam4 = zeros(7,5);
for row = 1:5, lam4(row,row) = 1; end       % lam4(:,row) = e_row
sigma4 = repmat(sigmaNode, 1, 5);
[~, ~, info4] = mee_primer_switch(X4, U4, lam4, sigma4, dLnode, par);
for row = 1:5
    assert(max(abs(info4.pel(:,row) - Bexpected(row,:).')) < tol, ...
        sprintf('Test4: B row %d mismatch', row));
end
fprintf('Test 4 (B(X) numeric-extraction vs closed-form, 5 rows) PASSED\n');

fprintf('ALL TESTS PASSED (test_verify_pmp_mee)\n');
