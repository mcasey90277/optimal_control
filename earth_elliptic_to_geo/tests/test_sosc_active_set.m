% n=3 vars; rows: 1 eq (active), 2 ineq. Ineq A active+strong (slack 0, |lam|
% big); ineq B active+WEAK (slack 0, |lam|~0). Expect nWeak==1, A has eq+strong.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
R.n=3; R.m=3;
R.A_all = sparse([1 0 0; 0 1 0; 0 0 1]);
R.gval  = [0; 0; 0];                  % all active (slack 0)
R.lam_g = [5; 3; 1e-12];              % eq; strong ineq; weak ineq
R.creg = struct('label',{'eqA','strongB','weakC'},'kind',{'eq','ineqHi','ineqHi'}, ...
    'rows',{1,2,3},'bound',{0,0,0},'node',{[],137,204});
K.sign = 1;
AS = sosc_active_set(R, K, sosc_defaults());
assert(AS.nEq==1 && AS.nStrong==1 && AS.nWeak==1, 'counts');
assert(AS.m_active==2, 'A = eq + strong only');
assert(any(contains(AS.weakLabels,'204')), 'weak label names the node');
assert(AS.licq, 'independent rows -> LICQ ok');

% Case 2: inactive ineqLo row excluded from both nWeak and AS.A, plus an
% ineqLo slack branch exercised (slack = g-bound) alongside a weak ineqHi row.
R2.n=4; R2.m=4;
R2.A_all = sparse(eye(4));
R2.gval  = [0; 0; 0.5; 0];             % row3 ineqLo INACTIVE (slack = g-bound = 0.5)
R2.lam_g = [1; 4; 0; 1e-11];           % eq; strong ineqLo; inactive; weak ineqHi
R2.creg = struct('label',{'eqA','loStrong','loInactive','hiWeak'}, ...
    'kind',{'eq','ineqLo','ineqLo','ineqHi'},'rows',{1,2,3,4}, ...
    'bound',{0,0,0,0},'node',{[],[],[],204});
K2.sign=1;
AS2 = sosc_active_set(R2, K2, sosc_defaults());
assert(AS2.nEq==1 && AS2.nStrong==1 && AS2.nWeak==1, 'case2: counts (eq/strong/weak)');
assert(AS2.m_active==2, 'case2: A = eq + strong only (inactive ineqLo row3 excluded)');
assert(numel(AS2.weakLabels)==1 && any(contains(AS2.weakLabels,'204')), 'case2: only the weak row labeled, names its node');
assert(AS2.licq, 'case2: independent rows -> licq true');

% Case 3: LICQ FALSE via rank-deficient active Jacobian (duplicate eq rows).
R3.n=2; R3.m=2;
R3.A_all = sparse([1 0; 1 0]);         % duplicate rows -> structural rank 1
R3.gval  = [0; 0];
R3.lam_g = [1; 1];
R3.creg = struct('label',{'eq1','eq2'},'kind',{'eq','eq'}, ...
    'rows',{1,2},'bound',{0,0},'node',{[],[]});
K3.sign=1;
AS3 = sosc_active_set(R3, K3, sosc_defaults());
assert(AS3.nEq==2 && AS3.nStrong==0 && AS3.nWeak==0, 'case3: two eq rows, no ineq');
assert(AS3.m_active==2, 'case3: both eq rows active');
assert(~AS3.licq, 'case3: rank-deficient active Jacobian -> licq false');

fprintf('test_sosc_active_set PASSED\n');
