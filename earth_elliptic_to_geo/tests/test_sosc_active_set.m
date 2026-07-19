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
fprintf('test_sosc_active_set PASSED\n');
