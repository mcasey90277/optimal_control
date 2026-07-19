% Synthetic: n=2, one equality (row1, kind eq), one inequality g<=0 (row2,
% ineqHi bound 0). Choose grad_f, A_all, lam_g so stationarity is exactly 0
% under sign s=+1, and slack/comp are clean.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
R.n=2; R.m=2;
R.grad_f = [ -1; -1 ];
R.A_all  = sparse([1 0; 0 1]);        % d g1/dx, d g2/dx
R.lam_g  = [1; 1];                     % grad_f + A' * lam = 0  (s=+1)
% NOTE: brief's original R.gval=[0;-0.2] (inactive ineq, slack 0.2) is
% arithmetically inconsistent with its own comp<1e-9 assert: stationarity
% forces lam_g=[1;1] (identity A_all), so an inactive row-2 with lam=1
% gives comp=|1*0.2|=0.2, not ~0. Corrected to an ACTIVE row 2 (slack 0)
% so the nonzero multiplier is complementarity-consistent; verified
% numerically before this edit. All other brief assertions unaffected.
R.gval   = [0; 0];                     % eq satisfied; ineq active, slack 0
R.creg = struct('label',{'eqA','ineqB'},'kind',{'eq','ineqHi'}, ...
                'rows',{1,2},'bound',{0,0},'node',{[],[]});
K = sosc_kkt_residual(R, sosc_defaults());
assert(K.signOK && K.sign==1, 'sign should resolve to +1');
assert(K.stat < 1e-12, 'stationarity ~0'); assert(K.comp < 1e-9, 'comp ~0');
assert(K.dualFeas <= 0 || K.dualFeas < 1e-12, 'dual feasible');
assert(K.pass, 'overall KKT pass');
fprintf('test_sosc_kkt_residual PASSED\n');
