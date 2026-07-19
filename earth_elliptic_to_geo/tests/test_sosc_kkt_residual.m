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

% Case 2 -- discriminating PASS: exercises ineqLo (active) + an inactive
% ineqHi row with nonzero slack, so a swapped kind/row mask or a mishandled
% ineqLo branch would show up as a nonzero comp or a spurious violation.
R2.n=3; R2.m=3;
R2.grad_f = [-1; -2; 0];
R2.A_all  = sparse(eye(3));
R2.lam_g  = [1; 2; 0];                 % eq; active ineqLo; inactive ineqHi (lam 0)
R2.gval   = [0; 0; -0.5];              % row3 inactive: slack 0.5
R2.creg = struct('label',{'eqA','loB','hiC'},'kind',{'eq','ineqLo','ineqHi'}, ...
                 'rows',{1,2,3},'bound',{0,0,0},'node',{[],[],[]});
K2 = sosc_kkt_residual(R2, sosc_defaults());
assert(K2.sign==1 && K2.stat < 1e-12, 'case2: sign +1, stationarity ~0');
assert(K2.primalIneq < 1e-12 && K2.comp < 1e-12, 'case2: feasible + complementary');
assert(K2.pass, 'case2: consistent KKT point must pass');

% Case 3 -- s=-1 resolved AND dual-infeasible FAIL path: pins down the
% resolved-sign branch that picks s=-1, and confirms K.pass goes false when
% dual feasibility is violated under that sign.
R3.n=1; R3.m=1;
R3.grad_f = 1;
R3.A_all  = sparse(1);
R3.lam_g  = 1;                          % with grad_f=1, A=1: rP=2, rM=0 -> s=-1
R3.gval   = 0;                          % active ineqHi
R3.creg = struct('label',{'hiOnly'},'kind',{'ineqHi'},'rows',{1},'bound',{0},'node',{[]});
K3 = sosc_kkt_residual(R3, sosc_defaults());
assert(K3.sign==-1, 'case3: s must resolve to -1 (rM<rP)');
assert(K3.dualFeas > sosc_defaults().dual, 'case3: s*lam=-1 is dual-infeasible');
assert(~K3.pass, 'case3: dual-infeasible point must NOT pass');

fprintf('test_sosc_kkt_residual PASSED\n');
