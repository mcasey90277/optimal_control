% test_sosc_10N  Integration: full SOSC certificate on the 10 N MEE min-fuel row.
%
% The payoff of Amendment B (reduced-Hessian inertia + WEAK_MIN verdict): the
% 10 N bang-bang min-fuel solution is a WEAK (non-strict) local minimum -- the
% reduced Hessian is PSD (no negative curvature) with a large flat manifold
% (~270 directions preserving the switching structure). The certificate must
% return verdict=WEAK_MIN, status=certified-weak-min, red.nneg==0, nFlat>0.
%
% Sec 11.4 hardening (2026-07-19): the 10 N KKT (nk=3885) is at/below
% tol.maxEigDim, so sosc_inertia now uses gold-standard eig(full(K)) instead
% of the unreliable ldl-pivot-sign inertia. This must reproduce the true
% reduced inertia (116,0,270) robustly (IN.method=='eig', IN.robust==true),
% with nFlat now the TRUE value (~270, not the ldl-workaround's ~265-269).
% This is a real IPOPT warm re-solve (~1-2 min).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
sosc = verify_sosc_mee(fullfile(module_root(),'results','MEE_M2_10N.mat'));
fprintf(['10 N verdict=%s reason="%s"\n  stat=%.2e drift=%.2e sign=%+d\n' ...
    '  primalEq=%.2e primalIneq=%.2e dualFeas=%.2e comp=%.2e\n' ...
    '  m_active=%d nEq=%d nStrong=%d nWeak=%d licq=%d\n' ...
    '  KKT inertia=[%d %d %d]  rankA=%d  method=%s robust=%d\n' ...
    '  reduced inertia red=[%d %d %d] redConsistent=%d  nFlat=%d\n'], ...
    sosc.verdict, sosc.reason, sosc.kkt.stat, sosc.drift, sosc.sign, ...
    sosc.kkt.primalEq, sosc.kkt.primalIneq, sosc.kkt.dualFeas, sosc.kkt.comp, ...
    sosc.active.m_active, sosc.active.nEq, sosc.active.nStrong, sosc.active.nWeak, ...
    sosc.active.licq, sosc.inertia.npos, sosc.inertia.nneg, sosc.inertia.nzero, ...
    sosc.inertia.rankA, sosc.inertia.method, sosc.inertia.robust, ...
    sosc.red.npos, sosc.red.nneg, sosc.red.nzero, ...
    sosc.inertia.redConsistent, sosc.nFlat);

% Orchestrator correctness
assert(isfield(sosc,'verdict') && ischar(sosc.verdict), 'well-formed verdict');
assert(any(strcmp(sosc.verdict,{'PASS','WEAK_MIN','FAIL','INCONCLUSIVE','ERROR'})), 'verdict domain');
assert(strcmp(sosc.status, sosc_decide(sosc.kkt,sosc.active,sosc.inertia).status), ...
    'status consistent with sosc_decide');
assert(sosc.kkt.stat  < sosc.thresholds.stat,  'stationarity machine-tight');
assert(sosc.drift     < sosc.thresholds.drift, 'warm-resolve drift machine-tight');

% The payoff (Amendment B acceptance gate): 10 N certifies as WEAK_MIN.
assert(strcmp(sosc.verdict,'WEAK_MIN'), ...
    '10 N expected WEAK_MIN, got %s (%s)', sosc.verdict, sosc.reason);
assert(strcmp(sosc.status,'certified-weak-min'), ...
    '10 N status expected certified-weak-min, got %s', sosc.status);
assert(sosc.red.nneg==0, '10 N reduced Hessian must have no negative curvature (red.nneg==0), got %d', sosc.red.nneg);
assert(sosc.nFlat>0, '10 N weak-min must have a flat manifold (nFlat>0), got %d', sosc.nFlat);

% Sec 11.4 hardening gate: KKT (nk=3885) is within tol.maxEigDim -> must use
% the gold-standard eig path, and must reproduce the true reduced inertia
% (116,0,270) -- the TRUE nFlat, not the ldl-workaround's ~265-269.
assert(strcmp(sosc.method,'eig'), '10 N must use eig inertia method, got %s', sosc.method);
assert(sosc.robust==true, '10 N inertia must be robust (eig, size-guarded), got %d', sosc.robust);
assert(abs(sosc.nFlat-270)<=3, ...
    '10 N nFlat expected ~270 (true eig value), got %d', sosc.nFlat);

fprintf('test_sosc_10N PASSED: certified-weak-min, method=%s robust=%d nFlat=%d\n', ...
    sosc.method, sosc.robust, sosc.nFlat);
