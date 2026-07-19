% TEST_MEE_THREADING  Task 3: opts.xf forwarding through homotopy_mee's eps-
% step sweep, plus the schema-older cache-fingerprint WARN (a cached step
% whose .fp lacks a field present in the CURRENT config's fingerprint -- e.g.
% the newly-added xf -- must WARN and be trusted, not hard-error). Task 4
% extends this file with further threading tests.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

%% (a) homotopy_mee forwards opts.xf into every eps-step casadi_lt_mee call,
% and records it in fp.xf. Behavioral, not source-grep: casadi_lt_mee
% computes termErr = norm(X(1:5,end) - xf), so recomputing that norm from
% the RETURNED trajectory against OUR custom (non-default) xf can only match
% out.termErr if xf actually reached casadi_lt_mee -- a stale/default xf
% inside the solver would produce a different termErr. maxIter=5 (mirrors
% test_mee_solver_smoke.m's non-convergence-gated smoke pattern): we only
% need the eps-step to RUN and cache, not converge.
par = kepler_lt_params(10, 1500, 2000);
seedOpts = struct('thr', 0.5, 'betaMode', 'transverse', 'N', 10, 'nRev', 1);
[sg, X0, U0, dL0, sinfo] = mee_seed(par, seedOpts);
x0 = X0(:, 1);
xfCustom = [0.80; 0.02; -0.01; 0.005; -0.005];

tmpA = tempname; mkdir(tmpA);
best = homotopy_mee(sg, X0, U0, dL0, struct('par', par, 'x0', x0, ...
    'tfTarget', 1.3*sinfo.tEnd, 'resDir', tmpA, 'tag', 'xffwd', ...
    'sched', [1], 'maxIter', 5, 'printLevel', 0, 'xf', xfCustom));

termErrCheck = norm(best.X(1:5, end) - xfCustom);
assert(abs(best.termErr - termErrCheck) < 1e-9, ...
    'homotopy_mee must forward opts.xf into casadi_lt_mee (termErr mismatch: got %.6e, expected %.6e)', ...
    best.termErr, termErrCheck);

S = load(fullfile(tmpA, 'xffwd_step01.mat'));
assert(isfield(S.fp, 'xf'), 'homotopy_mee must record opts.xf in fp.xf');
assert(isequal(S.fp.xf, xfCustom(:)), 'cached fp.xf must equal the passed opts.xf');

fprintf('test_mee_threading (a) xf-forwarding PASSED (termErr=%.3e)\n', best.termErr);

%% (b) schema-older cache: a step .mat whose fp lacks a field present in the
% CURRENT fp (here: xf) must WARN (id 'homotopy_mee:fpSchemaOlder') and be
% trusted as compatible, NOT hard-error with 'fingerprintMismatch' -- mirrors
% run_mintime_mee.m's check_cache_fp_mt pattern (already shipped there). The
% cache's fp.sched is set to match the resolved schedule exactly so the ONLY
% discrepancy versus the current fp is the missing xf field (a field present
% on both sides with a differing value must still hard-error -- not tested
% here, that is the pre-existing behavior this task must not weaken).
tmpB = tempname; mkdir(tmpB);
schedB = [0];                                   % single step -> ..._step01.mat
fp = struct('sched', schedB);                   % SCHEMA-OLDER: no .xf field
o  = struct('maxDefect', 1e-12, 'switches', 0, 'edge', 1, 'm_f_kg', 1400); %#ok<NASGU>
ok = true; e = 0;                                                          %#ok<NASGU>
Xk = zeros(7, 2); Uk = zeros(4, 2); dLk = 1;                               %#ok<NASGU>
save(fullfile(tmpB, 'thr_step01.mat'), 'o', 'ok', 'Xk', 'Uk', 'dLk', 'e', 'fp');

lastwarn('');
caughtErr = [];
try
    homotopy_mee((0:1).', zeros(7, 2), zeros(4, 2), 1, struct('par', par, ...
        'x0', zeros(7, 1), 'tfTarget', 30, 'resDir', tmpB, 'tag', 'thr', ...
        'sched', schedB, 'xf', [1; 0; 0; 0; 0]));
catch ME
    caughtErr = ME;
end
assert(isempty(caughtErr) || ~contains(caughtErr.identifier, 'fingerprintMismatch'), ...
    'schema-older cache (missing xf) must not hard-error with fingerprintMismatch');
[~, wid] = lastwarn();
assert(strcmp(wid, 'homotopy_mee:fpSchemaOlder'), ...
    'expected homotopy_mee:fpSchemaOlder WARN, got ''%s''', wid);

fprintf('test_mee_threading (b) schema-older-WARN PASSED\n');

%% Task 4: default cfg reuses the existing certified 10 N caches with
% NO fingerprint error (schema-older WARN only), and cfg.xf is carried into fp.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
src1 = fileread(fullfile(module_root(),'drivers','run_transfer_mee.m'));
src2 = fileread(fullfile(module_root(),'drivers','run_mintime_mee.m'));
assert(contains(src1,'xf') && contains(src1,'initElems'), 'run_transfer_mee must thread xf/initElems');
assert(contains(src2,'xf') && contains(src2,'initElems'), 'run_mintime_mee must thread xf/initElems');
% default-preserving reuse: loading the cached 10 N mintime anchor with the
% new (xf-bearing) fp must NOT throw fingerprintMismatch.
lastwarn('');
try
  out = run_mintime_mee(10, 25);      % reuses results/MEE_mintime_T100.mat
  assert(abs(out.tfmin-22.2206) < 1e-2, 'cached 10 N anchor tfmin preserved');
catch ME
  assert(isempty(strfind(ME.identifier,'fingerprintMismatch')), ...
     'default 10 N reuse must not fingerprint-error after adding xf/initElems');
  rethrow(ME);
end
fprintf('test_mee_threading (Task4) PASSED\n');

%% Task 4, Step 5 regression note: a default run_transfer_mee cfg must not
% throw a fingerprint error either -- it may reuse results/MEE_M2_10N.mat (or
% its _seed/_seed_probe/_step* intermediates) OR re-solve to m_f ~ 1377 kg.
lastwarn('');
try
  resT = run_transfer_mee(struct('thrustN', 10, 'ctf', 1.5, 'tfMinAnchor', 22.2248));
  assert(resT.report.certified, 'default run_transfer_mee cfg must certify');
catch ME
  assert(isempty(strfind(ME.identifier,'fingerprintMismatch')), ...
     'default run_transfer_mee reuse must not fingerprint-error after adding xf/initElems');
  rethrow(ME);
end
fprintf('test_mee_threading (Task4 Step5 regression) PASSED (m_f=%.2f kg)\n', resT.report.m_f_kg);

fprintf('test_mee_threading: ALL PASS\n');
