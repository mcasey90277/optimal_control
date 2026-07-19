% TEST_WARMSTART_MEE  No-solve-where-possible coverage for the warm-start
% mesh-refine machinery (review finding, Fix 2): cfg.warmStart
% (run_transfer_mee.m) and cfg.warmStartAnchor (run_mintime_mee.m) are
% load-bearing for every remaining thrust-ladder rung but had ZERO tests
% before this file. Both call sites now share interp_warmstart.m (factored
% out of what were two near-identical inline interp1 blocks), so this file
% tests that ONE pure function directly instead of driving either caller's
% full NLP machinery.
%
% Part 1 (no solve): interp_warmstart.m unit tests -- sizes on the
% destination grid, endpoint preservation, throttle-row bang-bang
% preservation (nearest interp), and the beta-row unit-norm renormalization
% fix (LATENT BUG check, re-review finding: the pre-refactor inline code did
% NOT renormalize after linearly interpolating the RTN thrust-direction unit
% vectors, so an interpolated midpoint between two well-separated unit
% directions came out sub-unit -- a real violation of casadi_lt_mee.m's own
% |beta|=1 constraint, confirmed live below, NOT a hypothetical).
%
% Part 2 (one cheap live solve, explicitly budgeted): warm-starts a 10 N
% fuel re-solve at N=147 (~20/rev; the brief's illustrative "N=155" target
% assumed a slightly different source revolution count than the certified
% MEE_M2_10N.mat's actual 7.326 revs actually gives at 20 nodes/rev) from
% the certified 25/rev fuel solution via cfg.warmStart, maxIter=5 -- mirrors
% test_mee_solver_smoke.m's non-convergence-gated pattern: assert the full
% out struct returns with no exception and the warm-start seed actually
% entered the solve (N/report fields populated), NOT that it certifies.
%
% REFERENCES: [1] interp_warmstart.m (function under test in Part 1).
%   [2] run_transfer_mee.m (cfg.warmStart, exercised live in Part 2).
%   [3] run_mintime_mee.m (cfg.warmStartAnchor, the OTHER caller of
%   interp_warmstart.m -- covered by Part 1 since it is the same function).
%   [4] test_mee_solver_smoke.m (the no-certification-required smoke
%   pattern Part 2 mirrors). [5] task-7-report.md review findings (Fix 2).

root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;

% =============================================================================
% Part 1: interp_warmstart.m -- pure function, no solve
% =============================================================================

% --- synthetic source trajectory: 5 nodes, beta rotating through 4 distinct
% RTN directions (deliberately NOT collinear, so LINEAR interpolation between
% adjacent columns produces a genuinely sub-unit midpoint -- this is what
% makes the renormalization check below a real test, not a vacuous one) ------
sigmaSrc = [0; 0.25; 0.50; 0.75; 1.00];
betaSrc  = [ 1  0  0 -1  0;      % R
             0  1  0  0 -1;      % T
             0  0  1  0  0];     % N
assert(max(abs(sqrt(sum(betaSrc.^2,1)) - 1)) < 1e-14, ...
    'test_warmstart_mee: fixture bug -- betaSrc columns must be unit vectors');
thrSrc   = [1 1 0 0 1];          % bang-bang, only {0,1}
Xsrc = zeros(7,5);
for r = 1:7
    Xsrc(r,:) = 10*sigmaSrc.' + r;    % distinct, easily-checked per-row values
end
Usrc = [betaSrc; thrSrc];
dLsrc = 5.0;

% destination grid: 9 nodes, chosen to land EXACTLY on the source midpoints
% (0.125, 0.375, ...) so the beta-dip case is sampled deterministically.
sigmaDst = linspace(0, 1, 9).';

W = interp_warmstart(Xsrc, Usrc, dLsrc, sigmaSrc, sigmaDst);

% --- (a) sizes on the destination grid --------------------------------------
assert(isequal(size(W.X), [7 9]), 'test_warmstart_mee: W.X size mismatch, got %s', mat2str(size(W.X)));
assert(isequal(size(W.U), [4 9]), 'test_warmstart_mee: W.U size mismatch, got %s', mat2str(size(W.U)));

% --- (b) endpoint preservation: X(:,1)/X(:,end) exact -----------------------
assert(max(abs(W.X(:,1)   - Xsrc(:,1)))   < 1e-12, 'test_warmstart_mee: W.X(:,1) must exactly match Xsrc(:,1)');
assert(max(abs(W.X(:,end) - Xsrc(:,end))) < 1e-12, 'test_warmstart_mee: W.X(:,end) must exactly match Xsrc(:,end)');

% --- (c) throttle row (nearest): every destination value is a MEMBER of the
% source's value set -- bang-bang preservation, no blurred intermediate
% throttle values ------------------------------------------------------------
srcThrSet = unique(thrSrc);
for k = 1:numel(W.U(4,:))
    d_ = min(abs(W.U(4,k) - srcThrSet));
    assert(d_ < 1e-12, ['test_warmstart_mee: destination throttle value %.6f at node %d ' ...
        'is not a member of the source value set %s (nearest interp should never blur it)'], ...
        W.U(4,k), k, mat2str(srcThrSet));
end

% --- (d) LATENT BUG CHECK: confirm the un-renormalized linear interpolant
% really would be sub-unit at the deliberate dip (sigmaDst=0.125, between
% betaSrc columns 1=[1;0;0] and 2=[0;1;0] -- linear midpoint [0.5;0.5;0],
% norm 1/sqrt(2)~0.7071), THEN confirm interp_warmstart's actual output at
% that same node IS unit-norm (the fix). This proves the bug was real (not
% hypothetical) and that the fix addresses it. ------------------------------
idxDip = find(abs(sigmaDst - 0.125) < 1e-12);
assert(isscalar(idxDip), 'test_warmstart_mee: fixture bug -- expected sigmaDst to contain 0.125 exactly');
rawBetaDip = interp1(sigmaSrc, betaSrc.', sigmaDst(idxDip), 'linear').';
rawNormDip = norm(rawBetaDip);
assert(rawNormDip < 0.99, ['test_warmstart_mee: LATENT BUG CHECK fixture failed to reproduce a ' ...
    'sub-unit raw linear interpolant (got norm=%.4f, expected < 0.99) -- the dip is not being ' ...
    'exercised, this test is not actually checking the renormalization fix'], rawNormDip);
fprintf(['test_warmstart_mee: LATENT BUG CONFIRMED REAL -- raw (un-renormalized) linear ' ...
    'interpolation of two unit RTN thrust directions at the dip node gives |beta|=%.4f ' ...
    '(sub-unit, violates casadi_lt_mee.m''s beta(1,k)^2+beta(2,k)^2+beta(3,k)^2==1 constraint); ' ...
    'interp_warmstart.m''s renormalization fixes this (checked next).\n'], rawNormDip);

% --- (e) beta rows (1-3) of the ACTUAL W.U output are unit-norm at EVERY
% destination node (the fix, applied globally, not just at the dip) --------
betaNorms = sqrt(sum(W.U(1:3,:).^2, 1));
assert(max(abs(betaNorms - 1)) < 1e-10, ['test_warmstart_mee: interp_warmstart.m output beta ' ...
    'rows are not unit-norm at every destination node (max deviation %.3e) -- renormalization ' ...
    'fix is not working'], max(abs(betaNorms - 1)));

% --- (f) dL passthrough ------------------------------------------------------
assert(W.dL == dLsrc, 'test_warmstart_mee: W.dL must exactly pass through dLsrc (got %.6f, expected %.6f)', ...
    W.dL, dLsrc);

fprintf(['test_warmstart_mee: Part 1 (interp_warmstart.m, no solve) ALL PASS -- sizes correct, ' ...
    'endpoints exact, throttle row bang-bang-preserved, beta rows renormalized to unit norm ' ...
    '(latent bug confirmed real and fixed), dL passthrough correct\n']);

% =============================================================================
% Part 2: ONE cheap live smoke -- run_transfer_mee.m's cfg.warmStart path,
% maxIter=5, budget mirrors test_mee_solver_smoke.m (no certification
% required, just "returns without exception, seed actually entered")
% =============================================================================
resDir = fullfile(module_root(), 'results');
srcFile = fullfile(resDir, 'MEE_M2_10N.mat');
assert(isfile(srcFile), ['test_warmstart_mee: prerequisite %s not found -- this is the ' ...
    'certified 25/rev 10 N fuel anchor Part 2 warm-starts FROM; if it was deleted, re-run ' ...
    'run_transfer_mee(struct(''thrustN'',10)) to regenerate it first'], srcFile);
Ssrc = load(srcFile);
srcRes = Ssrc.res;

cfg = struct();
cfg.thrustN     = 10;
cfg.ctf         = 1.5;
cfg.tfMinAnchor = 22.2248;
cfg.tag         = 'TEST_warmstart_smoke_T10N';
cfg.nodesPerRev = 20;
cfg.maxIter     = 5;
cfg.warmStart   = struct('sigma', srcRes.sigma, 'X', srcRes.fuel.X, ...
                          'U', srcRes.fuel.U, 'dL', srcRes.fuel.dL);

tStart = tic;
out = run_transfer_mee(cfg);
wallSec = toc(tStart);

budget_s = 120;   % same budget as test_mee_solver_smoke.m
assert(wallSec < budget_s, ['test_warmstart_mee: Part 2 warm-start smoke took %.1fs, ' ...
    'budget %.0fs'], wallSec, budget_s);

% --- the full out struct returned with no exception -- check it actually
% has the expected shape and that the warm-start seed genuinely entered the
% solve (N derived from warmStart.dL/(2*pi) at cfg.nodesPerRev=20, NOT a cold
% mee_seed N, and revs close to the source's 7.326, not some unrelated value)
assert(isfield(out, 'report') && isfield(out, 'fuel') && isfield(out, 'seed'), ...
    'test_warmstart_mee: run_transfer_mee output missing expected top-level fields');
assert(isfield(out.report, 'certified'), 'test_warmstart_mee: out.report missing .certified');
Nachieved = size(out.fuel.X, 2) - 1;
assert(Nachieved > 100, ['test_warmstart_mee: warm-started N=%d looks too small for a ' ...
    '~7.3-rev source at 20 nodes/rev -- warm-start seed may not have entered correctly'], Nachieved);
assert(abs(out.report.revs - 7.326) < 0.5, ['test_warmstart_mee: warm-started revs=%.3f is far ' ...
    'from the source solution''s 7.326 revs -- warm-start seed may not have entered correctly'], ...
    out.report.revs);

fprintf(['test_warmstart_mee: Part 2 (run_transfer_mee.m cfg.warmStart, live, maxIter=5) ALL ' ...
    'PASS -- wallSec=%.2f, N achieved=%d (target ~155 per the brief, actual %d from the source''s ' ...
    'true 7.326 revs at 20/rev), report.certified=%d (not required), revs=%.3f, no exception\n'], ...
    wallSec, Nachieved, Nachieved, out.report.certified, out.report.revs);

fprintf('test_warmstart_mee: ALL PASS (Part 1 + Part 2)\n');
