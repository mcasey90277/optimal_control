function ok = test_lift_space_dim()
%% Purpose:
%
%   Tests lift_space_dim -- the numerical-rank rule behind the abnormal-lift
%   gate (dim S) of mintime_hypothesis_gates.  The catalog-scale pass of
%   2026-09-07 flagged 512 entries with dim S = 0: their smallest singular
%   value sat between the fixed rank tolerance (1e-8 * sv(1)) and the
%   accuracy to which the accepted normal lift itself satisfies the lift
%   constraints (nullResid ~ 1e-7 .. 2.5e-6).  A null space cannot be
%   resolved finer than its known member's own residual, so the rule must
%   be: tol = max(rankTol*sv(1), 10*nullResid), capped at 1e-3*sv(1) so a
%   noisy lift can never inflate the count.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rankTol = 1e-8;
sv6 = [1; 0.7; 0.3; 0.1; 0.05; 0.01];

% (1) golden-cell shape: sv7 far below rankTol -> dim 1
[d, tol] = lift_space_dim([sv6; 1e-10], 1e-10, rankTol);
ok = chk(ok, d == 1 && tol >= 1e-8, sprintf('clean null space: dim = %d (tol %.1e)', d, tol));

% (2) the flagged shape: sv7 = 3e-7 > rankTol, nullResid 4e-7 -> dim 1
[d, tol] = lift_space_dim([sv6; 3e-7], 4e-7, rankTol);
ok = chk(ok, d == 1, sprintf('flagged shape (sv7 3e-7, nullResid 4e-7): dim = %d (tol %.1e)', d, tol));

% (3) two-dimensional null space -> dim 2 (a genuine abnormal lift)
[d, ~] = lift_space_dim([sv6(1:5); 2e-9; 1e-10], 1e-10, rankTol);
ok = chk(ok, d == 2, sprintf('two small singular values: dim = %d', d));

% (4) noisy lift cannot inflate: nullResid 1e-2, sv7 = 1e-3 -> unresolved (0)
[d, tol] = lift_space_dim([sv6; 1e-3], 1e-2, rankTol);
ok = chk(ok, d == 0 && tol <= 1e-3, sprintf('noisy lift: dim = %d (tol capped at %.1e)', d, tol));

% (5) the gap is reported: sv7/sv6 for shape (2)
[~, ~, gap] = lift_space_dim([sv6; 3e-7], 4e-7, rankTol);
ok = chk(ok, abs(gap - 3e-5) < 1e-12, sprintf('gap sv7/sv6 = %.1e', gap));

if ok, fprintf('TEST_LIFT_SPACE_DIM: ALL PASS\n');
else,  fprintf('TEST_LIFT_SPACE_DIM: FAILURE (see lines above)\n');
end
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
