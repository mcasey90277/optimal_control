% TEST_MESH_ORDER  Unit test for observed-order / Richardson estimation.
%
% Synthetic sequences with analytically known answers, so a failure here is a
% bug in the estimator and never an artifact of a solve.
%
% The cases that matter most are the DEGENERATE ones: the estimator must
% refuse to report an order when the differences flip sign or stop shrinking
% (a basin/topology change between levels), rather than emitting a
% plausible-looking number. A confidently wrong p is worse than a NaN here,
% because the study's conclusions are read off these values.
root = fileparts(fileparts(mfilename('fullpath')));  cd(root);  addpath(root);

% ---- (a) exact 2nd-order sequence about a known limit: f_k = L + C*h^2 ---
L0 = 7;  C = 3;  h = [1, 1/2, 1/4];
o = mesh_order(L0 + C*h.^2, [1 2 4]);
assert(abs(o.p - 2) < 1e-9, 'a: p = 2 expected, got %.9f', o.p);
assert(abs(o.rich - L0) < 1e-9, 'a: Richardson should recover %g, got %.9f', L0, o.rich);
assert(o.monotone, 'a: differences must be shrinking');

% ---- (b) exact 1st-order sequence ---------------------------------------
o = mesh_order(L0 + C*h, [1 2 4]);
assert(abs(o.p - 1) < 1e-9, 'b: p = 1 expected, got %.9f', o.p);
assert(abs(o.rich - L0) < 1e-9, 'b: Richardson should recover %g, got %.9f', L0, o.rich);

% ---- (c) two levels only -> order undefined, deltas still reported -------
o = mesh_order([1 1.5], [1 2]);
assert(isnan(o.p) && isnan(o.rich), 'c: an order needs 3 levels');
assert(abs(o.dLast - 0.5) < 1e-12, 'c: dLast = 0.5 expected, got %.15g', o.dLast);
assert(abs(o.rel - 0.5/1.5) < 1e-12, 'c: rel = dLast/|last|, got %.15g', o.rel);

% ---- (d) non-monotone (noise / basin change) is FLAGGED, not ordered -----
o = mesh_order([1, 1.5, 1.4], [1 2 4]);
assert(~o.monotone, 'd: must flag non-monotone');
assert(isnan(o.p), 'd: must REFUSE an order across a sign flip, got p = %.6f', o.p);
assert(isnan(o.rich), 'd: must refuse Richardson across a sign flip');

% ---- (e) differences shrinking but NOT converging is still an order ------
% A clean 4th-order sequence: p should come back as 4, not be clamped.
o = mesh_order(L0 + C*h.^4, [1 2 4]);
assert(abs(o.p - 4) < 1e-9, 'e: p = 4 expected, got %.9f', o.p);

% ---- (f) a stalled sequence (zero last difference) -> no order ----------
% v3 == v2 makes the difference ratio infinite; the estimator must return NaN
% rather than Inf-propagating into the report.
o = mesh_order([1, 1.5, 1.5], [1 2 4]);
assert(isnan(o.p), 'f: a zero last difference cannot support an order');
assert(o.dLast == 0, 'f: dLast should be exactly 0, got %.15g', o.dLast);

% ---- (g) four levels: sliding three-level slopes are reported ------------
% The plan requires four levels so that the CONSISTENCY of successive slopes
% is itself the evidence of an asymptotic regime. A clean 2nd-order sequence
% must give p = 2 on both windows.
h4 = [1, 1/2, 1/4, 1/8];
o = mesh_order(L0 + C*h4.^2, [1 2 4 8]);
assert(abs(o.p - 2) < 1e-9, 'g: overall p = 2 expected, got %.9f', o.p);
assert(numel(o.pWindows) == 2, 'g: expected 2 sliding windows, got %d', numel(o.pWindows));
assert(max(abs(o.pWindows - 2)) < 1e-9, 'g: both windows should read 2, got [%.6f %.6f]', o.pWindows);
assert(o.windowsConsistent, 'g: identical window slopes must read as consistent');

% ---- (h) non-constant refinement ratio is rejected ----------------------
% p = log(ratio)/log(r) assumes ONE r. Factors [1 2 3] have r = 2 then 1.5.
ok = false;
try
    mesh_order([3 2 1.5], [1 2 3]);
catch err
    ok = contains(err.identifier, 'ratio') || contains(lower(err.message), 'ratio');
end
assert(ok, 'h: a non-constant refinement ratio must be rejected, not averaged');

% ---- (i) a window spanning a BRANCH CHANGE is refused ------------------
% The regression this guards: the 1 N study computed a window across a
% 173 -> 171 switch-count change and then cited its agreement with a valid
% window as evidence of an asymptotic regime. The policy existed only as a
% header comment and nothing enforced it. Now it is arithmetic.
vals = L0 + C*h4.^2;                      % a clean 2nd-order series
o = mesh_order(vals, [1 2 4 8], [false true true]);   % first difference invalid
assert(isnan(o.pWindows(1)), 'i: a window using an invalid difference must be NaN');
assert(abs(o.pWindows(2) - 2) < 1e-9, 'i: the clean window should still read 2');
assert(o.nValidWindows == 1, 'i: expected 1 valid window, got %d', o.nValidWindows);
assert(~o.windowsConsistent, ...
    'i: ONE window cannot be "consistent" -- it agrees with itself trivially');
assert(~o.monotone, 'i: monotonicity is not claimable across a branch change');
assert(abs(o.p - 2) < 1e-9, 'i: p must come from the finest VALID window');

% (j) every window invalid -> no order at all
o = mesh_order(vals, [1 2 4 8], [false true false]);
assert(isnan(o.p) && o.nValidWindows == 0, 'j: no valid window -> no order');

fprintf('test_mesh_order: ALL PASS\n');
