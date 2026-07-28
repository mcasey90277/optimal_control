% TEST_MESH_MATCH_SWITCHES  Unit test for cross-mesh switch pairing.
%
% Switches must be paired across refinement levels by PHYSICAL TIME, because
% indices and counts both differ between meshes: a finer mesh may resolve a
% short arc the coarse one missed entirely. The property that matters is that
% an unmatched switch is REPORTED rather than quietly dropped -- a switch that
% appears or vanishes under refinement is a finding about the solution's
% topology, not a nuisance to be absorbed by the matcher.
%
% All cases are analytic. No solves, no campaign dependencies.
root = fileparts(fileparts(mfilename('fullpath')));  cd(root);  addpath(root);
tol = 1e-9;

% ---- (a) identical sets -> all matched, zero dt -------------------------
m = mesh_match_switches([1 2 3], [1 2 3], 0.1);
assert(isequal(size(m.pairs), [3 2]), 'a: expected 3 pairs, got %dx%d', size(m.pairs,1), size(m.pairs,2));
assert(isempty(m.unmatchedA) && isempty(m.unmatchedB), 'a: nothing should be unmatched');
assert(m.maxAbsDt < tol, 'a: dt should be zero, got %.3e', m.maxAbsDt);

% ---- (b) small perturbation within tol -> matched, SIGNED dt recovered ---
% dt is defined B - A (fine minus coarse), so the signs must come back as
% [+0.01, -0.01, +0.02], not as absolute values.
m = mesh_match_switches([1 2 3], [1.01 1.99 3.02], 0.1);
assert(isempty(m.unmatchedA) && isempty(m.unmatchedB), 'b: nothing should be unmatched');
assert(max(abs(m.dt - [0.01 -0.01 0.02])) < 1e-12, ...
    'b: signed dt = [0.01 -0.01 0.02] expected, got [%.4g %.4g %.4g]', m.dt);
assert(abs(m.maxAbsDt - 0.02) < 1e-12, 'b: maxAbsDt = 0.02 expected, got %.15g', m.maxAbsDt);

% ---- (c) EXTRA switch on the fine mesh -> reported, not silently dropped -
% This is the case the study cares about most: refinement resolved a switch
% the coarse mesh did not have. It must surface as unmatched on B, and the
% unmatched one must be the MIDDLE switch (t = 2), not an endpoint.
m = mesh_match_switches([1 3], [1 2 3], 0.1);
assert(isequal(size(m.pairs), [2 2]), 'c: expected 2 pairs, got %d', size(m.pairs,1));
assert(numel(m.unmatchedB) == 1, 'c: expected 1 unmatched on B, got %d', numel(m.unmatchedB));
assert(isempty(m.unmatchedA), 'c: nothing should be unmatched on A');
assert(m.unmatchedB == 2, 'c: the unmatched B index should be 2 (the t=2 switch), got %d', m.unmatchedB);

% ---- (d) beyond tol -> NOT matched (no stealing a distant partner) -------
m = mesh_match_switches(1, 5, 0.1);
assert(isempty(m.pairs), 'd: no pair should form across a 4.0 gap with tol 0.1');
assert(numel(m.unmatchedA) == 1 && numel(m.unmatchedB) == 1, 'd: both should be unmatched');

% ---- (e) GLOBAL nearest matching, not sequential -------------------------
% tA = [1 2], tB = [1.9 2.15], tol 0.5. A sequential matcher pairs A1 with B1
% (|1-1.9| = 0.9, inside tol) and then A2 with B2. The correct global-minimum
% matcher takes the closest pair first -- A2 with B1 at 0.1 -- after which A1
% has no partner within tol. Each switch is used at most once.
m = mesh_match_switches([1 2], [1.9 2.15], 0.5);
assert(isequal(size(m.pairs), [1 2]), 'e: expected exactly 1 pair, got %d', size(m.pairs,1));
assert(isequal(m.pairs, [2 1]), 'e: expected pair (A2,B1), got (A%d,B%d)', m.pairs(1), m.pairs(2));
assert(isequal(m.unmatchedA, 1), 'e: A1 should be unmatched');
assert(isequal(m.unmatchedB, 2), 'e: B2 should be unmatched');

% ---- (f) empty inputs are a legal degenerate case -----------------------
m = mesh_match_switches([], [1 2], 0.1);
assert(isempty(m.pairs), 'f: no pairs possible against an empty A');
assert(isempty(m.unmatchedA) && numel(m.unmatchedB) == 2, 'f: both B entries unmatched');
assert(isnan(m.maxAbsDt), 'f: maxAbsDt undefined with no pairs -> NaN');

fprintf('test_mesh_match_switches: ALL PASS\n');
