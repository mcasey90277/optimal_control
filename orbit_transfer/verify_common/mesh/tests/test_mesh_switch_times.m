% TEST_MESH_SWITCH_TIMES  Unit test for sub-grid switch-time extraction.
%
% The property under test is the one the study's validity rests on: a switch
% must be located STRICTLY INSIDE its bracketing interval, not snapped to a
% node. If extraction quantizes to the grid, every switch-time order the study
% reports is p ~ 1 by construction, confirming hypothesis H1 circularly and
% proving nothing (see the plan's Global Constraint on switch extraction).
%
% All cases are analytic -- no solves, no campaign dependencies.
root = fileparts(fileparts(mfilename('fullpath')));  cd(root);  addpath(root);

% ---- (a) linear Shat crossing at a KNOWN sub-grid point ------------------
% 5 nodes at sg = 0 .. 1 (h = 0.25). Put the zero of Shat at sg = 0.4, which
% is 60% of the way through the bracket [0.25, 0.50] -- deliberately not at a
% node, not at a midpoint, and not at a quarter point.
sg    = (0:0.25:1).';
tNode = sg;                                  % domain IS time here
Shat  = (sg.' - 0.4);                         % + for sg > 0.4
burn  = Shat < 0;                             % burn while Shat < 0 -> flips at 0.4
s = mesh_switch_times(sg, tNode, burn, Shat);

assert(s.n == 1, 'a: expected 1 switch, got %d', s.n);
assert(abs(s.sgSw - 0.4) < 1e-12, 'a: sgSw = 0.4 expected, got %.15g', s.sgSw);
assert(abs(s.tSw  - 0.4) < 1e-12, 'a: tSw = 0.4 expected, got %.15g', s.tSw);
assert(strcmp(s.method, 'switchfn'), 'a: method should be switchfn, got %s', s.method);

% THE ANTI-QUANTIZATION ASSERTION: the recovered location must sit clear of
% BOTH bracketing nodes. A thresholding implementation returns a node and
% fails here.
h = 0.25;
assert(min(abs(s.sgSw - [0.25 0.50])) > 0.1*h, ...
    'a: switch collapsed onto a node (sgSw = %.15g); extraction is grid-quantized', s.sgSw);
assert(s.subGridFrac > 0.05 && s.subGridFrac < 0.95, ...
    'a: subGridFrac = %.4f is at a bracket edge -- estimator collapsed to a node', s.subGridFrac);
assert(abs(s.subGridFrac - 0.6) < 1e-12, 'a: subGridFrac = 0.6 expected, got %.15g', s.subGridFrac);

% ---- (b) NON-UNIFORM physical time -> tSw must follow tNode, not sg ------
% Same sigma grid and same crossing at sg = 0.4, but physical time is
% quadratic in sigma. The crossing is 60% through the bracket, so the
% physical time is the 60% interpolant of [t(2), t(3)] = [0.0625, 0.25]:
%   t = 0.0625 + 0.6*(0.25 - 0.0625) = 0.175
tNodeQ = sg.^2;
s = mesh_switch_times(sg, tNodeQ, burn, Shat);
assert(abs(s.sgSw - 0.4) < 1e-12, 'b: sgSw unchanged at 0.4, got %.15g', s.sgSw);
assert(abs(s.tSw - 0.175) < 1e-12, ...
    'b: tSw = 0.175 expected (interpolated in tNode), got %.15g', s.tSw);
assert(abs(s.tSw - 0.4) > 1e-6, 'b: tSw must NOT equal sgSw when time is nonuniform');

% ---- (c) no switching function -> documented throttle fallback -----------
% Throttle crosses 0.5 between nodes. thr = [1 1 0.2 0 0] puts the 0.5
% crossing between node 2 (1.0) and node 3 (0.2), at fraction
% (1 - 0.5)/(1 - 0.2) = 0.625 -> sg = 0.25 + 0.625*0.25 = 0.40625
thr = [1 1 0.2 0 0];
s = mesh_switch_times(sg, tNode, thr > 0.5, [], thr);
assert(s.n == 1, 'c: expected 1 switch, got %d', s.n);
assert(strcmp(s.method, 'throttle'), 'c: method should be throttle, got %s', s.method);
assert(abs(s.sgSw - 0.40625) < 1e-12, 'c: sgSw = 0.40625 expected, got %.15g', s.sgSw);

% ---- (d) no transition -> empty, not an error ---------------------------
s = mesh_switch_times(sg, tNode, true(1,5), ones(1,5));
assert(s.n == 0, 'd: expected 0 switches, got %d', s.n);
assert(isempty(s.tSw) && isempty(s.sgSw), 'd: arrays must be empty');

% ---- (e) two adjacent switches must NOT be merged -----------------------
% burn = [0 1 0 0 0]: transitions between nodes 1-2 and 2-3. Shat must change
% sign across each bracket for the switchfn path to take both.
sg5   = (0:0.25:1).';
burn2 = [false true false false false];
Shat2 = [ 1 -1  1  2  3];
s = mesh_switch_times(sg5, sg5, burn2, Shat2);
assert(s.n == 2, 'e: expected 2 adjacent switches, got %d', s.n);
assert(all(diff(s.tSw) > 0), 'e: switch times must be strictly increasing');
assert(abs(s.sgSw(1) - 0.125) < 1e-12, 'e: first at 0.125, got %.15g', s.sgSw(1));
assert(abs(s.sgSw(2) - 0.375) < 1e-12, 'e: second at 0.375, got %.15g', s.sgSw(2));

fprintf('test_mesh_switch_times: ALL PASS\n');
