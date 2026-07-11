function test_warmstart_on_mesh()
% TEST_WARMSTART_ON_MESH  Originals exact; insert straddling a switch holds left.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
sigma = linspace(0, 1, 6).';                 % 5 intervals, 6 nodes
X = [ (1:6);  (11:16);  (21:26);  (31:36); ...
      (41:46); (51:56);  linspace(1,0.9,6);  linspace(0,1,6) ];   % [8x6]
% throttle switches between node 3 and node 4 (burn -> coast)
s   = [1 1 1 0 0 0];
al  = repmat([1;0;0], 1, 6);                 % unit +x direction
U   = [al; s];                                % [4x6]
out = struct('X', X, 'U', U);

% refine the interval that straddles the switch (interval 3: nodes 3-4)
sigmaNew = sort([sigma; 0.5*(sigma(3)+sigma(4))]);
isNew    = false(size(sigmaNew));  isNew(abs(sigmaNew - 0.5*(sigma(3)+sigma(4))) < 1e-15) = true;

[X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew);

assert(isequal(size(X0), [8 7]) && isequal(size(U0), [4 7]), 'sizes');
% originals preserved exactly
origCols = find(~isNew);
assert(max(abs(X0(:, origCols) - X), [], 'all') < 1e-15, 'original X preserved');
assert(max(abs(U0(:, origCols) - U), [], 'all') < 1e-15, 'original U preserved');
% inserted throttle holds the LEFT (pre-switch, burn) value = 1, never averaged to 0.5
ins = find(isNew);
assert(abs(U0(4, ins) - 1) < 1e-15, 'insert throttle must step-hold left value 1, got %g', U0(4, ins));
% inserted direction is unit-norm
assert(abs(norm(U0(1:3, ins)) - 1) < 1e-12, 'insert direction must be unit norm');
% inserted state is between its neighbors (pchip monotone here)
assert(X0(1, ins) > X(1,3) && X0(1, ins) < X(1,4), 'insert state interpolated');
fprintf('ALL PASS\n');
end
