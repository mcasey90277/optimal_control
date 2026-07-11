function test_refine_sigma()
% TEST_REFINE_SIGMA  Bisection preserves originals, marks inserts, honors guards.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
sigma = linspace(0, 1, 11).';        % 10 intervals
score = zeros(1, 10);  score(3) = 5;  score(7) = 9;   % two hot intervals
opts  = struct('K', 2, 'hFloor', 1e-9, 'maxAdd', 40);
[sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts);

assert(numel(sigmaNew) == 13, 'expected 11+2 nodes, got %d', numel(sigmaNew));
assert(nnz(isNew) == 2, 'exactly 2 inserted');
assert(issorted(sigmaNew), 'sigmaNew sorted');
assert(nDropped == 0, 'no drops expected');
% every original node still present (exact)
for k = 1:numel(sigma)
    assert(any(abs(sigmaNew - sigma(k)) < 1e-15), 'original node %d lost', k);
end
% inserts are the midpoints of intervals 3 and 7
mids = sort([(sigma(3)+sigma(4))/2; (sigma(7)+sigma(8))/2]);
got  = sort(sigmaNew(isNew));
assert(max(abs(got - mids)) < 1e-15, 'inserts must be interval midpoints');

% guard: hFloor skips a hot but too-thin interval
sigma2 = [0; 1e-10; 0.5; 1];  score2 = [9 0 0];   % interval 1 is 1e-10 wide
opts2  = struct('K', 1, 'hFloor', 1e-9, 'maxAdd', 40);
[sn2, in2, nd2] = refine_sigma(sigma2, score2, opts2);
assert(nnz(in2) == 0 && nd2 == 1, 'sub-hFloor interval must be dropped');

fprintf('ALL PASS\n');
end
