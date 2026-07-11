function test_prep_refine_seed()
% TEST_PREP_REFINE_SEED  Prepared 1.15x seed carries duals + required fields.
%
% NOTE: runs one eps=0 re-solve (~1-3 min).
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   % ms_band/setup_paths does NOT add ms_band itself
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();  cd(old);

src = fullfile(here, '..', 'sundman_minfuel_certified.mat');
out = fullfile(tempdir, 'refine_seed_test.mat');
outFile = prep_refine_seed(src, out);

S = load(outFile);
assert(isfield(S, 'out') && isfield(S.out, 'lamDef') && ~isempty(S.out.lamDef), 'lamDef present');
assert(isequal(size(S.out.lamDef, 1), 8), 'lamDef is 8xN');
assert(isfield(S, 'factor') && abs(S.factor - 1.15) < 1e-9, 'factor = 1.15');
assert(all(isfield(S, {'tauf0', 'sigma', 'rv0', 'rvf'})), 'required fields present');
assert(S.out.maxDefect < 1e-6, 'seed re-solve converged tight, got %.2e', S.out.maxDefect);
fprintf('ALL PASS (switches=%d, defect=%.2e)\n', S.out.switches, S.out.maxDefect);
end
