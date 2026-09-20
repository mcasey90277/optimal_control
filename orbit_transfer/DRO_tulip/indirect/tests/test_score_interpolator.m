function ok = test_score_interpolator()
% TEST_SCORE_INTERPOLATOR  The interpolator's scorer, with the SOLVER REPLACED by
% a fake (opts.solver) so that the bookkeeping is tested in seconds: the setup
% it shares with every consumer (catalog_blend_setup, checked against the
% sensitivities computed and stored on 09-19 by another code path), the
% queries, resume, how a solve that does not finish is counted, and the
% summary's arithmetic.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
recDir = fullfile(here, 'results', 'library_70mN_24x24_final');
catMat = fullfile(recDir, 'costate_catalog_dro_tulip_70mN.mat');
assert(isfile(catMat), 'the library of record is missing: %s', catMat);

% ---- the shared setup ---------------------------------------------------------
K = catalog_blend_setup(catMat);
P = load(fullfile(recDir, 'phase_transversality.mat'));  fp = fieldnames(P);  X = P.(fp{1});
both = isfinite(K.GD) & isfinite(X.dTf_dsD);
% (arrival: bit for bit. Departure: the stored map used a 1e-6 difference quotient whose round-off error is 5.8e-4;
%  this one uses phase_state's analytic derivative -- see phase_transversality_check/departureTangent)
ok = chk(ok, nnz(both) == 576 && max(abs(K.GD(both) - X.dTf_dsD(both))) < 1e-3 && max(abs(K.GA(both) - X.dTf_dsA(both))) < 1e-9, ...
         sprintf('the sensitivity maps agree with the ones phase_transversality_check stored (arrival %.1e, departure %.1e)', ...
                 max(abs(K.GA(both) - X.dTf_dsA(both))), max(abs(K.GD(both) - X.dTf_dsD(both)))));
ok = chk(ok, abs(K.gapD - 1/24) < 1e-12 && abs(K.gapA - 1/24) < 1e-12 && max(abs(K.edgeCap*K.phys.tStar/60 - 300)) < 1e-9, 'grid steps 1/24 and the 300 min cap that goes with them');
[zG, blend] = catalog_blend_guess(K, 0.77, 0.26);
[zL, ~] = phase_catalog_interp(K.sheet, 0.77, 0.26, struct('dTdsD', K.GD, 'dTdsA', K.GA, 'maxEdgeResidual', K.edgeCap));
ok = chk(ok, strcmp(blend.tier, 'bilinear') && isequal(zG(1:7), zL(1:7)) && abs(zG(8) - zL(8))*K.phys.tStar/60 > 60, ...
         'catalog_blend_guess: the interpolator''s costates, with t_f from the corners'' sensitivities (hours from the plain blend here)');
[zP, ~] = catalog_blend_guess(K, 0.77, 0.26, struct('tfFromSensitivities', false));
ok = chk(ok, isequal(zP, zL), 'and the plain blend when asked');

% ---- the scorer, solver faked: "converges" exactly when the blend is bilinear ---------
outDir = fullfile(tempdir, sprintf('score_interp_test_%d', feature('getpid')));  if isfolder(outDir), rmdir(outDir, 's'); end
calls = 0;
    function [finished, z, info] = fakeSolver(zGuess, ~, ~, ~)
        calls = calls + 1;
        finished = true;  z = NaN(8, 1);  info = struct('converged', false, 'normR', 1, 'iters', 600, 'wall', 1, 'guessMissKm', 5e4);
        if abs(zGuess(8) - round(zGuess(8), 3)) < 0.25e-3          % an arbitrary, deterministic half of the guesses "converge"
            z = zGuess;  info.converged = true;  info.iters = 10;  info.guessMissKm = 100;
        end
    end
o = struct('nQuery', 12, 'seed', 3, 'outDir', outDir, 'solver', @fakeSolver, 'pool', [], 'print', false);
R = score_interpolator(K, o);
ok = chk(ok, numel(R.rows) == 12 && all(isfield(R.rows, {'sD', 'sA', 'tier', 'converged', 'sameBranch', 'usable', 'guessMissKm', 'iters', 'base'})), 'one row per query, with the blend''s and the baseline''s outcome');
nSolved = calls;
ok = chk(ok, nSolved == 2*nnz(~strcmp({R.rows.tier}, 'none')), sprintf('two solves per query that has a guess: the blend and the baseline (%d)', nSolved));
ok = chk(ok, R.summary.nUsable == nnz([R.rows.usable]) && abs(R.summary.hitRate - R.summary.nUsable/12) < 1e-15 && ...
             sum([R.summary.byTier.n]) == 12, 'the summary is the rows'' arithmetic, tier by tier');
ok = chk(ok, all(~[R.rows(~[R.rows.converged]).usable]) && all(isnan([R.rows(~[R.rows.converged]).q2])), 'a query that did not converge is never usable and has no root to compare');
ok = chk(ok, isfile(fullfile(outDir, 'SCORE_VERDICT.txt')) && isfile(fullfile(outDir, 'score_rows.mat')), 'a verdict FILE and the rows are written');
R2 = score_interpolator(K, o);
ok = chk(ok, calls == nSolved && isequal([R2.rows.sD], [R.rows.sD]), 'run again: nothing is solved twice (resume), the same queries (seeded)');
o.nQuery = 15;  R3 = score_interpolator(K, o);
ok = chk(ok, numel(R3.rows) == 15 && isequal([R3.rows(1:12).sD], [R.rows.sD]) && calls > nSolved, 'asked for more: the first twelve are kept, only the new three are solved');
% ---- a solve the fence stopped is a FAILURE, not a skipped row --------------------------
rmdir(outDir, 's');
o = struct('nQuery', 4, 'seed', 3, 'outDir', outDir, 'pool', [], 'print', false, 'doBaseline', false, ...
           'solver', @(varargin) deal(false, [], []));
R = score_interpolator(K, o);
ok = chk(ok, numel(R.rows) == 4 && ~any([R.rows.usable]) && all(strcmp({R.rows(~strcmp({R.rows.tier}, 'none')).outcome}, 'stopped by the cap')), 'a solve the fence stopped counts against the hit rate, by name');
rmdir(outDir, 's');
if ok, fprintf('test_score_interpolator: ALL PASS\n'); else, fprintf('test_score_interpolator: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
