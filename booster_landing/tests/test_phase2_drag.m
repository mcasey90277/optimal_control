% TEST_PHASE2_DRAG  Drag-on collocation solve (coarse, warm-started from
% vacuum) converges; fuel differs from vacuum by a NONZERO but sane amount
% (drag helps braking: expect LESS fuel, sanity band 0 < dfuel < 40% of
% vacuum fuel); certify gates G1/G2/G5 pass with G3/G4 skipped. Also
% smoke-tests viz/plot_vacuum_vs_drag.m on the same solC/solD (task-11
% close-out review, Important-5: this campaign's fast suite otherwise had
% zero coverage of the Phase-2 branch's own viz product).
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
Pd   = P;  Pd.drag.on = true;
solD = solve_pdg_colloc(Pd, struct('N', 30, 'init', solC));
assert(solD.stats.success, 'drag solve failed');
fuelV = P.m0 - solC.mf;  fuelD = P.m0 - solD.mf;
assert(fuelD < fuelV, 'drag should reduce fuel (braking), got +%.1f kg', ...
       fuelD - fuelV);
assert(fuelV - fuelD < 0.4*fuelV, 'drag effect implausibly large');
rep = certify_pdg(solD, [], Pd);
assert(rep.G1_pass && rep.G2_pass && rep.G5_pass, 'drag gates failed');
assert(isequal(rep.G3_pass, 'skipped'), 'G3 should be skipped without twin');

% Fast coverage for viz/plot_vacuum_vs_drag.m (task-11 close-out review,
% Important-5: the new Phase-2 branch/viz had no fast test at all before
% this). Reuses the solC/solD already solved above -- no extra solve.
vizfile = fullfile(tempdir, 'test_phase2_vac_vs_drag.png');
if isfile(vizfile), delete(vizfile); end
fig = plot_vacuum_vs_drag(solC, solD, vizfile);
close(fig);
assert(isfile(vizfile), 'plot_vacuum_vs_drag did not write its PNG');

fprintf('test_phase2_drag PASS  fuel vac=%.1f drag=%.1f kg\n', fuelV, fuelD);
