% TEST_TVLQR_RICCATI  Riccati solution health along the nominal trajectory:
% P(t) symmetric positive semidefinite everywhere (PSD, not PD: the mass
% row/column is unweighted), gains finite, terminal condition P(tf)=Qf.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
assert(sol.stats.success);
ctrl = tvlqr_design(sol, P);
M    = numel(ctrl.tgrid);
for k = 1:M
    Pk = ctrl.Pt(:,:,k);
    assert(max(max(abs(Pk - Pk.'))) < 1e-6 * max(1,max(abs(Pk(:)))), ...
           'P not symmetric at k=%d', k);
    assert(min(eig((Pk+Pk.')/2)) > -1e-8 * max(1,max(abs(Pk(:)))), ...
           'P not PSD at k=%d', k);
end
assert(all(isfinite(ctrl.K(:))), 'non-finite gains');
% Compare against the weight tvlqr_design actually used (ctrl.Qf) rather
% than a hardcoded literal: task-7b made Qf(6,6) a DERIVED quantity (the
% phase-B vertical velocity weight, so the Riccati boundary layer does not
% collapse the terminal vertical gain -- see tvlqr_design). The property
% under test is "the terminal condition P(tf)=Qf is honored", which is what
% this now checks, independent of what the default weights happen to be.
assert(isfield(ctrl,'Qf'), 'tvlqr_design no longer reports the Qf it used');
assert(max(max(abs(ctrl.Pt(:,:,end) - ctrl.Qf))) < 1e-9*max(1,max(abs(ctrl.Qf(:)))), ...
       'terminal P ~= Qf');
fprintf('test_tvlqr_riccati PASS\n');
