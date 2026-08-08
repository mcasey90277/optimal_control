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
Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);
assert(max(max(abs(ctrl.Pt(:,:,end) - Qf))) < 1e-9, 'terminal P ~= Qf');
fprintf('test_tvlqr_riccati PASS\n');
