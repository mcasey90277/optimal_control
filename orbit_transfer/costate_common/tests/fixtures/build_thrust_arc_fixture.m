function build_thrust_arc_fixture()
% BUILD_THRUST_ARC_FIXTURE  One-off: distil the 2026-09-08 thrust-arclength
% diagnostic run (arclength_thrust, normal chart, 75.5 -> 72.0 mN, 401
% roots) into a small fixture for test_arclength_ms_thrust: the first root
% as a full ms unknown vector, the problem constants, and the reference
% (T, t_f) curve. The unknown packing is ms_bvp's:
%   p = [lam0(7); Y_2 .. Y_K (14 each); t_f],   Y(:,1) = [rv0; 1; lam0].
%
% INPUTS:  none (reads DRO_tulip/indirect/results/mintime_arclength_diag.mat
%          and mintime_70mN_direct.mat)
% OUTPUTS: none (writes fixtures/thrust_arc_75mN.mat, force-added to git)
here = fileparts(mfilename('fullpath'));
res = fullfile(here, '..', '..', '..', 'DRO_tulip', 'indirect', 'results');
S = load(fullfile(res, 'mintime_arclength_diag.mat'));  A = S.Aout;
D = load(fullfile(res, 'mintime_70mN_direct.mat'));
F = struct();
F.p0 = A.p{1};  F.TN0 = A.T_N(1);
F.rv0 = D.rv0(1:6);  F.rvf = D.rvf(1:6);  F.cnd = D.cnd;  F.muStar = D.muStar;
F.m0kg = 150;
F.refT_N = A.T_N;  F.reftf_nd = A.tf_nd;  F.reflam0 = A.lam0;
F.note = 'arclength_thrust diagnostic arc, normal chart, 2026-09-08 (Dx = max(|p0|,1e-2), sT = TN0, ds 0.05, dsMax 3.0)';
save(fullfile(here, 'thrust_arc_75mN.mat'), '-struct', 'F');
fprintf('fixture written: n = %d, K = %d, %d reference roots\n', numel(F.p0), (numel(F.p0)-8)/14+1, numel(F.refT_N));
end
