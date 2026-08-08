% TEST_COLLOC_SMOKE  Coarse-grid (N=15) nonconvex NLP solve converges and
% obeys physics: solved status, mass in (mdry, m0), thrust annulus and
% glideslope satisfied at nodes, terminal state at the pad at rest.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 15));
assert(sol.stats.success, 'IPOPT did not converge: %s', sol.stats.status);
assert(sol.mf > P.mdry && sol.mf < P.m0, 'final mass out of range');
Tmag = sqrt(sum(sol.U.^2, 1));
assert(all(Tmag >= P.Tmin - 1) && all(Tmag <= P.Tmax + 1), 'annulus violated');
rxy  = sqrt(sum(sol.X(1:2,:).^2, 1));
assert(all(rxy <= sol.X(3,:)/tand(P.gs_deg) + 1e-3), 'glideslope violated');
assert(max(abs(sol.X(1:6,end))) < 1e-3, 'terminal state not at rest on pad');
fprintf('test_colloc_smoke PASS  tf=%.2f s  mf=%.1f kg  fuel=%.1f kg\n', ...
        sol.tf, sol.mf, P.m0 - sol.mf);
