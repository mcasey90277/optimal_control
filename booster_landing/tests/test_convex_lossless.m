% TEST_CONVEX_LOSSLESS  Fixed-tf convexified PDG at a plausible tf: solved
% status, losslessness gap ||u||-sigma ~ 0 (the relaxation is TIGHT at the
% optimum -- the Acikmese/Ploen theorem, checked numerically), annulus in
% original variables, terminal conditions met (pad, descending at P.vf --
% ADJUDICATED 2026-08-08, was v(tf)=0, see booster_params.m).
%
% Uses solve_pdg_convex's single-tf mode (opts.tf fixed) to stay fast.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P   = booster_params();
sol = solve_pdg_convex(P, struct('tf', 25, 'Nconv', 60));
assert(sol.stats.success, 'convex solve failed: %s', sol.stats.status);
assert(sol.lossless_gap < 1e-4 * P.Tmax/P.m0, ...
       'relaxation not tight: gap=%.3e', sol.lossless_gap);
Tmag = sqrt(sum(sol.U.^2,1));
assert(all(Tmag >= P.Tmin*(1-1e-4)) && all(Tmag <= P.Tmax*(1+1e-4)), ...
       'annulus violated in original variables');
assert(max(abs(sol.X(1:6,end) - [0;0;0;P.vf])) < 1e-2, 'terminal not met');
fprintf('test_convex_lossless PASS  gap=%.2e  mf=%.1f kg\n', ...
        sol.lossless_gap, sol.mf);
