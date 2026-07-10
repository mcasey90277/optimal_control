function [Zseed, tJ] = seed_from_mintime(factor, M, epsSeed)
% SEED_FROM_MINTIME  MS seed for a near-min-time factor from min-time costates.
%
% Integrates the smoothed min-fuel PMP dynamics from the converged min-time
% initial costates over the LONGER horizon [0, factor*tfMin] and samples at
% tau-uniform joints. The end state misses the target (that is what the MS
% solve closes); at eps = 1 the basin is wide.
%
% INPUTS:
%   factor  - target tf factor (e.g. 1.01) [scalar]
%   M       - number of arcs [scalar]
%   epsSeed - smoothing value the seed is built at (use epsSchedule(1)=1)
%
% OUTPUTS:
%   Zseed - MS unknown seed [(14M-7)x1]
%   tJ    - tau-uniform joint times [1x(M+1)]

ref  = run_gto_tulip_indirect(false);
lam0 = ref.zSol(1:7);
prob = ms_problem(factor, epsSeed);
y0   = [prob.rv0; prob.m0; lam0];
sol  = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsSeed), [0 prob.tf], y0, prob.odeOpts);
tJ   = arc_boundaries_tau(sol.x, sol.y(1:3, :), M, prob.muStar);
yJ   = deval(sol, tJ);
Zseed = ms_pack(lam0, yJ(:, 2:M));
end
