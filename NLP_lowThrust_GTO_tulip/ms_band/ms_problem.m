function prob = ms_problem(factor, epsSmooth)
% MS_PROBLEM  Problem struct factory for the ms_band multiple-shooting solver.
%
% Centralizes constants, endpoints, and solver settings for the min-fuel
% GTO -> tulip transfer at t_f = factor x min-time. Arc joints prob.tJ are
% left empty; seed builders fill them (uniform in Sundman tau).
%
% INPUTS:
%   factor    - t_f as a multiple of the campaign min-time 6.290694 ND [scalar]
%   epsSmooth - Bertrand-Epenoy throttle smoothing parameter [scalar]
%
% OUTPUTS:
%   prob - struct: tfMin, factor, tf, rv0 [6x1], m0, rvf [6x1], Tmax, c,
%          muStar, p (full cr3bp_lt_params struct), epsSmooth,
%          tJ [] (joint times, 1x(M+1), set by caller), odeOpts

p          = cr3bp_lt_params(0.025, 15, 2100);
[rv0, rvf] = gto_tulip_endpoints(p);

prob.tfMin     = 6.290694;              % campaign constant (ND)
prob.factor    = factor;
prob.tf        = factor * prob.tfMin;
prob.rv0       = rv0(:);
prob.m0        = 1;
prob.rvf       = rvf(:);
prob.Tmax      = p.Tmax;
prob.c         = p.c;
prob.muStar    = p.muStar;
prob.p         = p;
prob.epsSmooth = epsSmooth;
prob.tJ        = [];
% RelTol 1e-13: costate rows (~2e2 magnitude through perigee arcs) set the
% MS residual floor at ~(costate scale x RelTol x sqrt(Nrows)) — measured
% 1.206e-9 at RelTol 1e-12, 2.253e-10 at 1e-13 (Task-5 floor analysis,
% review-verified). 1e-13 keeps the certified ||R|| <= 1e-9 gate honest.
prob.odeOpts   = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
end
