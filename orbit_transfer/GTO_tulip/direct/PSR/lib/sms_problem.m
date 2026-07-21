function prob = sms_problem(factor, epsSmooth)
% SMS_PROBLEM  Problem struct factory for the Sundman-domain MS solver.
%
% Sundman-domain clone of MS_PROBLEM: same constants/endpoints, plus the
% Sundman exponent pSund and the sigma-domain joint fields sigf / sJ
% (filled by seed builders; sigf is FIXED, never a decision variable).
% prob.resFun routes MS_SOLVE / EPS_MARCH to the 16-dim residual.
%
% INPUTS:
%   factor    - t_f as a multiple of the campaign min-time 6.290694 ND [scalar]
%   epsSmooth - Bertrand-Epenoy throttle smoothing parameter [scalar]
%
% OUTPUTS:
%   prob - struct: tfMin, factor, tf, rv0 [6x1], m0, rvf [6x1], Tmax, c,
%          muStar, p (full cr3bp_lt_params struct), epsSmooth, pSund,
%          sigf [] (total sigma length, set by seed builders), sJ []
%          (sigma joints, 1x(M+1), set by seed builders), resFun
%          (@sms_residual), odeOpts

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
prob.pSund     = 1.5;
prob.sigf      = [];
prob.sJ        = [];
prob.resFun    = @sms_residual;
% RelTol 1e-13 / AbsTol 1e-15: campaign-amended tolerances (Task-5 floor
% analysis) keeping the certified ||R|| <= 1e-9 gate honest.
prob.odeOpts   = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
end
