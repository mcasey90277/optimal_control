function [Zseed, prob, sol] = sms_seed_mintime(factor, M, epsSeed)
% SMS_SEED_MINTIME  Sundman-MS seed from the converged min-time costates.
%
% Builds the 16-dim initial state [rv0; 1; 0; lam0(7); lamT0] with
% lamT0 = -Ht(0) (so H_sigma = kappa*(Ht+lamT) = 0 at sigma = 0), then
% integrates SMS_EOM from sigma = 0 with an ode113 event on the carried
% time state Y(8) reaching tf: that event sigma IS sigf (fixed thereafter).
% Joints are sigma-uniform, samples come from the same run (deval). The
% end state misses the target for factor > 1 (that is what MS closes).
%
% INPUTS:
%   factor  - target tf factor (1.00 uses the converged min-time tf
%             ref.zSol(8) exactly, matching test_ms_reproduce_mintime)
%   M       - number of arcs [scalar]
%   epsSeed - smoothing value the seed is built at [scalar]
%
% OUTPUTS:
%   Zseed - MS unknown seed [(16M-8)x1]
%   prob  - problem struct from SMS_PROBLEM with tf, sigf, sJ set
%   sol   - ode113 solution struct of the seed integration (diagnostics)

ref  = run_gto_tulip_indirect(false);
lam7 = ref.zSol(1:7);
prob = sms_problem(factor, epsSeed);
if abs(factor - 1) < 1e-9
    prob.tf = ref.zSol(8);          % converged min-time tf exactly
end

% sigf estimate from a time-domain reference run (sets the event span)
solT = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsSeed), [0 prob.tf], ...
              [prob.rv0; prob.m0; lam7], prob.odeOpts);
r1     = sqrt(sum((solT.y(1:3, :) - [-prob.muStar; 0; 0]).^2, 1));
sigEst = trapz(solT.x, 1./r1.^prob.pSund);

% lamT0 = -Ht(0): Ht from SMS_EOM (independent of lamT)
y0 = [prob.rv0; prob.m0; 0; lam7; 0];
[~, Ht0] = sms_eom(0, y0, prob.Tmax, prob.c, prob.muStar, epsSeed, prob.pSund);
lam0   = [lam7; -Ht0];
y0(16) = -Ht0;

opts = odeset(prob.odeOpts, 'Events', @(s, y) tf_event(y, prob.tf));
sol  = ode113(@(s, y) sms_eom(s, y, prob.Tmax, prob.c, prob.muStar, ...
              epsSeed, prob.pSund), [0 1.5*sigEst], y0, opts);
if isempty(sol.xe)
    error('sms_seed_mintime:event', 't never reached tf within 1.5x sigma estimate');
end
prob.sigf = sol.xe(end);
prob.sJ   = linspace(0, prob.sigf, M+1);

yJ    = deval(sol, prob.sJ);
Zseed = sms_pack(lam0, yJ(:, 2:M));
end

% -------------------------------------------------------------------------
function [val, isterm, dir] = tf_event(y, tf)
% TF_EVENT  Terminal event: carried time state reaches tf (upward crossing).
%
% INPUTS:
%   y  - augmented state [16x1] (real; seed integrations only)
%   tf - target transfer time (ND) [scalar]
%
% OUTPUTS:
%   val    - event value t - tf [scalar]
%   isterm - 1 (stop integration)
%   dir    - +1 (upward crossings only; t is strictly increasing)
val    = y(8) - tf;
isterm = 1;
dir    = 1;
end
