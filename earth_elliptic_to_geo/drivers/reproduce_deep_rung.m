function best = reproduce_deep_rung(thrustN, warmMatPath, cfg)
% REPRODUCE_DEEP_RUNG  Certify a DEEP low-thrust rung (0.2 N, 0.1 N, ...) via the
% external-review recipe that first cracked 0.2 N (2026-07-20). Warm-chains from
% the certified rung above and applies the four levers that broke the deep-rung
% wall. This is the reproducible driver for the result recorded in
% process/DEEP_THRUST_LESSONS.md.
%
% THE RECIPE (each lever fixes a specific wall -- see the lessons note):
%   1. rung-adaptive dL bound   -> feasibility (a fixed dL<=2000 made 0.2/0.1 N
%                                  structurally infeasible: DeltaL~2168/4335)
%   2. liftDL = true            -> block-banded KKT (the scalar DeltaL is a dense
%                                  Jacobian column -> arrowhead KKT -> MUMPS/METIS
%                                  crash at n~30k+); liftDL caps col-nnz at ~O(1)
%   3. phase-correct beta       -> correct initial switching structure (a plain
%                                  sigma-interp warm-start beta phase-ALIASES
%                                  across rungs; warmstart_phase_beta fixes it)
%   4. eps-continuation 1->0    -> gentle energy->fuel path (NOT a direct eps=0
%      with GENEROUS maxIter        bang-bang solve) AND enough iterations so
%                                  every step converges tight (ok=1, defect~1e-13)
%                                  -- UNDER-ITERATION collapses the deep-eps tail.
%   NOT USED: scaleNLP (fights IPOPT gradient-based auto-scaling -> eps=1
%             restoration failure; proper explicit scaling needs a COMPLETE
%             user-scaling pass, still open). adaptiveEps bisection is armed as a
%             safety net (helps only if a step fails; with enough maxIter it did not).
%
% INPUTS:
%   thrustN     - target thrust [N] (e.g. 0.2, 0.1)
%   warmMatPath - path to the certified warm-source rung .mat (the rung ABOVE),
%                 a res-struct file (fields res.fuel.{X,U,dL}, res.sigma). e.g.
%                 'results/MEE_M2_0p5N.mat' for 0.2 N, 'results/MEE_M2_0p2N.mat'
%                 for 0.1 N.
%   cfg         - struct (optional): .nodesPerRev [8], .maxIter [3000],
%                 .ctf [1.5], .R0 [223.14] (R0-law tfmin ~ R0/thrustN [ND]),
%                 .m0kg [1500], .ispS [2000], .resDir [tempdir], .printLevel [5]
%
% OUTPUTS:
%   best - homotopy_mee best-struct: .certified .epsReached .m_f_kg .switches
%          .maxDefect .termErr .incDeg .dL .X .U (a certified fuel solution iff
%          best.certified == 1, i.e. epsReached==0 with maxDefect < 1e-8)
%
% REFERENCES:
%   [1] process/DEEP_THRUST_LESSONS.md (this recipe + the lessons that produced it).
%   [2] external review (GPT-5.6-terra + Gemini 3.1 Pro), 2026-07-19.
%   [3] warmstart_phase_beta.m, homotopy_mee.m (adaptiveEps), casadi_lt_mee.m (liftDL).
if nargin < 3, cfg = struct(); end
d = @(f,v) optdef(cfg, f, v);
nodesPerRev = d('nodesPerRev', 8);
maxIter     = d('maxIter', 3000);
ctf         = d('ctf', 1.5);
R0          = d('R0', 223.14);
m0kg        = d('m0kg', 1500);
ispS        = d('ispS', 2000);
resDir      = d('resDir', fullfile(tempdir, sprintf('deeprung_%gN', thrustN)));
printLevel  = d('printLevel', 5);

S = load(warmMatPath);
assert(isfield(S,'res'), 'reproduce_deep_rung: %s must hold a res-struct', warmMatPath);
src = S.res;  Tprev = src.fp.thrustN;
assert(Tprev > thrustN, 'warm source thrust (%g) must be ABOVE the target (%g)', Tprev, thrustN);

% C-law: DeltaL_guess = DeltaL_prev * (T_prev / T_new); N from the implied revs.
dL0 = src.fuel.dL * (Tprev / thrustN);
N   = round(nodesPerRev * dL0/(2*pi));
sigmaDst = linspace(0, 1, N+1).';
par = kepler_lt_params(thrustN, m0kg, ispS);
tf  = ctf * (R0 / thrustN);                          % anchor-free R0-law tfmin * ctf

% warm start: interp state from the rung above, rescale the time row to this tf,
% and RECOMPUTE beta phase-correctly (kills the sigma-interp aliasing).
W  = interp_warmstart(src.fuel.X, src.fuel.U, src.fuel.dL, src.sigma, sigmaDst);
X0 = W.X;
if X0(7,end) > 0, X0(7,:) = X0(7,:) * (tf / X0(7,end)); end
U0 = warmstart_phase_beta(X0, sigmaDst, dL0, par, W.U(4,:));

fprintf(['reproduce_deep_rung: %g N warm-chained from %g N | N=%d (n~%d), dL0=%.1f, ' ...
         'tf=%.1f | liftDL + phase-beta + eps-cont(maxIter=%d)\n'], ...
        thrustN, Tprev, N, 11*(N+1)+1, dL0, tf, maxIter);

sched = [1 0.7 0.5 0.35 0.25 0.18 0.12 0.08 0.05 0.03 0.02 0.012 0.007 0.004 0.002 0.001 0];
ho = struct('par', par, 'x0', X0(:,1), 'tfTarget', tf, 'maxIter', maxIter, ...
    'resDir', resDir, 'tag', sprintf('deep_%gN', thrustN), 'printLevel', printLevel, ...
    'fp', struct('thrustN', thrustN), 'xf', [1;0;0;0;0], 'sched', sched, ...
    'liftDL', true, 'scaleNLP', false, 'adaptiveEps', true, 'epsMinStep', 3e-4, ...
    'maxSolves', 60);
[best, ~] = homotopy_mee(sigmaDst, X0, U0, dL0, ho);

fprintf(['reproduce_deep_rung: %g N -> certified=%d epsReached=%.4g m_f=%.4f kg ' ...
         'switches=%d defect=%.2e termErr=%.2e revs=%.1f\n'], ...
        thrustN, best.certified, best.epsReached, best.m_f_kg, best.switches, ...
        best.maxDefect, best.termErr, best.dL/(2*pi));
end
