function best = eps_march(Zseed, probBase, epsSchedule, tolR)
% EPS_MARCH  Guarded Bertrand-Epenoy smoothing continuation at fixed tf.
%
% Solves the MS system down an eps schedule, warm-starting each step from
% the last converged solution. A failed step triggers geometric bisection
% between the last success and the failed target; a step that still fails
% is abandoned (guard: the warm start is never poisoned).
%
% Within-step relay (campaign amendment, 2026-07-10): if a step's solve
% returns non-success with flag==0 (iteration cap) AND its ||R|| fell by
% >=10% relative to that step's starting residual, the step is re-solved
% at the SAME eps from its own returned iterate (never from a different
% eps) — up to 4 relays. Relaying stops as soon as a relay converges, or
% as soon as a relay improves <10%; the step is then failed and its
% (non-converged) iterate is discarded per the normal guard, exactly as
% a step that never relayed at all.
%
% success requires reaching eps <= 1e-3 (target 1e-4).
%
% INPUTS:
%   Zseed       - seed at epsSchedule(1) [(14M-7)x1]
%   probBase    - problem struct with tJ set (epsSmooth overridden per step)
%   epsSchedule - decreasing values [1xQ], default
%                 [1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4]
%   tolR        - per-solve success threshold [scalar, default 1e-9]
%
% OUTPUTS:
%   best - struct: Z, eps (smallest converged), resNorm, success
%          (eps <= 1e-3), history (struct array per attempted step, incl.
%          relays actually used on that step)

if nargin < 3 || isempty(epsSchedule)
    epsSchedule = [1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4];
end
if nargin < 4 || isempty(tolR), tolR = 1e-9; end

MAX_RELAY   = 4;
MAX_LM_ITER = 200;

list = epsSchedule(:).';
best = struct('Z', [], 'eps', Inf, 'resNorm', Inf, 'success', false, ...
              'history', struct('eps', {}, 'resNorm', {}, 'converged', {}, ...
                                 'relays', {}));
Zwarm = Zseed(:);
idx = 1;
while idx <= numel(list)
    epsK = list(idx);
    prob = probBase;  prob.epsSmooth = epsK;

    Rstart = norm(ms_residual(Zwarm, prob));         % step's starting residual
    out    = ms_solve(Zwarm, prob, tolR, MAX_LM_ITER);

    nRelay = 0;
    while ~out.success && out.flag == 0 && out.resNorm <= 0.9*Rstart ...
            && nRelay < MAX_RELAY
        nRelay = nRelay + 1;
        fprintf('eps_march: relay %d/%d at eps=%.3g (||R||=%.3e, >=10%% cut vs %.3e) — resolving from own iterate\n', ...
                nRelay, MAX_RELAY, epsK, out.resNorm, Rstart);
        Rstart = out.resNorm;
        out    = ms_solve(out.Z, prob, tolR, MAX_LM_ITER);
    end

    best.history(end+1) = struct('eps', epsK, 'resNorm', out.resNorm, ...
                                 'converged', out.success, 'relays', nRelay); %#ok<AGROW>

    if out.success
        Zwarm = out.Z;
        best.Z = out.Z;  best.eps = epsK;  best.resNorm = out.resNorm;
        idx = idx + 1;
    elseif isfinite(best.eps) && sqrt(best.eps*epsK) < 0.9*best.eps
        list = [list(1:idx-1), sqrt(best.eps*epsK), list(idx:end)];  % bisect
    else
        fprintf('eps_march: abandoning at eps=%.3g (relays used=%d, ||R||=%.3e)\n', ...
                epsK, nRelay, out.resNorm);
        break;
    end
end
best.success = best.eps <= 1e-3;
end
