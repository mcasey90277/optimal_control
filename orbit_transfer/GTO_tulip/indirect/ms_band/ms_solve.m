function out = ms_solve(Zseed, prob, tolR, maxIter)
% MS_SOLVE  Guarded Levenberg-Marquardt solve of the MS system at fixed (tf, eps).
%
% Wraps lsqnonlin (levenberg-marquardt, analytic sparse Jacobian from
% MS_RESIDUAL). LM measured 4x better than dogleg on this problem family
% (campaign record). ScaleProblem='jacobian' (Marquardt scaling) added
% Task 7: probe-verified ~10x early-phase speedup absorbing the costate
% (~2e2) vs state (~1) row-magnitude mismatch (probe_scaled_lm.m/.log).
% The input seed is never modified; on failure the caller keeps its own
% warm start (guard discipline). Task S1: optional prob.resFun selects the
% residual (default @ms_residual; the Sundman-domain path sets
% @sms_residual via SMS_PROBLEM).
%
% INPUTS:
%   Zseed   - unknown-vector seed [(14M-7)x1] (16M-8 on the sms path)
%   prob    - problem struct with tJ set [1x(M+1)] (sJ on the sms path)
%   tolR    - success threshold on ||R||_2 [scalar, e.g. 1e-9]
%   maxIter - LM iteration cap [scalar]
%
% OUTPUTS:
%   out - struct: Z [(14M-7)x1], resNorm, flag (lsqnonlin exitflag),
%         success (resNorm <= tolR), iterations

opts = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'levenberg-marquardt', ...
    'ScaleProblem', 'jacobian', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', 20*maxIter);

if isfield(prob, 'resFun') && ~isempty(prob.resFun)
    resFun = prob.resFun;
else
    resFun = @ms_residual;
end
[Z, res2, ~, flag, outp] = lsqnonlin(@(zz) resFun(zz, prob), ...
                                     Zseed(:), [], [], opts);
out = struct('Z', Z, 'resNorm', sqrt(res2), 'flag', flag, ...
             'success', sqrt(res2) <= tolR, 'iterations', outp.iterations);
fprintf('ms_solve: eps=%.3g tf=%.6f ||R||=%.3e iters=%d flag=%d\n', ...
        prob.epsSmooth, prob.tf, out.resNorm, out.iterations, flag);
end
