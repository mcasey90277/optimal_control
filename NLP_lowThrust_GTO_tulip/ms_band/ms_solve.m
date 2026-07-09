function out = ms_solve(Zseed, prob, tolR, maxIter)
% MS_SOLVE  Guarded Levenberg-Marquardt solve of the MS system at fixed (tf, eps).
%
% Wraps lsqnonlin (levenberg-marquardt, analytic sparse Jacobian from
% MS_RESIDUAL). LM measured 4x better than dogleg on this problem family
% (campaign record). The input seed is never modified; on failure the
% caller keeps its own warm start (guard discipline).
%
% INPUTS:
%   Zseed   - unknown-vector seed [(14M-7)x1]
%   prob    - problem struct with tJ set [1x(M+1)]
%   tolR    - success threshold on ||R||_2 [scalar, e.g. 1e-9]
%   maxIter - LM iteration cap [scalar]
%
% OUTPUTS:
%   out - struct: Z [(14M-7)x1], resNorm, flag (lsqnonlin exitflag),
%         success (resNorm <= tolR), iterations

opts = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', 20*maxIter);

[Z, res2, ~, flag, outp] = lsqnonlin(@(zz) ms_residual(zz, prob), ...
                                     Zseed(:), [], [], opts);
out = struct('Z', Z, 'resNorm', sqrt(res2), 'flag', flag, ...
             'success', sqrt(res2) <= tolR, 'iterations', outp.iterations);
fprintf('ms_solve: eps=%.3g tf=%.6f ||R||=%.3e iters=%d flag=%d\n', ...
        prob.epsSmooth, prob.tf, out.resNorm, out.iterations, flag);
end
