function out = ifs_solve(Z0, prob, opts)
% IFS_SOLVE  Levenberg-Marquardt solve of the IFS multiple-shooting system.
%
% Wraps lsqnonlin (levenberg-marquardt, our sparse complex-step Jacobian from
% IFS_RESIDUAL, ScaleProblem='jacobian' to absorb the costate-vs-state row
% magnitude mismatch). The seed is never modified; the caller keeps its warm
% start on failure.
%
% INPUTS:
%   Z0   - seed unknown vector [(8+17k)x1]
%   prob - problem struct (see plan Shared data layout)
%   opts - struct: tolR success threshold on ||R||_2 [1e-8], maxIter [200]
% OUTPUTS:
%   out - struct: Z, resNorm, iterations, flag, success (resNorm<=tolR),
%         seedResNorm
%
% REFERENCES:
%   [1] ms_band/ms_solve.m (LM settings). [2] the IFS design spec.
if ~isfield(opts,'tolR'),    opts.tolR = 1e-8;  end
if ~isfield(opts,'maxIter'), opts.maxIter = 200; end
lmopts = optimoptions('lsqnonlin', ...
    'Display','iter', 'Algorithm','levenberg-marquardt', ...
    'ScaleProblem','jacobian', 'SpecifyObjectiveGradient',true, ...
    'FunctionTolerance',1e-24, 'StepTolerance',1e-14, ...
    'MaxIterations',opts.maxIter, 'MaxFunctionEvaluations',50*opts.maxIter);
seedResNorm = norm(ifs_residual(Z0, prob));
[Z, res2, ~, flag, op] = lsqnonlin(@(zz) ifs_residual(zz, prob), Z0(:), [], [], lmopts);
out = struct('Z',Z,'resNorm',sqrt(res2),'iterations',op.iterations,'flag',flag, ...
             'success', sqrt(res2) <= opts.tolR, 'seedResNorm', seedResNorm);
fprintf('ifs_solve: k=%d ||R0||=%.2e ||R||=%.2e iters=%d flag=%d\n', ...
        prob.k, seedResNorm, out.resNorm, out.iterations, flag);
end
