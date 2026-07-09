function J = ms_jacobian_cs(yJ, prob, M)
% MS_JACOBIAN_CS  Sparse MS Jacobian by per-arc complex-step STMs.
%
% For each arc k, perturbs each varying component of the arc's initial
% state by 1i*hCS*scale and re-integrates: column = imag(yEnd)/(hCS*scale).
% Arc 1 varies only its 7 costates (state fixed at rv0, m0); arcs 2..M vary
% all 14 components. Blocks: continuity row k has +Phi_k (w.r.t. arc-k
% unknowns) and -I (w.r.t. y_{k+1}); terminal rows take Phi_M rows [1:6,14].
%
% INPUTS:
%   yJ   - initial augmented state of each arc [14xM] (from MS_UNPACK)
%   prob - problem struct with tJ [1x(M+1)]
%   M    - number of arcs [scalar]
%
% OUTPUTS:
%   J - sparse Jacobian [(14M-7)x(14M-7)]
%
% REFERENCES:
%   [1] Martins, Sturdza, Alonso, ACM TOMS 29(3), 2003 (complex step).

hCS = 1e-20;
Phi = cell(1, M);
for k = 1:M            % NOTE: safe to change to parfor if a pool is open
    if k == 1, dirs = 8:14; else, dirs = 1:14; end
    Pk = zeros(14, numel(dirs));
    for d = 1:numel(dirs)
        idx      = dirs(d);
        scale    = max(1, abs(yJ(idx, k)));
        yp       = complex(yJ(:, k));
        yp(idx)  = yp(idx) + 1i*hCS*scale;
        [~, Yc]  = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob.Tmax, ...
                   prob.c, prob.muStar, prob.epsSmooth), ...
                   [prob.tJ(k) prob.tJ(k+1)], yp, prob.odeOpts);
        Pk(:, d) = imag(Yc(end, :).')./(hCS*scale);
    end
    Phi{k} = Pk;
end

colOf = @(k) 7 + 14*(k - 2);            % column offset of unknown y_k, k >= 2
ii = []; jj = []; vv = [];
for k = 1:M-1
    r0 = 14*(k-1);
    if k == 1, c0 = 0; else, c0 = colOf(k); end
    [ii, jj, vv] = addBlock(ii, jj, vv, r0, c0, Phi{k});
    [ii, jj, vv] = addBlock(ii, jj, vv, r0, colOf(k+1), -eye(14));
end
[ii, jj, vv] = addBlock(ii, jj, vv, 14*(M-1), colOf(M), Phi{M}([1:6 14], :));
J = sparse(ii, jj, vv, 14*M - 7, 14*M - 7);
end

% -------------------------------------------------------------------------
function [ii, jj, vv] = addBlock(ii, jj, vv, r0, c0, B)
[nR, nC] = size(B);
[br, bc] = ndgrid(1:nR, 1:nC);
ii = [ii; r0 + br(:)];
jj = [jj; c0 + bc(:)];
vv = [vv; B(:)];
end
