function [X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew)
% WARMSTART_ON_MESH  Build a no-resample warm start on a refined sigma mesh.
%
% Every original node's state/control is copied VERBATIM (the no-resample
% discipline); only inserted nodes are filled. States and thrust direction
% are pchip-interpolated (direction renormalized to unit norm); the throttle
% is STEP-held from the nearest original node to the LEFT (the pre-switch
% side), so an insert straddling a switch is never seeded with a smeared
% intermediate throttle -- the re-solve relocates the switch.
%
% INPUTS:
%   out      - struct with X [8xnN], U [4xnN] on the OLD mesh
%   sigma    - old normalized nodes [nNx1]
%   sigmaNew - refined normalized nodes [nN'x1], contains all of sigma
%   isNew    - logical mask of inserted nodes [nN'x1]
%
% OUTPUTS:
%   X0 - warm-start states on the refined mesh [8xnN']
%   U0 - warm-start controls on the refined mesh [4xnN']
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md
%   [2] OPTIMALITY_VERIFICATION_PLAN.md sec F.3 (no-resample requirement)

sigma = sigma(:);  sigmaNew = sigmaNew(:);  isNew = logical(isNew(:));
X = out.X;  U = out.U;
nNn = numel(sigmaNew);
X0 = zeros(8, nNn);  U0 = zeros(4, nNn);

% originals verbatim
X0(:, ~isNew) = X;
U0(:, ~isNew) = U;

ins = find(isNew).';
if ~isempty(ins)
    sv = sigmaNew(ins);
    % states + direction by pchip on the OLD (sigma, .) grid
    X0(:, ins)     = interp1(sigma, X.',        sv, 'pchip').';
    al             = interp1(sigma, U(1:3, :).', sv, 'pchip').';
    al             = al ./ sqrt(sum(al.^2, 1));
    U0(1:3, ins)   = al;
    % throttle: step-hold from the nearest original node to the left
    for q = ins
        kk = find(sigma <= sigmaNew(q), 1, 'last');
        kk = min(max(kk, 1), numel(sigma));
        U0(4, q) = U(4, kk);
    end
end
end
