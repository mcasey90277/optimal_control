function [Z0, sigma, info] = build_guess(mode, N, rv0, rvf, Tmax, c, muStar)
% BUILD_GUESS  Initial decision vector + density-matched mesh fractions.
%
% Returns both the warm-start decision vector and the NONUNIFORM normalized
% mesh sigma used by the transcription. The mesh is density-matched to the
% reference arc: node fractions are quantiles of the reference integrator's
% own adaptive time grid (equal number of adaptive steps per segment), so
% mesh resolution follows the dynamics -- dense through perigee passes,
% sparse on slow arcs. A uniform mesh at practical N under-resolves perigee
% (defects ~1) and the NLP wanders.
%
% Modes:
%   'indirect'   - solve the indirect (PMP shooting) problem with
%                  pumpkyn.cr3bp.tfMin, propagate the optimal arc, and
%                  interpolate states + primer-direction controls onto the
%                  mesh. High-quality warm start: the NLP then only has to
%                  close the transcription (defect) error.
%   'tangential' - propagate max-thrust velocity-aligned steering
%                  (w = v/||v||, the classic orbit-raising heuristic) and
%                  mesh it. Independent of the indirect machinery, but the
%                  terminal state does NOT match the tulip target -- the
%                  NLP must close a large rendezvous gap through a
%                  many-revolution spiral; may converge slowly, to a
%                  different local minimum, or not at all. For
%                  experimentation.
%
% INPUTS:
%   mode   - 'indirect' or 'tangential' [char]
%   N      - number of trapezoidal segments [scalar]
%   rv0    - initial position/velocity (ND) [1x6]
%   rvf    - target position/velocity (ND) [1x6]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   Z0     - initial decision vector [10*(N+1)+1 x 1]
%   sigma  - normalized node times [(N+1)x1], sigma(1)=0, sigma(end)=1
%   info   - struct with .tf (guess transfer time) and, for 'indirect',
%            .zIndirect (converged [lambda0(7); tf])
%
% REFERENCES:
%   [1] pumpkynPie Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m (source of
%       the converged indirect costate guess).
%   [2] Betts, SIAM 2010, Ch. 4 (mesh refinement rationale).

nNodes = N + 1;

switch lower(mode)
    case 'indirect'
        % Converged costate/tf guess published in the pumpkynPie demo:
        zSeed = [ 190.476497248065; -79.7064866984696; -0.430399154713168;
                    0.301159446575878; 0.586671892449694;
                   -0.00711582435720301; 4.32931089137559; 6.29081541876621];
        zInd = pumpkyn.cr3bp.tfMin(rv0, rvf, zSeed, Tmax, c, muStar);
        tf   = zInd(8);
        [tauArc, yArc] = pumpkyn.cr3bp.tfMinProp(tf, [rv0(:).', 1, zInd(1:7).'], ...
                                                 Tmax, c, muStar);
        [tauArc, keep] = unique(tauArc, 'stable');
        yArc  = yArc(keep, :);
        sigma = density_matched_mesh(tauArc, nNodes);

        tMesh = sigma.*tf;
        Xg    = interp1(tauArc, yArc(:,1:7), tMesh, 'pchip').';
        lamV  = interp1(tauArc, yArc(:,11:13), tMesh, 'pchip').';
        Wg    = unit_columns(-lamV);               % primer direction
        info  = struct('tf', tf, 'zIndirect', zInd);

    case 'tangential'
        tf   = 6.3;                                % ballpark transfer time
        opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
        steer = @(t, x) lt_dynamics(x, x(4:6)./sqrt(sum(x(4:6).^2)), ...
                                    Tmax, c, muStar);
        [tauArc, yArc] = ode113(steer, [0 tf], [rv0(:); 1], opts);
        [tauArc, keep] = unique(tauArc, 'stable');
        yArc  = yArc(keep, :);
        sigma = density_matched_mesh(tauArc, nNodes);

        tMesh = sigma.*tf;
        Xg    = interp1(tauArc, yArc, tMesh, 'pchip').';
        Vg    = Xg(4:6, :);
        Wg    = unit_columns(Vg);
        info  = struct('tf', tf);

    otherwise
        error('build_guess:mode', 'Unknown guess mode: %s', mode);
end

Z0 = [Xg(:); Wg(:); tf];
end

% -------------------------------------------------------------------------
function sigma = density_matched_mesh(tauArc, nNodes)
% Node fractions as quantiles of the adaptive integration grid: segment k
% spans the same number of integrator steps everywhere, so node density
% tracks where the integrator worked hardest.
tauArc = unique(tauArc(:), 'stable');    % defensive: drop duplicates
if ~all(isfinite(tauArc)) || any(diff(tauArc) <= 0)
    error('build_guess:mesh', 'integrator abscissae not finite/increasing');
end
sigFull = (tauArc - tauArc(1))./(tauArc(end) - tauArc(1));
sigma   = interp1(linspace(0, 1, numel(sigFull)).', sigFull, ...
                  linspace(0, 1, nNodes).');
sigma(1)   = 0;
sigma(end) = 1;
if any(diff(sigma) <= 0)
    error('build_guess:mesh', 'mesh fractions not strictly increasing');
end
end

% -------------------------------------------------------------------------
function U = unit_columns(V)
% Normalize columns to unit length, guarding near-zero columns by carrying
% the previous column's direction (arbitrary +x for a zero first column).
mags = sqrt(sum(V.^2, 1));
U    = V;
tol  = 1e-12;
for kCol = 1:size(V, 2)
    if mags(kCol) > tol
        U(:, kCol) = V(:, kCol)./mags(kCol);
    elseif kCol > 1
        U(:, kCol) = U(:, kCol-1);
    else
        U(:, kCol) = [1; 0; 0];
    end
end
end
