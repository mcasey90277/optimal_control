function [z, info] = ms_tfmin(rv0, rvf, seed, Tmax, c, muStar, opts)
% MS_TFMIN  Multiple-shooting solve of the CR3BP minimum-time PMP problem.
%
% Solves the same problem as pumpkyn.cr3bp.tfMin (terminal conditions
% r(tf)=rf, v(tf)=vf, lambda_m(tf)=0, H(tf)=0 with H = 1 + lambda'f), but
% with the arc split into K short segments whose junction states are extra
% unknowns. Short segments kill the ~1e3x Lyapunov amplification that makes
% single shooting from collocation-derived costates intractable (measured:
% catalog seeds miss by 36,000-560,000 km when single-shot). The seed is a
% full state+costate TRAJECTORY, e.g. from a direct collocation solution.
%
% All dynamics, STMs, and the Hamiltonian come from pumpkyn calls
% (tfMinProp / tfMinEoM) so every convention matches tfMin exactly; a
% converged z here is interchangeable with a tfMin solution.
%
% INPUTS:
%   rv0    - initial rotating-frame pos/vel, ND [6x1 or 1x6]
%   rvf    - final   rotating-frame pos/vel, ND [6x1 or 1x6]
%   seed   - struct:
%              .tf     initial time-of-flight guess, ND [scalar]
%              .tGrid  segment-boundary times, 0..tf [1 x K+1]
%              .Y      augmented states [r;v;m;lam_r;lam_v;lam_m] at the
%                      boundary times [14 x K+1] (col 1 state part is
%                      overwritten with [rv0; 1])
%   Tmax   - max thrust accel, ND [scalar]
%   c      - exhaust velocity, ND [scalar]
%   muStar - CR3BP mass ratio [scalar]
%   opts   - (optional) struct: .maxIter [30], .tolR [1e-10], .wallSec
%            [300], .verbose [false]
%
% OUTPUTS:
%   z    - [lambda0(7); tf] in tfMin's convention [8x1] (best iterate)
%   info - struct: .converged, .normR (final inf-norm), .iters, .wall (s),
%          .Y (14 x K junction START states; the final propagated endpoint
%          is not included), .tGrid
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMin / tfMinProp / tfMinEoM (D. Koblick, Coorbital) -
%       single-shooting original; conventions and dynamics reused verbatim.
%   [2] Betts, "Practical Methods for Optimal Control", ch. 3 (multiple
%       shooting structure).

if nargin < 7, opts = struct(); end
d = @(f,v) getfieldwithdefault(opts, f, v);
maxIter = d('maxIter', 100);
tolR    = d('tolR', 1e-10);
wallSec = d('wallSec', 300);
verbose = d('verbose', false);

rv0 = rv0(:);  rvf = rvf(:);
K   = numel(seed.tGrid) - 1;
sig = seed.tGrid(:)'/seed.tGrid(end);          % fixed normalized breakpoints
dsg = diff(sig);                               % [1 x K]

% unknown vector p = [lambda0(7); Y_2..Y_K (14 each); tf]
Y0 = seed.Y;
Y0(1:7,1) = [rv0; 1];
p = [Y0(8:14,1); reshape(Y0(:,2:K), 14*(K-1), 1); seed.tf];
n = numel(p);                                  % = 14K - 6

% Trust-region-dogleg (same globalization tfMin uses) on the multiple-
% shooting system, with the analytic block Jacobian. A plain Newton +
% backtracking loop stalls on rough seeds (measured: first step rejected on
% 2 of 3 pilot cells); the trust region handles those steps properly.
tStart = tic;
fopts = optimoptions('fsolve', 'Display', dispmode(verbose), ...
    'Algorithm','trust-region-dogleg', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', tolR^2, 'StepTolerance', 1e-13, ...
    'MaxIterations', maxIter, 'MaxFunctionEvaluations', 10*maxIter, ...
    'OutputFcn', @(x,ov,st) toc(tStart) > wallSec);   % wall budget
[p, ~, ~, fout] = fsolve(@residual, p, fopts);
[R, ~] = residual(p);
normR = norm(R, inf);
it = fout.iterations;

z = [p(1:7); p(end)];
info = struct('converged', normR < tolR, 'normR', normR, 'iters', it, ...
              'wall', toc(tStart), 'tGrid', sig*p(end), ...
              'Y', [ [rv0; 1; p(1:7)], reshape(p(8:end-1), 14, K-1)]);

% ------------------------------------------------------------------------
    function [R, J] = residual(p)
    % RESIDUAL  Multiple-shooting residual and (optionally) dense Jacobian.
    % INPUTS:  p [n x 1] unknown vector.  OUTPUTS: R [n x 1]; J [n x n].
    lam0 = p(1:7);
    tf   = p(end);
    needJ = nargout > 1;
    % Iterate sanity guard: a wild trust-region iterate can drive a segment
    % propagation into integrator collapse (step size -> eps, unbounded
    % memory -- observed killing the process). Reject it with a large
    % residual instead of propagating it.
    % The rejection residual must exceed ANY real accepted residual, else a
    % solver whose current ||R|| is larger would read it as an improvement
    % and ACCEPT the garbage iterate (review finding, Gemini 2026-08-04).
    if ~all(isfinite(p)) || tf < 0.3*seed.tf || tf > 3*seed.tf ...
            || max(abs(p)) > 1e5
        R = 1e8*ones(n,1);
        if needJ, J = eye(n); end   % placeholder; step is rejected on ratio
        return
    end
    Yj   = [ [rv0; 1; lam0], reshape(p(8:end-1), 14, K-1) ];   % junctions
    R = zeros(n,1);
    if needJ, J = zeros(n,n); end
    PHI0 = eye(14);
    for k = 1:K
        y0k = Yj(:,k);
        if needJ, y0k = [y0k; PHI0(:)]; end
        try
            [~, Yout] = pumpkyn.cr3bp.tfMinProp(dsg(k)*tf, y0k, Tmax, c, muStar);
        catch
            R = 1e8*ones(n,1);
            if needJ, J = eye(n); end
            return
        end
        yh = Yout(end,1:14)';
        if needJ, PHIk = reshape(Yout(end,15:210), 14, 14); end
        [Fh, Hh, dHdy] = pumpkyn.cr3bp.tfMinEoM(0, [yh; PHI0(:)], Tmax, c, muStar);
        Fh = Fh(1:14);
        colsK = colidx(k);
        if k < K
            rows = (k-1)*14 + (1:14);
            R(rows) = yh - Yj(:,k+1);
            if needJ
                J(rows, colsK)      = segJ(PHIk, k);
                J(rows, colidx(k+1)) = -eye(14);
                J(rows, n)          = J(rows, n) + Fh*dsg(k);
            end
        else
            rows = (K-1)*14 + (1:8);
            R(rows) = [yh(1:6) - rvf; yh(14); Hh];
            if needJ
                T = [PHIk(1:6,:); PHIk(14,:); dHdy*PHIk];      % [8 x 14]
                G = [Fh(1:6); Fh(14); dHdy*Fh];                % [8 x 1]
                J(rows, colsK) = tsel(T, k);
                J(rows, n)     = J(rows, n) + G*dsg(k);
            end
        end
    end
    end

    function c_ = colidx(k)
    % COLIDX  Columns of the unknowns feeding segment k's initial state.
    % INPUTS: k segment index. OUTPUTS: c_ column index vector.
    if k == 1, c_ = 1:7;                       % lambda0 only (state fixed)
    else,      c_ = 7 + (k-2)*14 + (1:14);
    end
    end

    function B = segJ(PHIk, k)
    % SEGJ  d(segment-k endpoint)/d(segment-k unknowns).
    % INPUTS: PHIk [14x14] segment STM; k index. OUTPUTS: B block.
    if k == 1, B = PHIk(:, 8:14);              % w.r.t. lambda0
    else,      B = PHIk;
    end
    end

    function B = tsel(T, k)
    % TSEL  Terminal-row block restricted to segment-K unknowns.
    % INPUTS: T [8x14]; k index. OUTPUTS: B block.
    if k == 1, B = T(:, 8:14);
    else,      B = T;
    end
    end
end

% ------------------------------------------------------------------------
function m = dispmode(verbose)
% DISPMODE  fsolve display mode from the verbose flag.
% INPUTS: verbose logical. OUTPUTS: m [char].
if verbose, m = 'iter'; else, m = 'off'; end
end

% ------------------------------------------------------------------------
function v = getfieldwithdefault(s, f, v0)
% GETFIELDWITHDEFAULT  s.(f) if present else v0.
% INPUTS: s struct; f field name; v0 default. OUTPUTS: v value.
if isfield(s, f), v = s.(f); else, v = v0; end
end
