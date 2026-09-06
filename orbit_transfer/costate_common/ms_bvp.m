function [p, info] = ms_bvp(prob, seed, opts)
%% Purpose:
%
%   GENERIC multiple-shooting two-point BVP engine -- the family- and
%   problem-agnostic core of ms_tfmin, moved to the shared library
%   (migration #3). The arc is split into K segments whose junction states
%   are extra unknowns; short segments kill the Lyapunov amplification that
%   makes single shooting from approximate seeds intractable (measured on
%   the CR3BP min-time problem: collocation seeds miss by 36,000-560,000 km
%   single-shot, converge in a few iterations here).
%
%   The PROBLEM is supplied as three closures (prob struct below), so the
%   engine owns only structure: the unknown vector p = [y1(freeIdx0);
%   y2..yK; tf], the block-bidiagonal-bordered Jacobian assembled from
%   segment STMs, the trust-region-dogleg solve, and the iterate sanity
%   guards. Breakpoints stay FIXED as fractions of tf.
%
%  ASSUMPTIONS / NOTES:
%
% • The rejection residual (opts.rejectR) must exceed ANY residual a real
%   iterate can reach, else a solver whose current ||R|| is larger reads
%   the guard as an improvement and ACCEPTS the garbage iterate (review
%   finding, Gemini 2026-08-04).
% • prob.prop must THROW on integrator collapse; the engine converts the
%   throw into a rejected iterate.
% • Final time tf is an unknown by default (free-tf problems). With
%   opts.fixedTf = true the trailing tf unknown is DROPPED, tf = seed.tf is
%   held constant, and prob.terminal must then supply exactly nf conditions
%   (min-energy / fixed-time catalogs). This is a real structural switch,
%   not a tight tf guard (which would leave a singular Jacobian column).
%
%% Inputs:
%
%  prob                     struct
%   .ny                     double                  State dimension (14 for
%                                                   CR3BP PMP: [rv; m; lam])
%   .freeIdx0               [1 x nf]                Indices of y1 that are
%                                                   UNKNOWNS (the fixed
%                                                   complement comes from
%                                                   seed.Y(:,1))
%   .prop                   fhandle                 [yh, PHI] = prop(dt, y0,
%                                                   needSTM): propagate y0
%                                                   for dt; PHI [ny x ny]
%                                                   when needSTM (else [])
%   .rhs                    fhandle                 F = rhs(y): dynamics at
%                                                   a point [ny x 1]
%   .terminal               fhandle                 [g, dgdy] = terminal(y,
%                                                   needJ): terminal
%                                                   residual [nT x 1] and
%                                                   its Jacobian [nT x ny]
%
%  seed                     struct
%   .tf                     double                  Time-of-flight guess
%   .tGrid                  [1 x K+1]               Segment boundary times,
%                                                   0..tf
%   .Y                      [ny x K+1]              States at the boundary
%                                                   times (col 1: fixed
%                                                   part used as-is)
%
%  opts                     struct (optional)       .maxIter [100], .tolR
%                                                   [1e-10], .wallSec [300],
%                                                   .verbose [false],
%                                                   .tfLo [0.3], .tfHi [3]
%                                                   (x seed.tf guards),
%                                                   .pMax [1e5], .rejectR
%                                                   [1e8], .keepSTMs [false],
%                                                   .fixedTf [false] (drop
%                                                   the tf unknown; tf =
%                                                   seed.tf), .polishMax
%                                                   [5] Newton polish steps
%                                                   after an early fsolve
%                                                   exit
%
%% Outputs:
%
%  p                        [n x 1]                 Best iterate:
%                                                   [y1(freeIdx0); y2..yK;
%                                                   tf], n = nf+ny(K-1)+1
%                                                   (fixedTf: no trailing
%                                                   tf, n = nf+ny(K-1))
%
%  info                     struct                  .converged, .normR,
%                                                   .iters, .wall (s),
%                                                   .tGrid (scaled by tf),
%                                                   .Y [ny x K junction
%                                                   starts], .PHI {1 x K}
%                                                   segment STMs at the
%                                                   final iterate (only if
%                                                   keepSTMs)
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  M. Casey  fixed-tf variant (opts.fixedTf) for min-energy    (c) 08/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: linear oscillator y'' = -y as a BVP: y(0)=0 fixed, y'(0) free,
   %hit y(tf)=1 with y'(tf)=0 at free tf. Exact answer: y=sin(t), tf=pi/2.
     A = [0 1; -1 0];
     pr = struct('ny',2, 'freeIdx0',2, ...
         'prop', @(dt,y0,nS) demo_prop(A,dt,y0,nS), ...
         'rhs',  @(y) A*y, ...
         'terminal', @(y,nJ) deal([y(1)-1; y(2)], eye(2)));
     K = 4;  tg = linspace(0, 1.2, K+1);
     sd = struct('tf', 1.2, 'tGrid', tg, 'Y', [sin(tg); cos(tg)]);
     [p_, inf_] = ms_bvp(pr, sd, struct('tolR',1e-8));
     fprintf('demo: y''(0) = %.9f (exact 1), tf = %.9f (exact pi/2 = %.9f)\n', ...
             p_(1), p_(end), pi/2);
     fprintf('      converged = %d, ||R|| = %.1e, iters = %d\n', ...
             inf_.converged, inf_.normR, inf_.iters);
   %Demo 2: FIXED tf = pi/2, y(0)=0 fixed, y'(0) free, hit y(tf)=1 only
   %(one terminal condition for one unknown). Exact answer: y'(0) = 1.
     prF = pr;  prF.terminal = @(y,nJ) deal(y(1)-1, [1 0]);
     tgF = linspace(0, pi/2, K+1);
     sdF = struct('tf', pi/2, 'tGrid', tgF, 'Y', 0.9*[sin(tgF); cos(tgF)]);
     [pF, infF] = ms_bvp(prF, sdF, struct('tolR',1e-8,'fixedTf',true));
     fprintf('demo (fixed tf): y''(0) = %.9f (exact 1), n = %d unknowns\n', ...
             pF(1), numel(pF));
     fprintf('      converged = %d, ||R|| = %.1e, iters = %d\n', ...
             infF.converged, infF.normR, infF.iters);
     return
end

if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
maxIter  = d('maxIter', 100);
tolR     = d('tolR', 1e-10);
wallSec  = d('wallSec', 300);
verbose  = d('verbose', false);
tfLo     = d('tfLo', 0.3);
tfHi     = d('tfHi', 3);
pMax     = d('pMax', 1e5);
rejectR  = d('rejectR', 1e8);
keepSTMs = d('keepSTMs', false);
fixedTf  = d('fixedTf', false);

ny  = prob.ny;
fi0 = prob.freeIdx0(:)';
nf  = numel(fi0);
K   = numel(seed.tGrid) - 1;
sig = seed.tGrid(:)'/seed.tGrid(end);          % fixed normalized breakpoints
dsg = diff(sig);                               % [1 x K]
y1fix = seed.Y(:,1);                           % fixed part of column 1

% unknown vector p = [y1(freeIdx0); Y_2..Y_K (ny each); tf]
% (fixedTf: the trailing tf is absent and tf = seed.tf throughout)
p = [seed.Y(fi0,1); reshape(seed.Y(:,2:K), ny*(K-1), 1)];
if ~fixedTf, p(end+1) = seed.tf; end
n = numel(p);                                  % = nf + ny(K-1) (+1 free tf)

% Trust-region-dogleg on the multiple-shooting system with the analytic
% block Jacobian. A plain Newton + backtracking loop stalls on rough seeds
% (measured: first step rejected on 2 of 3 pilot cells).
tStart = tic;
fopts = optimoptions('fsolve', 'Display', dispmode(verbose), ...
    'Algorithm','trust-region-dogleg', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', tolR^2, 'StepTolerance', 1e-13, ...
    'MaxIterations', maxIter, 'MaxFunctionEvaluations', 10*maxIter, ...
    'OutputFcn', @(x,ov,st) toc(tStart) > wallSec);   % wall budget
[p, ~, ~, fout] = fsolve(@residual, p, fopts);
[R, J] = residual(p);
normR = norm(R, inf);
iters = fout.iterations;

% NEWTON POLISH. fsolve's own exit test (first-order optimality ||J'R|| <
% 1e-6) fires early on a short, well-conditioned arc where ||J|| is small:
% measured on the min-energy binding, "Equation solved" at ||R|| = 3.9e-10
% with tolR = 1e-10 -- one Newton step short. Disabling that test instead
% makes fsolve grind at the residual floor (golden DRO cell: 1 -> 7 iters
% for the same 1.2e-13). So: leave fsolve alone, and when it returns above
% tolR take plain Newton steps with the analytic Jacobian while they help
% (at most polishMax, default 5). Already-converged runs are untouched.
polishMax = d('polishMax', 5);
kp = 0;
while normR > tolR && kp < polishMax && toc(tStart) < wallSec
    pTry = p - J\R;
    [Rt, Jt] = residual(pTry);
    nt = norm(Rt, inf);
    if ~(nt < normR), break, end             % no improvement (or NaN): stop
    p = pTry;  R = Rt;  J = Jt;  normR = nt;
    kp = kp + 1;
end
iters = iters + kp;

% Jacobian conditioning at the final iterate (converged or stalled): the
% fold-vs-grazing diagnostic for continuation walls (MfMax review
% 2026-09-06). A singular dS/dz at a turning point in the continuation
% parameter shows up as cond(J) blowing up as the wall is approached.
if all(isfinite(J(:))), condJ = cond(J); else, condJ = Inf; end
info = struct('converged', normR < tolR, 'normR', normR, ...
              'iters', iters, 'wall', toc(tStart), 'condJ', condJ, ...
              'tGrid', sig*tfOf(p), 'Y', junctions(p));
if keepSTMs
    [info.PHI, info.Yend] = collectSTMs(p);
end

% ------------------------------------------------------------------------
    function tf = tfOf(p)
    % TFOF  Final time of an iterate: the trailing unknown, or the fixed
    % seed value.  INPUTS: p [n x 1].  OUTPUTS: tf double.
    if fixedTf, tf = seed.tf; else, tf = p(end); end
    end

    function Yj = junctions(p)
    % JUNCTIONS  Junction START states from the unknown vector.
    % INPUTS: p [n x 1].  OUTPUTS: Yj [ny x K].
    y1 = y1fix;  y1(fi0) = p(1:nf);
    Yj = [y1, reshape(p(nf+(1:ny*(K-1))), ny, K-1)];
    end

    function [PHI, Yend] = collectSTMs(p)
    % COLLECTSTMS  Segment STMs at the final iterate (for the conjugate-
    % point test) and the terminal state y(tf) (so the test can sample the
    % final segment's flow).  INPUTS: p.  OUTPUTS: PHI {1 x K} of
    % [ny x ny]; Yend [ny x 1].
    Yj = junctions(p);  tf = tfOf(p);
    PHI = cell(1, K);
    Yend = [];
    for k = 1:K
        [Yend, PHI{k}] = prob.prop(dsg(k)*tf, Yj(:,k), true);
    end
    end

    function [R, J] = residual(p)
    % RESIDUAL  Multiple-shooting residual and (optionally) dense Jacobian.
    % INPUTS:  p [n x 1] unknown vector.  OUTPUTS: R [n x 1]; J [n x n].
    tf = tfOf(p);
    needJ = nargout > 1;
    % Iterate sanity guard: a wild trust-region iterate can drive a segment
    % propagation into integrator collapse (step size -> eps, unbounded
    % memory -- observed killing the process). Reject it with the large
    % residual instead of propagating it. (The tf window is vacuous when
    % tf is fixed.)
    if ~all(isfinite(p)) || tf < tfLo*seed.tf || tf > tfHi*seed.tf ...
            || max(abs(p)) > pMax
        R = rejectR*ones(n,1);
        if needJ, J = eye(n); end   % placeholder; step is rejected on ratio
        return
    end
    Yj = junctions(p);
    R = zeros(n,1);
    if needJ, J = zeros(n,n); end
    for k = 1:K
        try
            [yh, PHIk] = prob.prop(dsg(k)*tf, Yj(:,k), needJ);
        catch
            R = rejectR*ones(n,1);
            if needJ, J = eye(n); end
            return
        end
        Fh = prob.rhs(yh);
        colsK = colidx(k);
        if k < K
            rows = (k-1)*ny + (1:ny);
            R(rows) = yh - Yj(:,k+1);
            if needJ
                J(rows, colsK)       = segsel(PHIk, k);
                J(rows, colidx(k+1)) = -eye(ny);
                if ~fixedTf, J(rows, n) = J(rows, n) + Fh*dsg(k); end
            end
        else
            [g, dgdy] = prob.terminal(yh, needJ);
            rows = (K-1)*ny + (1:numel(g));
            R(rows) = g;
            if needJ
                J(rows, colsK) = segsel(dgdy*PHIk, k);
                if ~fixedTf, J(rows, n) = J(rows, n) + dgdy*Fh*dsg(k); end
            end
        end
    end
    end

    function c_ = colidx(k)
    % COLIDX  Columns of the unknowns feeding segment k's initial state.
    % INPUTS: k segment index. OUTPUTS: c_ column index vector.
    if k == 1, c_ = 1:nf;                      % free part of y1 only
    else,      c_ = nf + (k-2)*ny + (1:ny);
    end
    end

    function B = segsel(M, k)
    % SEGSEL  Restrict a d(...)/d(y_k) block to segment k's unknowns.
    % INPUTS: M [* x ny]; k index. OUTPUTS: B block.
    if k == 1, B = M(:, fi0);
    else,      B = M;
    end
    end
end

% ------------------------------------------------------------------------
function [yh, PHI] = demo_prop(A, dt, y0, needSTM)
% DEMO_PROP  Exact propagator for the linear demo.  INPUTS: A [2x2]; dt;
% y0 [2x1]; needSTM logical.  OUTPUTS: yh [2x1]; PHI [2x2] or [].
E = expm(A*dt);
yh = E*y0;
if needSTM, PHI = E; else, PHI = []; end
end

% ------------------------------------------------------------------------
function m = dispmode(verbose)
% DISPMODE  fsolve display mode from the verbose flag.
% INPUTS: verbose logical. OUTPUTS: m [char].
if verbose, m = 'iter'; else, m = 'off'; end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
