function R = pmp_residual(sg, X, U, lam, fFun, LFun, opts)
% PMP_RESIDUAL  Continuous two-point-BVP residuals of a direct-collocation solution.
%
% Answers a question none of the existing verification layers does: the solver
% certifies the NLP, the defect check certifies the DISCRETE dynamics, and
% foc_check certifies first-order conditions OF THE NLP -- but nothing asks
% whether the trajectory, read as a continuous object, satisfies the continuous
% necessary conditions of the indirect (Pontryagin) formulation.
%
% WHY THIS IS EASIER THAN A MESH LADDER. A mesh-convergence study is a
% COMPARATIVE measurement: it needs two or more solutions and therefore a proof
% that they are the same solution branch. That proof is what defeated the
% ladder (Phase 2, shelved). This is an ABSOLUTE measurement of ONE solution,
% so no branch question arises.
%
% PROPAGATION IS PER-INTERVAL, NEVER GLOBAL. The costate dynamics are
% exponentially unstable: a single forward shot from lam(0) across a
% many-revolution transfer diverges and measures nothing. This repo has hit
% that wall twice (the 40-rev forward-shooting wall; the IFS wall at
% ||R|| ~ 0.023). Here each interval is integrated from the SOLUTION'S OWN
% values at its left node and compared at its right node, so nothing
% accumulates and the output is a residual MAP over the trajectory -- which
% localizes where the discretization struggles, more useful than a scalar.
%
% TWO CONTROL VARIANTS, REPORTED SEPARATELY, because they answer different
% questions and conflating them hides the answer to both:
%   'interp' - the control is the transcription's own, linearly interpolated
%              across the interval. Isolates the accuracy of the PRIMAL
%              propagation: how much does the trapezoid miss between nodes?
%   'pmp'    - the control is recomputed from the costate by the minimum
%              condition u = argmin_u H. Additionally tests whether the stored
%              control IS the Pontryagin control. The difference between the
%              two variants is therefore attributable to the minimum condition,
%              not to the dynamics.
%
% THE COSTATE DYNAMICS ARE AD-DERIVED, NOT HAND-WRITTEN. With
% H = L(x,u) + lam'*f(x,u), the adjoint equation is
%
%     dlam/dsigma = -(dH/dx)' = -( dL/dx' + (df/dx)' * lam )
%
% and both Jacobians are taken by CasADi from the SAME Function handles the
% caller passes -- ideally the very ones the transcription used. The check must
% be independent in INTEGRATION, not in the definition of the physics: a
% re-derived right-hand side would test the author's algebra rather than the
% mesh.
%
% INPUTS:
%   sg    - independent-variable grid, monotonic [(N+1)x1]
%   X     - states at nodes [nx x (N+1)]
%   U     - controls at nodes [nu x (N+1)]
%   lam   - costates at nodes [nx x (N+1)], already mapped from the NLP duals
%   fFun  - casadi.Function (x,u) -> dX/dsigma, or (x,u,s) when the dynamics
%           depend EXPLICITLY on the independent variable [nx x 1]. The arity
%           is detected, not declared. The MEE transcription is the
%           non-autonomous case: its node longitude is L = pi + sigma*DeltaL,
%           so the Gauss equations carry sigma explicitly.
%   LFun  - casadi.Function (x,u) or (x,u,s) -> running-cost integrand
%           dJ/dsigma [1x1]; pass [] for a pure Mayer problem (dL/dx = 0)
%   opts  - struct (optional):
%           .relTol .absTol   integrator tolerances [default 1e-12, 1e-14]
%           .variant  'interp' | 'pmp' | 'both'    [default 'both']
%           .uPMP     function handle (x,lam) -> u implementing the minimum
%                     condition; REQUIRED for the 'pmp' variant. Problem
%                     specific, so it is supplied rather than inferred.
%           .defectTol  self-gate threshold on reproducing the transcription's
%                     own trapezoidal defect [default 1e-8]
%           .verbose  [default true]
%
% OUTPUTS:
%   R - struct:
%       .Rx_interp   [1xN] ||x_prop - X(:,k+1)|| per interval, interpolated u
%       .Rlam_interp [1xN] ||lam_prop - lam(:,k+1)|| per interval
%       .Rx_pmp .Rlam_pmp  same under the PMP control law (NaN if not run)
%       .H           [1x(N+1)] Hamiltonian at each node
%       .Hvar        max|H - mean(H)| / max(|H|)  -- constancy diagnostic
%       .dx          [nx x N] SIGNED primal residual per interval, oriented
%                    (discrete node) - (exact flow from the previous node).
%                    .Rx_interp is its column norm. The sign is required for
%                    the objective-sensitivity estimate, which dots it with the
%                    defect multiplier.
%       .hLocal      [1xN] independent-variable step per interval
%       .defectCell  [1xN] per-interval norm of the recomputed DISCRETE defect.
%                    Plotted against .Rx it is the study's most direct picture:
%                    how much better the solution satisfies the discrete
%                    equations than the continuous ones.
%       .selfDefect  max trapezoidal defect recomputed from the caller's fFun
%       .gatePass    logical: selfDefect <= opts.defectTol
%       .nx .nu .N
%
% THE SELF-VALIDATING GATE. Before any residual is reported, the trapezoidal
% defects are recomputed from the supplied fFun and the stored X, U. If they
% are not machine-tight, the caller's fFun is NOT the transcription's dynamics
% and every residual below is meaningless. R.gatePass records this and the
% function warns loudly. A residual computed with the wrong right-hand side is
% worse than no residual, because it looks like evidence.
%
% REFERENCES:
%   [1] Pontryagin et al., "The Mathematical Theory of Optimal Processes," 1962.
%   [2] Betts, "Practical Methods for Optimal Control," 2nd ed., SIAM 2010,
%       ch. 4 (discretize-then-optimize versus the adjoint system).
%   [3] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, PHASE 3.

if nargin < 7, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
relTol    = d('relTol', 1e-12);
absTol    = d('absTol', 1e-14);
variant   = d('variant', 'both');
uPMP      = d('uPMP', []);
defectTol = d('defectTol', 1e-8);
verbose   = d('verbose', true);

sg = sg(:).';
N1 = numel(sg);  N = N1 - 1;
nx = size(X,1);  nu = size(U,1);
assert(size(X,2) == N1 && size(U,2) == N1 && size(lam,2) == N1, ...
    'pmp_residual:size', 'X, U, lam must all have %d columns', N1);
assert(size(lam,1) == nx, 'pmp_residual:size', 'lam must have %d rows', nx);

import casadi.*

% Arity detection: a non-autonomous transcription passes (x,u,s). Wrapping the
% autonomous case rather than branching everywhere keeps one code path below.
fNonAut = (fFun.n_in() >= 3);
LNonAut = ~isempty(LFun) && (LFun.n_in() >= 3);

% --- AD-derived adjoint right-hand side --------------------------------------
xs = MX.sym('x', nx);  us = MX.sym('u', nu);  ls = MX.sym('l', nx);
ss = MX.sym('s', 1);
if fNonAut, fx = fFun(xs, us, ss); else, fx = fFun(xs, us); end
Hs = ls.' * fx;
if ~isempty(LFun)
    if LNonAut, Hs = Hs + LFun(xs, us, ss); else, Hs = Hs + LFun(xs, us); end
end
dHdx  = jacobian(Hs, xs).';                 % [nx x 1]
adjFn = Function('adj', {xs, us, ls, ss}, {-dHdx});
HFn   = Function('H',   {xs, us, ls, ss}, {Hs});
fEval = Function('fe',  {xs, us, ss}, {fx});

% --- SELF-VALIDATING GATE ----------------------------------------------------
% Recompute the transcription's own trapezoidal defects with the supplied fFun.
fN = zeros(nx, N1);
for k = 1:N1, fN(:,k) = full(fEval(X(:,k), U(:,k), sg(k))); end
defects = zeros(nx, N);
for k = 1:N
    h = sg(k+1) - sg(k);
    defects(:,k) = X(:,k+1) - X(:,k) - (h/2)*(fN(:,k) + fN(:,k+1));
end
R.defectCell = vecnorm(defects, 2, 1);      % per-interval discrete defect norm
R.selfDefect = max(abs(defects(:)));
R.gatePass   = R.selfDefect <= defectTol;
if ~R.gatePass
    warning('pmp_residual:rhsMismatch', ...
        ['recomputed trapezoidal defect is %.3e (> %.1e): the supplied fFun is ' ...
         'NOT the dynamics this solution was produced with. Every residual ' ...
         'below is meaningless -- fix the wiring before reading them.'], ...
        R.selfDefect, defectTol);
end

% --- Hamiltonian at the nodes ------------------------------------------------
H = zeros(1, N1);
for k = 1:N1, H(k) = full(HFn(X(:,k), U(:,k), lam(:,k), sg(k))); end
R.H    = H;
R.Hvar = max(abs(H - mean(H))) / max(max(abs(H)), realmin);

% --- per-interval propagation ------------------------------------------------
doInterp = any(strcmp(variant, {'interp','both'}));
doPMP    = any(strcmp(variant, {'pmp','both'}));
if doPMP && isempty(uPMP)
    warning('pmp_residual:noUPMP', ...
        'variant requests the PMP control but opts.uPMP was not supplied; skipping it');
    doPMP = false;
end

R.Rx_interp   = nan(1,N);  R.Rlam_interp = nan(1,N);
R.Rx_pmp      = nan(1,N);  R.Rlam_pmp    = nan(1,N);
R.dx          = nan(nx,N);          % SIGNED primal residual, per interval
R.hLocal      = diff(sg);

odeOpt = odeset('RelTol', relTol, 'AbsTol', absTol);
for k = 1:N
    a = sg(k);  b = sg(k+1);
    z0 = [X(:,k); lam(:,k)];
    zEnd = [X(:,k+1); lam(:,k+1)];

    if doInterp
        uOf = @(s) U(:,k) + (s - a)/(b - a) * (U(:,k+1) - U(:,k));
        [~, Z] = ode113(@(s,z) local_rhs(s, z, nx, uOf, fEval, adjFn), [a b], z0, odeOpt);
        zp = Z(end,:).';
        % SIGNED, oriented as (discrete) - (true flow): the amount by which the
        % discrete solution's next node departs from where the exact dynamics
        % would have carried it. This is the quantity the defect multiplier
        % prices, so the sign must be kept -- a norm cannot be dotted with a
        % costate, and cancellation between intervals is itself a result.
        R.dx(:,k)        = zEnd(1:nx) - zp(1:nx);
        R.Rx_interp(k)   = norm(R.dx(:,k));
        R.Rlam_interp(k) = norm(zp(nx+1:end)  - zEnd(nx+1:end));
    end
    if doPMP
        uOf = @(s) [];   % placeholder; the PMP law needs z, handled in local_rhs
        [~, Z] = ode113(@(s,z) local_rhs_pmp(s, z, nx, uPMP, fEval, adjFn), [a b], z0, odeOpt);
        zp = Z(end,:).';
        R.Rx_pmp(k)   = norm(zp(1:nx)     - zEnd(1:nx));
        R.Rlam_pmp(k) = norm(zp(nx+1:end) - zEnd(nx+1:end));
    end
end

R.nx = nx;  R.nu = nu;  R.N = N;

if verbose
    fprintf('\n=== pmp_residual: N = %d intervals, nx = %d ===\n', N, nx);
    fprintf('  self-gate defect      : %.3e   %s\n', R.selfDefect, ...
        local_tf(R.gatePass, 'PASS', 'FAIL -- residuals below are meaningless'));
    fprintf('  Hamiltonian mean/var  : %.6e  /  %.3e (relative)\n', mean(H), R.Hvar);
    if doInterp
        fprintf('  R_x   (interp u)      : max %.3e   median %.3e\n', ...
            max(R.Rx_interp), median(R.Rx_interp));
        fprintf('  R_lam (interp u)      : max %.3e   median %.3e\n', ...
            max(R.Rlam_interp), median(R.Rlam_interp));
    end
    if doPMP
        fprintf('  R_x   (PMP u)         : max %.3e   median %.3e\n', ...
            max(R.Rx_pmp), median(R.Rx_pmp));
        fprintf('  R_lam (PMP u)         : max %.3e   median %.3e\n', ...
            max(R.Rlam_pmp), median(R.Rlam_pmp));
    end
end
end

% ---------------------------------------------------------------------------
function dz = local_rhs(s, z, nx, uOf, fEval, adjFn)
% LOCAL_RHS  Augmented [state; costate] flow with a PRESCRIBED control.
% INPUTS: s - independent var; z [2nx x 1]; nx; uOf(s) control; fFun; adjFn
% OUTPUTS: dz [2nx x 1]
x = z(1:nx);  l = z(nx+1:end);  u = uOf(s);
dz = [full(fEval(x,u,s)); full(adjFn(x,u,l,s))];
end

% ---------------------------------------------------------------------------
function dz = local_rhs_pmp(s, z, nx, uPMP, fEval, adjFn)
% LOCAL_RHS_PMP  Augmented flow with the control recomputed from the costate by
% the minimum condition, rather than taken from the stored solution.
% INPUTS: z [2nx x 1]; nx; uPMP(x,lam); fFun; adjFn   OUTPUTS: dz [2nx x 1]
x = z(1:nx);  l = z(nx+1:end);  u = uPMP(x, l);
dz = [full(fEval(x,u,s)); full(adjFn(x,u,l,s))];
end

% ---------------------------------------------------------------------------
function s = local_tf(c, a, b)
% LOCAL_TF  Pick a label by condition. INPUTS: c,a,b  OUTPUTS: s [char]
if c, s = a; else, s = b; end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
