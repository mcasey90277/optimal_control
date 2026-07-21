function ic = psr_ipopt_certify(solFile, opts)
% PSR_IPOPT_CERTIFY  Local-minimality certificate from IPOPT's NATIVE inertia.
%
% The robust second-order certifier for this problem. psr_second_order.m
% reconstructs the KKT Hessian and factorizes it OURSELVES -- and over a 40-rev
% spiral that matrix is so ill-conditioned (cond ~1e9-1e16) that the LDL pivot
% signs near the noise floor are unreliable, producing hundreds of spurious
% tiny-magnitude "negative" eigenvalues (see the tolEig sweep in that file's
% history). IPOPT sidesteps this entirely: its inertia-controlled linear solver
% (MUMPS) checks the KKT inertia at EVERY iteration on the well-scaled system,
% and adds a Hessian regularization delta_w ONLY when the reduced Hessian is
% indefinite. If IPOPT reaches convergence with delta_w = 0 on the final
% iterations, the reduced Hessian is positive definite WITHOUT correction --
% i.e. the second-order sufficient conditions hold. That is a clean,
% conditioning-robust local-minimality certificate.
%
% Method: reconstruct the exact NLP, warm-start IPOPT AT the given solution (it
% converges in a handful of iterations), and read the per-iteration Hessian
% regularization (delta_w) from the solver stats. Certify if delta_w is
% negligible over the converged tail.
%
% SCOPE: delta_w=0 at convergence certifies the barrier problem's 2nd-order
% sufficient conditions at the final (small) mu. For eps>0 (strictly convex in
% the throttle) this is a STRICT local minimum of the NLP. For eps=0 (bang-bang,
% linear in the throttle) the barrier supplies the throttle-direction curvature,
% so it certifies a (non-strict / weak) local minimum -- strictness there lives
% in the switching times (Maurer-Osmolovskii), a separate question. Either way,
% delta_w=0 means "no descent direction of negative curvature" = a local min.
%
% INPUTS:
%   solFile - solution .mat in seed layout: out.X/out.U, sigma, tauf0, rv0,
%             rvf, factor
%   opts    - (optional) struct:
%             eps      - homotopy eps the solution was solved at [0]
%             maxIter  - IPOPT cap for the warm re-solve [100]
%             tailN    - # of final iterations that must have ~0 regularization
%                        [5]
%             regTol   - delta_w below this counts as ZERO (no regularization)
%                        [1e-8]
%             verbose  - [true]
%
% OUTPUTS:
%   ic - struct:
%     .certLocalMin  - logical: converged AND delta_w<=regTol over the tail
%     .converged     - the warm re-solve reached optimality
%     .nIter         - iterations the warm re-solve took
%     .regTail       - the last tailN regularization values (delta_w)
%     .regTailMax    - max delta_w over the tail
%     .defect        - max dynamics defect of the re-solved point (feasibility)
%     .verdict       - one-line human-readable verdict
%
% REFERENCES:
%   [1] Wachter & Biegler, "On the implementation of an interior-point ...
%       (IPOPT)," Math. Prog. 106 (2006) -- inertia correction (delta_w).
%   [2] Nocedal & Wright, Numerical Optimization 2e, Ch. 19 (interior point,
%       second-order conditions via inertia).
%   [3] PSR/psr_second_order.m (the ill-conditioned reconstruction this replaces
%       for certification).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'eps'),     opts.eps = 0;      end
if ~isfield(opts,'maxIter'), opts.maxIter = 100; end
if ~isfield(opts,'tailN'),   opts.tailN = 5;    end
if ~isfield(opts,'regTol'),  opts.regTol = 1e-8; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
vp = @(varargin) fprintf(varargin{:});
if ~opts.verbose, vp = @(varargin) []; end

S = load(solFile);
p = cr3bp_lt_params(0.025, 15, 2100);
tf = S.out.X(8, end);

% PREFERRED path: read IPOPT's regularization history from the solution's OWN
% solve (captured by casadi_minfuel_sundman). No re-solve, and it is the delta_w
% of the actual converged solve. Only re-solve if the file predates the capture.
if isfield(S.out,'regHistory') && ~isempty(S.out.regHistory)
    reg = S.out.regHistory(:).';
    ic.converged = true;      % a saved solution is a converged one
    if isfield(S.out,'maxDefect'), ic.defect = S.out.maxDefect; else, ic.defect = NaN; end
    ic.fromResolve = false;
    vp('psr_ipopt_certify: read delta_w history from the solution''s own solve (%d iters)\n', numel(reg));
else
    % Fallback: warm-start IPOPT AT the solution and re-converge (warmTight
    % honours the near-bang bounds, so no mint-cliff slide at eps=0).
    vp('psr_ipopt_certify: no stored regHistory -> warm re-solve at eps=%.3g...\n', opts.eps);
    o = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
            S.out.X, S.out.U, S.tauf0, 1.5, opts.maxIter, opts.eps, true);
    ic.converged = o.success && o.maxDefect < 1e-6;
    ic.defect = o.maxDefect;
    reg = o.regHistory(:).';
    ic.fromResolve = true;
end
ic.nIter = numel(reg);

if isempty(reg)
    ic.certLocalMin = false;
    ic.regTail = []; ic.regTailMax = NaN;
    ic.verdict = ['NO SIGNAL: IPOPT regularization history unavailable from the ' ...
        'solver stats (CasADi build without regularization_size) -- cannot read the ' ...
        'native inertia; fall back to psr_second_order or re-solve with logging.'];
    vp('  %s\n', ic.verdict);
    return
end

tailN = min(opts.tailN, numel(reg));
ic.regTail = reg(end-tailN+1:end);
ic.regTailMax = max(ic.regTail);

% Require at least tailN iterations: IPOPT initializes delta_w=0 BEFORE any
% inertia factorization, so a solve that terminates in ~0 iters (e.g. a warm
% re-solve that accepts the point immediately) returns regHistory=[0] -- the
% DEFAULT, not an inertia verdict. Certifying off that would be spurious.
ic.certLocalMin = ic.converged && (ic.regTailMax <= opts.regTol) && (ic.nIter >= opts.tailN);
if ~ic.converged
    ic.verdict = sprintf(['NOT CERTIFIED: the warm re-solve did not converge tight ' ...
        '(defect %.2e) -- the seed is not a clean KKT point.'], ic.defect);
elseif ic.nIter < opts.tailN
    ic.verdict = sprintf(['NO SIGNAL: only %d iteration(s) -- IPOPT accepted the point ' ...
        'before running >=%d inertia checks, so delta_w=%.1e is its pre-factorization ' ...
        'default, not a verdict. Re-solve from a slightly perturbed/looser warm start, ' ...
        'or read regHistory from the original solve.'], ic.nIter, opts.tailN, ic.regTailMax);
elseif ic.certLocalMin
    ic.verdict = sprintf(['LOCAL MIN (IPOPT native inertia): converged in %d iters with ' ...
        'delta_w = 0 over the final %d iterations (max %.1e <= %.1e). The reduced ' ...
        'Hessian is positive definite without regularization on the well-scaled system.'], ...
        ic.nIter, tailN, ic.regTailMax, opts.regTol);
else
    ic.verdict = sprintf(['NOT CERTIFIED: IPOPT added Hessian regularization delta_w up to ' ...
        '%.2e over the final %d iterations (> %.1e) -- the reduced Hessian needed ' ...
        'correction (indefinite/degenerate at convergence).'], ...
        ic.regTailMax, tailN, opts.regTol);
end
vp('  converged=%d  nIter=%d  regTailMax=%.2e  ->  %s\n', ...
   ic.converged, ic.nIter, ic.regTailMax, ic.verdict);
end
