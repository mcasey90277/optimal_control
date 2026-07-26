function ic = foc_ipopt_inertia(regHistory, opts)
% FOC_IPOPT_INERTIA  Local-minimality verdict from IPOPT's NATIVE inertia (delta_w).
%
% Campaign-agnostic interpreter for the per-iteration Hessian-regularization
% history CasADi/IPOPT report via opti.stats().iterations.regularization_size.
% IPOPT's inertia-controlled linear solver (MUMPS) checks the KKT inertia at
% EVERY iteration on the (well-scaled) barrier subproblem and adds a Hessian
% regularization delta_w ONLY when the reduced Hessian is indefinite. If the
% solve converges with delta_w = 0 over its final iterations, the reduced
% Hessian is positive definite WITHOUT correction there -- i.e. the
% second-order sufficient conditions hold at that point. That is a clean,
% conditioning-robust local-minimality certificate: it needs no ill-conditioned
% KKT-Hessian reconstruction of our own (contrast psr_second_order.m, which
% factorizes an explicit KKT matrix and is unreliable over a many-rev spiral,
% cond ~1e9-1e16).
%
% Pure interpretation: this function does NOT solve or re-solve anything. It
% only reads a regHistory vector the caller already has (from a solver's own
% out.regHistory, captured at solve time). Interpretation logic (tail window,
% tolerance, verdict phrasing) is ported from
% GTO_tulip/direct/certify/psr_ipopt_certify.m -- see that file's REFERENCES for
% the IPOPT inertia-correction citation. That function additionally handles
% the no-regHistory-stored fallback (warm re-solve); this generic version
% assumes the caller already has (or does not have) the history in hand and
% reports NO-DATA rather than re-solving.
%
% INPUTS:
%   regHistory - per-iteration delta_w values from the solve, or [] if the
%                solver predates capture / the CasADi build lacks
%                regularization_size [1 x nIter or []]
%   opts       - (optional) struct:
%                tailN - # of final iterations that must have ~0
%                        regularization [5]
%                tol   - delta_w below this counts as ZERO (no
%                        regularization) [1e-8]
%
% OUTPUTS:
%   ic - struct:
%     .certLocalMin - logical: regHistory nonempty AND max(tail) <= tol
%     .maxTailDw    - max delta_w over the tail window (NaN if regHistory empty)
%     .tail         - the last min(tailN,nIter) delta_w values ([] if empty)
%     .nIter        - numel(regHistory)
%     .verdict      - one-line human-readable verdict [char]
%
% REFERENCES:
%   [1] Wachter & Biegler, "On the implementation of an interior-point ...
%       (IPOPT)," Math. Prog. 106 (2006) -- inertia correction (delta_w).
%   [2] Nocedal & Wright, Numerical Optimization 2e, Ch. 19 (interior point,
%       second-order conditions via inertia).
%   [3] GTO_tulip/direct/certify/psr_ipopt_certify.m (the campaign-specific
%       precedent this interpreter is ported from; two of its verdict
%       phrasings are mirrored verbatim below).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'tailN') || isempty(opts.tailN), opts.tailN = 5;    end
if ~isfield(opts,'tol')   || isempty(opts.tol),   opts.tol   = 1e-8; end

reg = regHistory(:).';
ic.nIter = numel(reg);

if isempty(reg)
    ic.certLocalMin = false;
    ic.tail = [];
    ic.maxTailDw = NaN;
    ic.verdict = 'NO-DATA: regHistory absent (solver predates capture)';
    return
end

tailN = min(opts.tailN, ic.nIter);
ic.tail = reg(end-tailN+1:end);
ic.maxTailDw = max(ic.tail);
ic.certLocalMin = ic.maxTailDw <= opts.tol;

if ic.certLocalMin
    ic.verdict = sprintf(['LOCAL MIN (IPOPT native inertia): delta_w = 0 over the final ' ...
        '%d iterations (max %.1e <= %.1e). The reduced Hessian is positive definite ' ...
        'without regularization on the well-scaled system.'], tailN, ic.maxTailDw, opts.tol);
else
    ic.verdict = sprintf(['NOT CERTIFIED: IPOPT added Hessian regularization delta_w up to ' ...
        '%.2e over the final %d iterations (> %.1e) -- the reduced Hessian needed ' ...
        'correction (indefinite/degenerate at convergence).'], ic.maxTailDw, tailN, opts.tol);
end
end
