function so = psr_switch_hessian(solFile, opts)
% PSR_SWITCH_HESSIAN  Switching-time (Maurer-Osmolovskii) 2nd-order test.
%
% The correct local-minimality certificate for a BANG-BANG solution -- the one
% psr_second_order.m cannot give because the raw NLP Hessian is degenerate for
% a linear-in-throttle objective. A bang-bang extremal's second-order behaviour
% lives in its SWITCHING TIMES, so this reduces the (infinite-dim) problem to a
% finite k-dimensional one over the k switch times and tests the projected
% Hessian there.
%
% Reduced problem. Freeze the arc structure (which arcs burn/coast) and the
% thrust-DIRECTION profile in arc-normalized coordinates (the direction already
% satisfies its own pointwise optimality -- the primer / strengthened
% Legendre-Clebsch condition, checked separately by the verifier). The only
% free parameters are the k switch locations sigma = (sigma_1,...,sigma_k) in
% the Sundman mesh variable. Integrate the state from the fixed start X0 =
% [rv0; m=1; t=0] arc by arc (throttle u_a on arc a, frozen direction alpha_a),
% giving the terminal state as a function of sigma:
%     J(sigma)     = propellant fraction = 1 - m(sigma_f)     (minimize)
%     c(sigma) = 0 : [r,v](sigma_f) = [rf,vf] and t(sigma_f) = t_f   (7 eqns)
% (sigma_f and sigma_0 are FIXED, so the arc spans always sum to the full
% Sundman length -- no separate total-time constraint.) At the optimum the
% reduced gradient dJ/dsigma + nu' dc/dsigma vanishes (this is exactly the
% switching condition S=0 at each switch, used here as a SELF-CHECK). Local
% minimality  <=>  the reduced Hessian of the Lagrangian L = J + nu' c,
% projected onto the feasible switch-time variations Z = null(dc/dsigma), is
% positive definite:  Z' (d^2 L / dsigma^2) Z  > 0.
%
% Derivatives use the complex-step method through an arc-by-arc integrator
% (each arc reparameterized to xi in [0,1] with a complex span factor, so the
% switch times can be perturbed complexly): dc/dsigma and dJ/dsigma are exact
% (imag/h); the Hessian is a finite-difference of the complex-step gradient of
% L. The 40-rev flow sensitivity is genuinely ill-conditioned (the same wall
% the indirect shooting hits) -- cond(dc/dsigma) is reported, and a large value
% weakens (does not invalidate) the projected-Hessian eigenvalues.
%
% SCOPE / HONESTY: this certifies local minimality of the reduced SWITCHING-TIME
% problem with the direction profile frozen at its (optimal) nominal. Together
% with (i) strict bang-bang Sdot != 0 at switches and (ii) primer/Legendre
% optimality of the direction -- both reported by the first-order verifier --
% a positive-definite projected switching-time Hessian is the Maurer-Osmolovskii
% sufficient condition for a bang-bang local minimum. It is NOT a global claim,
% and it inherits the base-point feasibility of the direct solution (reported).
%
% *** FINDING (2026-07-12, tested on the 1.15x refined solution) ***
% This FORWARD-FLOW formulation is BLOCKED by the 40-rev conditioning wall --
% the same instability that defeated the indirect finishing solve (IFS). The
% trapezoidal-collocation solution satisfies the DISCRETE defects but does NOT
% correspond to a nearby CONTINUOUS trajectory under forward integration: a
% high-accuracy re-integration of the exact direct control from X0 diverges by
% ||r||~3, ||v||~5, |t|~10 over the 40 revs. So the base point of this reduced
% problem is infeasible (||c|| ~ 20), the reduced gradient is not ~0, and any
% Hessian computed here is meaningless. The function DETECTS this (baseFeas
% self-check) and returns a BLOCKED verdict rather than a bogus certificate.
% The FIX (next build) is a MULTIPLE-SHOOTING / STM variational formulation:
% keep the collocation trajectory as the base (feasible), compute the
% switch-time sensitivities dPsi/dsigma_i = Phi(sigma_f,sigma_i)*[f_burn-f_coast]
% via the state-transition matrix integrated along that base, and form the
% second variation with segment matching -- so no quantity is forward-propagated
% across all 40 revs. That defeats the conditioning wall the way multiple
% shooting always does.
%
% INPUTS:
%   solFile - solution .mat in seed layout: out.X [8xnN], out.U [4xnN],
%             sigma, tauf0, rv0, rvf, factor
%   opts    - (optional) struct:
%             odeRelTol / odeAbsTol - flow tolerances [1e-9 / 1e-11]
%             hCS   - complex step [1e-20]
%             hFD   - Hessian finite-difference step (on switch times) [1e-6]
%             tolEig- reduced-Hessian PD threshold [1e-10 * scale]
%             gTol  - max reduced-gradient norm accepted in the self-check [1e-3]
%             verbose - [true]
%
% OUTPUTS:
%   so - struct:
%     .certLocalMin  - logical: projected switching-time Hessian is PD
%     .k             - number of switches
%     .redHessEig    - [1x(k-7)] eigenvalues of the projected Hessian (sorted)
%     .redHessMinEig - smallest eigenvalue (> 0 required)
%     .redGradNorm   - ||dJ/dsigma + nu' dc/dsigma|| (self-check, ~0 expected)
%     .baseFeas      - ||c(sigma*)|| at the nominal (direct-solution feasibility
%                      under continuous integration; ~collocation error)
%     .condDc        - cond(dc/dsigma) (40-rev sensitivity conditioning)
%     .nu            - [7x1] terminal multipliers (Lagrange, = terminal costates)
%     .verdict       - one-line human-readable verdict
%
% REFERENCES:
%   [1] Osmolovskii & Maurer, "Applications to Regular and Bang-Bang Control:
%       SSC in Calculus of Variations and Optimal Control," SIAM, 2012.
%   [2] Maurer, Buskens, Kim, Kaya, "Optimization methods for the verification
%       of second-order sufficient conditions for bang-bang controls," 2005.
%   [3] PSR/psr_second_order.m (the NLP-level test this complements).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'odeRelTol'), opts.odeRelTol = 1e-9;  end
if ~isfield(opts,'odeAbsTol'), opts.odeAbsTol = 1e-11; end
if ~isfield(opts,'hCS'),   opts.hCS   = 1e-20; end
if ~isfield(opts,'hFD'),   opts.hFD   = 1e-6;  end
if ~isfield(opts,'gTol'),  opts.gTol  = 1e-3;  end
if ~isfield(opts,'baseFeasTol'), opts.baseFeasTol = 1e-2; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
vp = @(varargin) fprintf(varargin{:});
if ~opts.verbose, vp = @(varargin) []; end

% ---- load + structure -------------------------------------------------------
S = load(solFile);
p = cr3bp_lt_params(0.025, 15, 2100);
Tmax = p.Tmax; c = p.c; muStar = p.muStar; pSund = 1.5;
rv0 = S.rv0(:); rvf = S.rvf(:); tauf0 = S.tauf0;
tf  = S.out.X(8,end);
[sigSw, arcU, sigGrid] = psr_switch_times(S);
k = numel(sigSw);
alphaNom = S.out.U(1:3, :);                     % frozen direction field (in sigma)
nomEdges = [sigGrid(1), sigSw, sigGrid(end)];   % NOMINAL arc edges (real, fixed)
X0 = [rv0; 1; 0];
odeO = odeset('RelTol', opts.odeRelTol, 'AbsTol', opts.odeAbsTol);
ws = warning('off','MATLAB:ode113:IntegrationTolNotMet');
cleanup = onCleanup(@() warning(ws));
vp('psr_switch_hessian: k=%d switches\n', k);

% ---- flow: [c(7); J] as a function of the (complex) switch times ------------
    function [cc, J] = flow(sigVec)
        edges = [sigGrid(1), sigVec(:).', sigGrid(end)];   % perturbed spans
        X = cast(X0, 'like', sigVec(1)+0);                 % promote if complex
        for a = 1:k+1
            dSig = edges(a+1) - edges(a);
            sigL = nomEdges(a);  sigR = nomEdges(a+1);      % nominal geom (real)
            u    = arcU(a);
            rhs  = @(xi, XX) arc_eom(xi, XX, dSig, u, sigL, sigR, ...
                       sigGrid, alphaNom, tauf0, Tmax, c, muStar, pSund);
            [~, Xo] = ode113(rhs, [0 1], X, odeO);
            X = Xo(end, :).';
        end
        cc = [X(1:6) - rvf; X(8) - tf];
        J  = 1 - X(7);
    end

% ---- nominal base feasibility (GUARD: 40-rev forward-shooting divergence) ---
[c0, ~] = flow(sigSw);
so.baseFeas = norm(c0);
so.k = k;
if so.baseFeas > opts.baseFeasTol
    so.certLocalMin = false;
    so.condDc = NaN; so.redGradNorm = NaN; so.nu = [];
    so.redHessEig = []; so.redHessMinEig = NaN;
    so.verdict = sprintf(['BLOCKED (40-rev forward-shooting divergence): re-integrating ' ...
        'the direct control from X0 misses the terminal state by ||c||=%.2e (>> %.1e). ' ...
        'The collocation solution does not correspond to a nearby continuous trajectory ' ...
        'over 40 revs, so this forward-flow switching-time Hessian is ill-posed. ' ...
        'Fix: multiple-shooting / STM variational formulation (see header FINDING).'], ...
        so.baseFeas, opts.baseFeasTol);
    vp('  %s\n', so.verdict);
    return
end

% ---- complex-step Jacobian --------------------------------------------------
dc = zeros(7, k);  dJ = zeros(1, k);
for j = 1:k
    sp = complex(sigSw);  sp(j) = sp(j) + 1i*opts.hCS;
    [cj, Jj] = flow(sp);
    dc(:, j) = imag(cj)/opts.hCS;
    dJ(j)    = imag(Jj)/opts.hCS;
end
so.condDc = cond(dc);

% terminal multipliers nu: minimize ||dJ + nu' dc|| -> nu = -(dc') \ dJ'
nu = -(dc.') \ (dJ.');
gRed = dJ + (nu.'*dc);
so.redGradNorm = norm(gRed);  so.nu = nu;
vp('  base feasibility ||c||=%.2e   reduced-gradient ||dJ+nu''dc||=%.2e   cond(dc)=%.2e\n', ...
   so.baseFeas, so.redGradNorm, so.condDc);
if so.redGradNorm > opts.gTol
    warning('psr_switch_hessian:notStationary', ...
        ['reduced gradient %.2e > %.1e: the switch times are not first-order ' ...
         'optimal for this (frozen-direction) reduced problem -- Hessian verdict ' ...
         'is unreliable (base solution / direction-freezing mismatch)'], ...
        so.redGradNorm, opts.gTol);
end

% ---- reduced Hessian of L = J + nu'c via FD of the complex-step gradient -----
    function g = gradL(sigVec)
        g = zeros(k, 1);
        for jj = 1:k
            sp2 = complex(sigVec);  sp2(jj) = sp2(jj) + 1i*opts.hCS;
            [cc, JJ] = flow(sp2);
            g(jj) = imag(JJ + nu.'*cc)/opts.hCS;
        end
    end

vp('  building %dx%d switching-time Hessian (FD of complex-step gradient)...\n', k, k);
g0 = gradL(sigSw);
H = zeros(k, k);
for j = 1:k
    sp = sigSw;  sp(j) = sp(j) + opts.hFD;
    H(:, j) = (gradL(sp) - g0)/opts.hFD;
end
H = (H + H.')/2;

% ---- project onto feasible switch-time variations + eigen-test --------------
Z = null(dc);                                   % k x (k - rank(dc))
so.k = k;
if isempty(Z)
    % k <= #terminal-constraints (7): the terminal constraints pin every
    % switch-time direction, so there are NO free directions to test. The
    % switching-time problem is (vacuously) locally minimal iff first-order
    % optimality holds. No Hessian eig-test is defined; guard the scalar-&&.
    so.certLocalMin  = (so.redGradNorm <= opts.gTol);
    so.redHessEig    = [];
    so.redHessMinEig = Inf;
    so.verdict = sprintf(['VACUOUS (k=%d <= 7 terminal constraints): no free ' ...
        'switch-time directions after the terminal rendezvous; first-order ' ...
        'optimality (||g||=%.2e) is the whole certificate.'], k, so.redGradNorm);
    vp('  %s\n', so.verdict);
    return
end
Hr = Z.'*H*Z;
Hr = (Hr + Hr.')/2;
ev = sort(real(eig(Hr)));
so.redHessEig = ev(:).';
so.redHessMinEig = min(ev);
scaleH = max(median(abs(diag(Hr))), 1);
tolEig = 1e-10 * scaleH;  if isfield(opts,'tolEig'), tolEig = opts.tolEig; end

so.certLocalMin = (so.redHessMinEig > tolEig) && (so.redGradNorm <= opts.gTol);
if so.certLocalMin
    so.verdict = sprintf(['LOCAL MIN (switching-time SSC): projected Hessian PD, ' ...
        'min eig %.3e over %d free switch-time directions (cond(dc)=%.1e)'], ...
        so.redHessMinEig, size(Z,2), so.condDc);
elseif so.redGradNorm > opts.gTol
    so.verdict = sprintf(['UNRELIABLE: switch times not first-order optimal for the ' ...
        'frozen-direction reduced problem (||g||=%.2e)'], so.redGradNorm);
elseif so.redHessMinEig <= tolEig && so.redHessMinEig > -tolEig
    so.verdict = sprintf(['INCONCLUSIVE: projected Hessian is degenerate (min eig ' ...
        '%.3e ~ 0); cond(dc)=%.1e may be masking curvature'], so.redHessMinEig, so.condDc);
else
    so.verdict = sprintf(['NOT a local min (switching-time): projected Hessian has a ' ...
        'negative eigenvalue %.3e over the feasible switch-time directions'], so.redHessMinEig);
end
vp('  %s\n', so.verdict);
end

% =============================================================================
function dX = arc_eom(xi, X, dSig, u, sigL, sigR, sigGrid, alphaNom, ...
                      tauf0, Tmax, c, muStar, pSund)
% Arc dynamics in the normalized coordinate xi in [0,1]: dX/dxi = dSig * dX/dsigma.
% The thrust direction is FROZEN (nominal), sampled at the nominal sigma mapped
% from xi (real), so complex perturbations of the switch times flow only through
% the complex span factor dSig -- keeping the complex-step derivative exact.
r = X(1:3);  v = X(4:6);  m = X(7);
dd = [r(1)+muStar; r(2); r(3)];
rr = [r(1)-1+muStar; r(2); r(3)];
d1 = sqrt(sum(dd.^2));
d3 = d1^3;  r3 = sqrt(sum(rr.^2))^3;
gr = [r(1); r(2); 0] - (1-muStar)*dd./d3 - muStar*rr./r3;
hv = [2*v(2); -2*v(1); 0];
kap = d1^pSund;
if u ~= 0
    nomSig = sigL + xi*(sigR - sigL);           % nominal (real) location
    al = interp1(sigGrid.', alphaNom.', nomSig, 'linear', 'extrap').';
    accel = gr + hv + u*(Tmax/m).*al;
else
    accel = gr + hv;
end
dX = dSig * tauf0 * kap * [v; accel; -u*Tmax/c; 1];
end
