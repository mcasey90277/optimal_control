function [Z, prob, meta] = ifs_seed_adjoint(matFile, opts)
% IFS_SEED_ADJOINT  IFS seed via a backward adjoint sweep / smoother (Rung A).
%
% Two methods (opts.method):
%
%   'sweep'  - pure backward sweep: init the costates at tau_f from min-fuel
%              terminal transversality + the direct NLP's terminal-constraint
%              multipliers, integrate the 8 costate ODEs backward tau_f->tau_0
%              with the state frozen to the direct solution. RESULT (2026-07-12,
%              1.12x): FALSIFIED as a seed -- the adjoint flow amplifies
%              ~1.2e12 over the 40-rev backward span, so a 1-2% terminal-dual
%              error rides to ~1e10 seed residual even though the terminal init
%              itself is good (q = 0.98 at the switch nearest tau_f). Kept for
%              diagnostics.
%
%   'smooth' - adjoint SMOOTHER (default): fit the terminal costate to the
%              WHOLE dual-map costate history, not just tau_f. Unknowns
%              z = [lam_f(rv 6); lam_f(t); s] (lam_m(tau_f)=0 pinned exact;
%              s = scale of the dual data, since beta is itself uncertain).
%              Residual = [component-scaled misfit to the dual history at
%              subsampled nodes; weighted switch rows q_s - 1] where
%              q = ||lamV||c/m + lamM must equal 1 at every switch (S=0).
%              Solved by Gauss-Newton with complex-step sweep sensitivities
%              and a truncated-SVD step (the ~1e12-amplified directions are
%              determined by mid-trajectory data where observable, truncated
%              where not) + backtracking line search. Node costates are then
%              BLENDED: sweep value where the local fit misfit < blendTol,
%              (s-scaled) dual value elsewhere -- worst case equals the dual
%              map instead of inheriting sweep blowup.
%
% The sweep RHS is the costate half of IFS_EOM evaluated at [x_dir(tau); lam]
% (pchip interpolation of out.X), so it is exactly consistent with the IFS
% residual's EOM by construction. Switch structure comes from the DUAL-map
% switching function's zero crossings (as in IFS_SEED -- the certified count;
% raw throttle crossings overcount, e.g. 12 vs the certified 10 at 1.12x, and
% are reported in meta only as a cross-check).
%
% INPUTS:
%   matFile - direct solution .mat (out.X [8xnN], out.U [4xnN], out.lamDef
%             [8xN], optionally out.lamAll, factor, tauf0, sigma)
%   opts    - struct (all optional):
%             method    - 'smooth' (default) | 'sweep'
%             tauParam  - 'sigmoid' (default) | 'direct' switch-time param
%             odeRelTol - sweep integrator RelTol [1e-11]
%             odeAbsTol - sweep integrator AbsTol [1e-13]
%             nData     - subsampled data nodes for the fit [256]
%             fitScale  - fit the data scale s (else pinned at beta0) [false]
%                         (free s drifted 5.6x on 1.12x -- switch rows too weak
%                         to pin it against the data rows; see results doc)
%             wS        - switch-row weight [5]
%             relTrunc  - SVD truncation (rel sigma_max) in the GN step [1e-8]
%             maxIter   - GN iterations [15]
%             blendTol  - relative misfit below which a node trusts the sweep [0.05]
%             hCS       - complex-step size [1e-20]
%             betaFit   - ('sweep' only) fit beta via fzero on q(beta)-1 [true]
%             betaTolX  - ('sweep' only) fzero rel tolerance [1e-3]
%             verbose   - print progress [false]
%
% OUTPUTS:
%   Z    - seed unknown vector [(8+17k)x1]
%   prob - IFS problem struct (rendezvous terminal; see ifs_residual)
%   meta - struct: k, tauSwitch, uArc, beta/beta0 (scale actually used / dual
%          fit), seedResNorm, signAgree, qSw, growth/growthTau (backward
%          amplification diagnostic), dualSrc, nSweeps, kThrottle (raw
%          throttle-crossing count), and for 'smooth': fitResHist, sFit,
%          misfit [1xnN] (per-node relative misfit), fracSweep (fraction of
%          seed costate nodes taken from the sweep vs dual fallback)
%
% REFERENCES:
%   [1] PLAN_OF_ATTACK_2.md, Rung A (this folder); RESULTS_RUNG01_RUNG2.md.
%   [2] CONSULT_GPT56_response.md Q4 (sweep rationale; its accuracy prediction
%       ignored the backward amplification this smoother variant repairs).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'method'),    opts.method    = 'smooth';  end
if ~isfield(opts,'tauParam'),  opts.tauParam  = 'sigmoid'; end
if ~isfield(opts,'odeRelTol'), opts.odeRelTol = 1e-11;     end
if ~isfield(opts,'odeAbsTol'), opts.odeAbsTol = 1e-13;     end
if ~isfield(opts,'nData'),     opts.nData     = 256;       end
if ~isfield(opts,'fitScale'),  opts.fitScale  = false;     end
if ~isfield(opts,'wS'),        opts.wS        = 5;         end
if ~isfield(opts,'relTrunc'),  opts.relTrunc  = 1e-8;      end
if ~isfield(opts,'maxIter'),   opts.maxIter   = 15;        end
if ~isfield(opts,'blendTol'),  opts.blendTol  = 0.05;      end
if ~isfield(opts,'hCS'),       opts.hCS       = 1e-20;     end
if ~isfield(opts,'betaFit'),   opts.betaFit   = true;      end
if ~isfield(opts,'betaTolX'),  opts.betaTolX  = 1e-3;      end
if ~isfield(opts,'verbose'),   opts.verbose   = false;     end

% ---- constants, node grid, direct states/throttle, dual costates -----------
[~, sd, info] = sms_seed_duals(matFile, 40, 1e-4, 'd');
tauN = info.tauN;  X = info.X;  U = info.U;  nN = numel(tauN);  N = nN - 1;
Tmax = sd.Tmax;  c = sd.c;  muStar = sd.muStar;  pSund = sd.pSund;
beta0 = info.beta;
lamDual = info.Y16(9:16, :);                 % dual-map costates, beta0 scale

% ---- terminal-BC duals ------------------------------------------------------
% lam_g layout of casadi_minfuel_sundman (verified on f1120/f1140):
%   defects 8N | unit nN | bounds 24*nN | init BC 8 | terminal BC 7 (rv 6, tf 1)
D = load(matFile);
expected = 8*N + nN + 24*nN + 8 + 7;
if isfield(D.out,'lamAll') && numel(D.out.lamAll) == expected
    L = D.out.lamAll(:);
    g_rv = L(end-6:end-1);  g_tf = L(end);
    dualSrc = 'lamAll';
else
    g_rv = -D.out.lamDef(1:6, end);  g_tf = -D.out.lamDef(8, end);
    dualSrc = 'lamDef';
end
gTerm = [g_rv(:); 0; g_tf];      % [8x1]; lam_m(tau_f)=0 EXACT

% ---- switch structure from the dual-map S (certified count) ----------------
Sn = 1 - sqrt(sum(lamDual(4:6,:).^2,1)).*c./X(7,:) - lamDual(7,:);
cr = find(diff(sign(Sn)) ~= 0);
tauSw = zeros(1, numel(cr));
for q = 1:numel(cr)
    kk = cr(q);
    tauSw(q) = tauN(kk) + (0 - Sn(kk))*(tauN(kk+1) - tauN(kk))/(Sn(kk+1) - Sn(kk));
end
k = numel(tauSw);
edges = [tauN(1), tauSw, tauN(end)];
uArc = zeros(1, k+1);
for a = 1:k+1
    uArc(a) = double(interp1(tauN, Sn, 0.5*(edges(a) + edges(a+1)), 'linear') < 0);
end
kThrottle = nnz(diff(sign(U(4,:) - 0.5)) ~= 0);      % cross-check only

% ---- frozen-state interpolant + sweep machinery -----------------------------
ppX = pchip(tauN, X);       % 8-comp state, C^1 at kinks; ppval -> [8 x nq]
odeO = odeset('RelTol', opts.odeRelTol, 'AbsTol', opts.odeAbsTol);
nSweeps = 0;

    function xq = xAt(tt)
        % frozen direct state at times tt -> [8 x numel(tt)]
        xq = ppval(ppX, tt(:).');
    end

    function dlam = adjointRhs(tt, lam, u)
        % costate half of ifs_eom at the frozen direct state
        Yf = [xAt(tt); lam];
        dY = ifs_eom([], Yf, Tmax, c, muStar, pSund, u);
        dlam = dY(9:16);
    end

    function sols = sweep(lamEnd)
        % backward arc-by-arc integration tau_f -> tau_0 from lam(tau_f)=lamEnd
        nSweeps = nSweeps + 1;
        sols = cell(1, k+1);
        for aa = k+1:-1:1
            sols{aa} = ode113(@(tt,y) adjointRhs(tt, y, uArc(aa)), ...
                              [edges(aa+1), edges(aa)], lamEnd, odeO);
            lamEnd = sols{aa}.y(:, end);              % value at arc start
        end
    end

    function lam = lamAt(sols, tt)
        % sample the sweep at times tt (row) -> [8 x numel(tt)]
        lam = zeros(8, numel(tt), 'like', sols{end}.y);
        for qq = 1:numel(tt)
            aa = find(tt(qq) >= edges(1:end-1) & tt(qq) <= edges(2:end), 1);
            lam(:, qq) = deval(sols{aa}, tt(qq));
        end
    end

    function [qm, qs] = qSwitch(sols)
        % switching condition at the switch times: want q = 1 (S = 0)
        lamS = lamAt(sols, tauSw);
        xS   = xAt(tauSw);
        qs   = sqrt(sum(lamS(4:6,:).^2, 1)).*c./xS(7,:) + lamS(7,:);
        qm   = mean(qs);
    end

    % --- smoother ('smooth' method) helpers; zIdx/Ddual/sc/tauD are set in the
    %     'smooth' case before these are first called ---------------------------
    function lamEnd = lamEndOf(zv)
        lamEnd = zeros(8,1,'like',zv);
        lamEnd(zIdx) = zv(1:7);               % lam_m(tau_f) stays 0
    end

    function [R, sols] = fitRes(zv)
        % scaled residual: [data misfit rows; switch rows]
        sols = sweep(lamEndOf(zv));
        lamD = lamAt(sols, tauD);                       % [8 x nd]
        Rd   = (lamD - zv(8)*Ddual)./sc;                % component-scaled
        [~, qs] = qSwitch(sols);
        R = [Rd(:); opts.wS*(qs(:) - 1)];
    end

meta = struct('k', k, 'tauSwitch', tauSw, 'uArc', uArc, 'beta0', beta0, ...
              'dualSrc', dualSrc, 'kThrottle', kThrottle, 'method', opts.method);

switch opts.method
% =============================================================================
case 'sweep'   % pure backward sweep (diagnostic; falsified as a seed at 1.12x)
% =============================================================================
    beta = beta0;
    if opts.betaFit && k > 0
        f = @(b) qSwitch(sweep(b*gTerm)) - 1;
        f0 = f(beta0);
        if opts.verbose, fprintf('beta0=%.6e  q-1=%+.3e\n', beta0, f0); end
        if f0 > 0, bLo = beta0/1.5; bHi = beta0; fLo = f(bLo); fHi = f0;
        else,      bLo = beta0; bHi = beta0*1.5; fLo = f0; fHi = f(bHi);
        end
        nExp = 0;
        while fLo*fHi > 0 && nExp < 8
            if f0 > 0, bLo = bLo/1.5; fLo = f(bLo);
            else,      bHi = bHi*1.5; fHi = f(bHi);
            end
            nExp = nExp + 1;
        end
        if fLo*fHi < 0
            beta = fzero(f, [bLo, bHi], optimset('TolX', opts.betaTolX*beta0));
        else
            warning('ifs_seed_adjoint:betaBracket', ...
                    'no sign change bracketing q(beta)=1; keeping beta0=%.4e', beta0);
        end
    end
    sols = sweep(beta*gTerm);
    lamN = real(lamAt(sols, tauN));
    meta.beta = beta;

% =============================================================================
case 'smooth'  % adjoint smoother: fit lam_f to the whole dual history
% =============================================================================
    % unknowns z = [lam_f(1:6); lam_f(8); s]; lam_m(tau_f) = 0 pinned exact
    zIdx = [1:6, 8];                          % lam_f components in z(1:7)
    idxD = unique(round(linspace(1, nN, opts.nData)));
    Ddual = lamDual(:, idxD)/beta0;           % dual history, unit (beta=1) scale
    sc = max(median(abs(lamDual), 2), 1e-3*max(median(abs(lamDual), 2)));  % [8x1]
    tauD = tauN(idxD);

    z = [beta0*gTerm(zIdx); beta0];           % init: terminal duals, s = beta0

    [R, sols] = fitRes(z);
    rn = norm(R);
    fitResHist = rn;
    if opts.verbose, fprintf('smoother it 0: ||R||=%.6e\n', rn); end

    for it = 1:opts.maxIter
        % Jacobian: complex-step over z(1:7) (a sweep each), analytic s-column
        J = zeros(numel(R), 8);
        for p = 1:7
            zp = complex(z);  zp(p) = zp(p) + 1i*opts.hCS;
            Rp = fitRes(zp);
            J(:, p) = imag(Rp)/opts.hCS;
        end
        Jd = -(Ddual./sc);                                   % d/ds of data rows
        J(:, 8) = [Jd(:); zeros(k, 1)];
        if ~opts.fitScale, J(:, 8) = 0; end                  % s pinned at beta0

        % truncated-SVD GN step on the column-equilibrated system
        cs = max(sqrt(sum(J.^2, 1)), 1e-300);
        [Us, Ss, Vs] = svd(J./cs, 'econ');
        sv = diag(Ss);
        keep = sv >= opts.relTrunc*sv(1);
        dz = -(Vs(:, keep)*(diag(1./sv(keep))*(Us(:, keep).'*R)))./cs.';

        % backtracking line search
        alpha = 1;  ok = false;
        for ls = 1:12
            [Rt, solsT] = fitRes(z + alpha*dz);
            rt = norm(Rt);
            if rt < (1 - 1e-4*alpha)*rn, ok = true; break; end
            alpha = alpha/2;
        end
        if ~ok
            if opts.verbose, fprintf('smoother it %d: no descent, stop\n', it); end
            break
        end
        z = z + alpha*dz;  R = Rt;  sols = solsT;  rnPrev = rn;  rn = rt;
        fitResHist(end+1) = rn; %#ok<AGROW>
        if opts.verbose
            fprintf('smoother it %d: ||R||=%.6e  alpha=%.3g  kept %d/8 sv\n', ...
                    it, rn, alpha, nnz(keep));
        end
        if rn > (1 - 1e-3)*rnPrev, break; end               % plateau
    end

    % per-node relative misfit + blended node costates
    lamSw = real(lamAt(sols, tauN));                        % [8 x nN]
    sFit  = real(z(8));
    mis   = sqrt(mean(((lamSw - sFit*lamDual/beta0)./sc).^2, 1));  % [1 x nN]
    trust = mis < opts.blendTol;
    lamN  = lamSw;
    lamN(:, ~trust) = sFit*lamDual(:, ~trust)/beta0;
    meta.beta = sFit;  meta.sFit = sFit;
    meta.fitResHist = fitResHist;  meta.misfit = mis;
    meta.fracSweep = mean(trust);

otherwise
    error('ifs_seed_adjoint:method', 'unknown method %s', opts.method);
end

% ---- diagnostics -------------------------------------------------------------
SnSw = 1 - sqrt(sum(lamN(4:6,:).^2,1)).*c./X(7,:) - lamN(7,:);
meta.signAgree = 100*mean((SnSw < 0) == (Sn < 0));
lamMag = sqrt(sum(lamN.^2, 1));
[meta.growth, gi] = max(lamMag./max(lamMag(end), realmin));
meta.growthTau = tauN(gi);
meta.nCrossSweep = nnz(diff(sign(SnSw)) ~= 0);
meta.nSweeps = nSweeps;
xS = xAt(tauSw);
lamS = zeros(8, k);
for q = 1:k
    [~, qn] = min(abs(tauN - tauSw(q)));  lamS(:, q) = lamN(:, qn);
end
meta.qSw = sqrt(sum(lamS(4:6,:).^2, 1)).*c./xS(7,:) + lamS(7,:);

% ---- assemble Z + prob (mirrors ifs_seed 'full') ----------------------------
lam0 = lamN(:, 1);
N16  = [xAt(tauSw); interp1(tauN.', lamN.', tauSw(:), 'linear').'];  % [16 x k]
odeOpts = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
prob = struct('rv0', sd.rv0(:), 'm0', 1, 't0', 0, 'tau0', tauN(1), ...
    'Tmax', Tmax, 'c', c, 'muStar', muStar, 'pSund', pSund, ...
    'tauf', tauN(end), 'k', k, 'uArc', uArc, ...
    'termMode', 'rendezvous', 'rvf', sd.rvf(:), 'tf', sd.tf, ...
    'odeOpts', odeOpts, 'tauParam', opts.tauParam);
Z = ifs_pack(lam0, N16, ifs_gseed(tauSw(:), tauN(1), tauN(end), opts.tauParam));

meta.seedResNorm = norm(ifs_residual(Z, prob));
end
