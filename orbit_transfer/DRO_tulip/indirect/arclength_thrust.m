function A = arclength_thrust(rv0, rvf, seed0, TN0, cnd, muStar, opts)
%% Purpose:
%
%   Pseudo-arclength continuation of the minimum-time multiple-shooting
%   solution in (X, T): follow the solution CURVE through a turning point
%   instead of stepping in thrust, which cannot pass one.
%
%   Why this exists (FINDINGS section 32): thrust continuation of the
%   DRO -> tulip family stalls at 75.0 mN with |lam_0| diverging (5.4 -> 45.9)
%   and cond(J) reaching 1e12, while the conjugate test still PASSES and
%   dim S stays 1 -- a singular Jacobian, not a loss of optimality. A direct
%   solve at 70 mN then converges, but at 26.44 d against 15.15 d at
%   75.5 mN. Either the short family FOLDS near 75 mN (and 26.44 d is the
%   answer) or it continues and a ~16 d solution exists unfound. Those give
%   OPPOSITE verdicts on the cislunar abstract, so the branch must be
%   followed, not guessed.
%
%   A simple fold shows: tangent thrust-component -> 0 then CHANGES SIGN,
%   R_X loses exactly one rank, and the augmented [R_X R_T] stays full rank.
%
%  ASSUMPTIONS / NOTES:
%
% • Everything is done in SCALED coordinates. The unknowns mix costates
%   (O(10-50)), junction states (O(1)) and t_f (O(3)); an unscaled Euclidean
%   arclength would be dominated by whichever block happens to be largest,
%   and would change meaning when K changes (GPT-6 Astra, 2026-09-08).
% • R_T is a central finite difference at FIXED unknowns -- differentiating
%   short independent segment propagations, not a re-optimized solution.
% • The engine's residual is reached through ms_bvp's assembleOnly mode, so
%   this shares the production residual and Jacobian rather than copying it.
%
%% Inputs:
%
%  rv0, rvf                 [6 x 1]                 departure / arrival states
%
%  seed0                    struct                  CONVERGED start: .tf,
%                                                   .tGrid [1 x K+1],
%                                                   .Y [14 x K+1]
%
%  TN0                      double                  starting thrust, NEWTONS
%
%  cnd, muStar              double                  ND exhaust speed, mass ratio
%
%  opts                     struct (optional)       .ds [0.05] initial arclength
%                                                   step, .dsMin [1e-4],
%                                                   .dsMax [0.4], .nStep [120],
%                                                   .newtonTol [1e-9],
%                                                   .newtonMax [12],
%                                                   .newtonTarget [4] Newton
%                                                   iterations the step
%                                                   controller aims for,
%                                                   .Tstop [0.05] N,
%                                                   .m0kg [150], .logFile '',
%                                                   .direction [-1] initial
%                                                   thrust sense (+1 = up),
%                                                   .binding [@ms_tfmin]
%                                                   or @ms_tfmin_hom (rho
%                                                   free on the sphere)
%
%% Outputs:
%
%  A                        struct                  .T_N .tf_nd .tf_days
%                                                   .tauT (tangent thrust
%                                                   component) .sminX
%                                                   .sminAug .normR .p
%                                                   .fold (indices where
%                                                   tauT changes sign)
%                                                   .Tmin (deepest thrust
%                                                   reached)
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 7, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
ds      = d('ds', 0.05);
dsMin   = d('dsMin', 1e-4);
dsMax   = d('dsMax', 0.4);
nStep   = d('nStep', 120);
nTol    = d('newtonTol', 1e-9);
nMax    = d('newtonMax', 12);
nTarget = d('newtonTarget', 4);
direction = d('direction', -1);
Tstop   = d('Tstop', 0.05);
Tmax_   = d('Tmax', 1.0);
m0kg    = d('m0kg', 150);
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
lStar = 389703.264829278;  tStar = 382981.289129055;
% THE shared propulsion conversion; c arrives as an argument here
ndp = nd_propulsion([], [], m0kg, lStar, tStar);
ndT = ndp.ndT;

%% The binding: normal chart (ms_tfmin, rho = 1) or the homogeneous one
%  (ms_tfmin_hom, rho free on the sphere). The unknown vector is packed by
%  the ENGINE and read back from its assembleOnly info, so this driver never
%  assumes where t_f or the extras sit.
binding = d('binding', @ms_tfmin);
K  = numel(seed0.tGrid) - 1;
ny = 14;
[~, inf0] = resHandle(binding, rv0, rvf, seed0, ndT(TN0), cnd, muStar);
p0  = inf0.p;  n = numel(p0);
ctf = inf0.ctf;                                % index of t_f in p
nX  = inf0.nExtra;                             % extras (rho) after t_f
lamIdx = [1:7, reshape(bsxfun(@plus,(8:14)',(0:(K-2))*ny),1,[]) + 7];

% residual factory: one call per thrust value
mkRes = @(TN) resHandle(binding, rv0, rvf, seed0, ndT(TN), cnd, muStar);

%% Scaling ----------------------------------------------------------------
Dx = max(abs(p0), 1e-2);                 % per-unknown scale, floored
sT = TN0;                                % thrust scale
pOf = @(x) x .* Dx;
x   = p0 ./ Dx;       tpar = TN0/sT;

[R, Jp] = feval(mkRes(TN0), p0);
lg('arclength start: T = %.2f mN, tf = %.4f ND (%.3f d), ||R|| = %.2e, n = %d', ...
   TN0*1000, seed0.tf, seed0.tf*tStar/86400, norm(R,inf), n);
assert(norm(R, inf) < 1e-6, 'start point is not a root (||R|| = %.2e)', norm(R,inf));

A = struct('T_N', [], 'tf_nd', [], 'tf_days', [], 'tauT', [], ...
           'sminX', [], 'sminAug', [], 'normR', [], 'dsUsed', [], ...
           'dsNext', [], 'nNewton', [], 'tanResid', [], 'lam0', [], ...
           'pMaxAbs', [], 'fracCostate', [], 'rho', [], 'p', {{}}, 'fold', [], 'Tmin', NaN);

tau = [];  dsUsed = 0;  nNewtonLast = 0;                                 % previous tangent, for orientation
for step = 0:nStep
    TN = tpar*sT;
    Rf = mkRes(TN);
    [R, Jp] = feval(Rf, pOf(x));
    Jx = Jp .* Dx(:)';                    % R_x = R_p * D
    Rt = fdRT(mkRes, pOf(x), TN, sT);     % R_t = R_T * sT (central difference)

    % --- diagnostics at this root -----------------------------------------
    sX = svd(Jx);      sminX  = sX(end);
    sA = svd([Jx Rt]); sminAug = sA(end);

    % --- tangent: right null vector of the FULL augmented matrix ---------
    % NOT v = -Jx\Rt: that chart is exactly singular AT a fold, which is the
    % one place the tangent matters (Astra review 2026-09-08). The full SVD
    % (not 'econ' -- a wide matrix's extra right-null column is dropped by
    % the economy form) gives the structural null vector directly.
    [~, ~, VA] = svd([Jx Rt]);
    tauNew = VA(:, end);
    tanResid = norm([Jx Rt]*tauNew) / max(norm([Jx Rt], 'fro'), realmin);
    if isempty(tau)
        % initial orientation: opts.direction = -1 (default) heads toward
        % LOWER thrust, +1 toward higher -- both ends of a branch are needed
        % to map it (Astra: two ends beat 500 one-sided steps)
        if sign(tauNew(end)) ~= direction, tauNew = -tauNew; end
    elseif tauNew'*tau < 0
        tauNew = -tauNew;                              % keep the arc oriented
    end
    tau = tauNew;

    pNow = pOf(x);
    A.T_N(end+1)   = TN;
    A.tf_nd(end+1) = pNow(ctf);               % t_f, wherever the engine put it
    if nX > 0, A.rho(end+1) = pNow(end); else, A.rho(end+1) = NaN; end
    A.tf_days(end+1) = A.tf_nd(end)*tStar/86400;
    A.tauT(end+1) = tau(end);  A.sminX(end+1) = sminX;
    A.sminAug(end+1) = sminAug;  A.normR(end+1) = norm(R,inf);
    A.dsUsed(end+1) = dsUsed;                 % the step that was ACCEPTED,
    A.dsNext(end+1) = ds;                     % not the next proposal
    A.nNewton(end+1) = nNewtonLast;
    A.tanResid(end+1) = tanResid;
    A.lam0(end+1) = norm(pNow(1:7));          % is the chart running away?
    A.pMaxAbs(end+1) = max(abs(pNow));
    % how much of the last step was costates vs state/time? If the arclength
    % is being eaten by multiplier growth, the fold reading is wrong.
    if numel(A.p) >= 1
        % in SCALED coordinates -- the metric the arclength actually uses.
        % Unscaled, O(50) costates dominate O(1) states by magnitude alone
        % and the fraction reads 1.000 whatever the geometry is doing.
        dx = (pNow - A.p{end}) ./ Dx;
        iL = lamIdx(lamIdx < ctf);
        A.fracCostate(end+1) = norm(dx(iL))/max(norm(dx(1:ctf-1)), realmin);
    else
        A.fracCostate(end+1) = NaN;
    end
    A.p{end+1} = pNow;
    lg(['  step %3d: T = %10.6f mN  tf = %10.6f d  tau_T = %+.4e  ' ...
        'sminX = %.2e  sminAug = %.2e  |R| = %.1e  dsUsed = %.2e  ' ...
        'nNw = %2d  |lam0| = %9.3f  fracCostate = %.3f  tanRes = %.1e  rho = %.6f'], ...
        step, TN*1000, A.tf_days(end), tau(end), sminX, sminAug, norm(R,inf), ...
        dsUsed, nNewtonLast, A.lam0(end), A.fracCostate(end), tanResid, A.rho(end));

    if numel(A.tauT) > 1 && sign(A.tauT(end)) ~= sign(A.tauT(end-1))
        A.fold(end+1) = numel(A.tauT);
        lg('  *** TANGENT THRUST COMPONENT CHANGED SIGN -> FOLD at T = %.3f mN ***', TN*1000);
    end
    if TN < Tstop, lg('  reached Tstop = %.1f mN', Tstop*1000); break, end
    if TN > Tmax_,  lg('  reached Tmax = %.1f mN',  Tmax_*1000); break, end
    if step == nStep, break, end

    % --- predictor / corrector -------------------------------------------
    ok = false;
    while ~ok && ds >= dsMin
        w0 = [x; tpar];  wp = w0 + ds*tau;  w = wp;
        for it = 1:nMax
            Ti = w(end)*sT;
            if Ti <= 0, break, end
            [Ri, Jpi] = feval(mkRes(Ti), pOf(w(1:end-1)));
            if ~all(isfinite(Ri)) || norm(Ri,inf) > 1e6, break, end
            Jxi = Jpi .* Dx(:)';
            Rti = fdRT(mkRes, pOf(w(1:end-1)), Ti, sT);
            if ~all(isfinite(Jxi(:))) || ~all(isfinite(Rti)), break, end
            F  = [Ri; tau'*(w - wp)];
            if norm(F, inf) < nTol, ok = true; nNewtonLast = it - 1; break, end
            JF = [Jxi Rti; tau'];
            dw = -(JF \ F);
            if ~all(isfinite(dw)), break, end
            w = w + dw;
            if it == nMax
                % test convergence AFTER the final update: a point that
                % converges on the last iterate was being reported as a
                % failure and forcing a needless halving.
                Tl = w(end)*sT;
                if Tl > 0
                    Rl = feval(mkRes(Tl), pOf(w(1:end-1)));
                    if all(isfinite(Rl)) && ...
                       norm([Rl; tau'*(w - wp)], inf) < nTol
                        ok = true;  nNewtonLast = nMax;
                    end
                end
            end
        end
        if ok
            x = w(1:end-1);  tpar = w(end);
            dsUsed = ds;
            % scale by the effort actually spent, not by mere success: a
            % fixed 1.3x inflation after every step produces overshoot/halve
            % cycles that look like an arclength collapse.
            fac = sqrt(nTarget/max(nNewtonLast, 1));
            ds  = min(dsMax, ds*min(max(fac, 0.5), 2.0));
        else
            ds = ds/2;
            if ds < dsMin, lg('  arclength stalled: ds < dsMin at T = %.3f mN', TN*1000); end
        end
    end
    if ~ok, break, end
end
A.Tmin = min(A.T_N);
lg('arclength done: %d points, thrust %.3f -> %.3f mN, %d fold(s)', ...
   numel(A.T_N), A.T_N(1)*1000, A.Tmin*1000, numel(A.fold));
end

% ------------------------------------------------------------------------
function [h, inf_] = resHandle(binding, rv0, rvf, seed0, Tnd, cnd, muStar)
% RESHANDLE  The production ms residual for one thrust value, through the
% chosen binding.  INPUTS: binding (fhandle); rv0; rvf; seed0; Tnd; cnd;
% muStar.  OUTPUTS: h (fhandle p -> [R,J]); inf_ (assembleOnly info).
sd = seed0;  sd.Y(1:7,1) = [rv0(:); 1];
[~, inf_] = binding(rv0(:), rvf(:), sd, Tnd, cnd, muStar, ...
    struct('assembleOnly', true, 'handleOnly', true, 'tfLo', 0.05, 'tfHi', 20, 'pMax', 1e7));
h = inf_.residual;
end

function Rt = fdRT(mkRes, p, TN, sT)
% FDRT  Central difference of the residual in the thrust parameter, at
% FIXED unknowns, returned in the SCALED parameter (R_T * sT).
% INPUTS: mkRes; p; TN; sT.  OUTPUTS: Rt [n x 1].
% Step-size study 2026-09-08 (rt_study): the derivative is on a flat
% plateau for h_rel in [1e-2, 1e-4] (successive changes ~1e-10 relative),
% and drifts to ~1e-7 by h_rel = 1e-6 -- the OLD choice sat at the noisy
% end. sigma_min([R_x R_T]) was identical to 4 digits across the whole
% range, so the augmented-regularity evidence stands.
h = 1e-3*max(TN, 1e-3);
Rp = feval(mkRes(TN + h), p);
Rm = feval(mkRes(TN - h), p);
Rt = (Rp - Rm)/(2*h) * sT;
end

function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function logmsg(f, s)
% LOGMSG  Append to a log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
