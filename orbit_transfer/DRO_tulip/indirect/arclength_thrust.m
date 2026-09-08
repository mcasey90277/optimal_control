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
%                                                   .Tstop [0.05] N,
%                                                   .m0kg [150], .logFile ''
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
Tstop   = d('Tstop', 0.05);
m0kg    = d('m0kg', 150);
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
lStar = 389703.264829278;  tStar = 382981.289129055;
ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);

%% Pack the starting unknown vector, exactly as ms_bvp orders it ----------
K  = numel(seed0.tGrid) - 1;
ny = 14;
p0 = [seed0.Y(8:14,1); reshape(seed0.Y(:,2:K), [], 1); seed0.tf];
n  = numel(p0);

% residual factory: one call per thrust value
mkRes = @(TN) resHandle(rv0, rvf, seed0, ndT(TN), cnd, muStar);

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
           'sminX', [], 'sminAug', [], 'normR', [], 'ds', [], ...
           'p', {{}}, 'fold', [], 'Tmin', NaN);

tau = [];                                 % previous tangent, for orientation
for step = 0:nStep
    TN = tpar*sT;
    Rf = mkRes(TN);
    [R, Jp] = feval(Rf, pOf(x));
    Jx = Jp .* Dx(:)';                    % R_x = R_p * D
    Rt = fdRT(mkRes, pOf(x), TN, sT);     % R_t = R_T * sT (central difference)

    % --- diagnostics at this root -----------------------------------------
    sX = svd(Jx);      sminX  = sX(end);
    sA = svd([Jx Rt]); sminAug = sA(end);

    % --- tangent: [Jx Rt] tau = 0, ||tau|| = 1 ----------------------------
    v = -(Jx \ Rt);
    tauNew = [v; 1];  tauNew = tauNew/norm(tauNew);
    if isempty(tau)
        if tauNew(end) > 0, tauNew = -tauNew; end     % head toward LOWER thrust
    elseif tauNew'*tau < 0
        tauNew = -tauNew;                              % keep the arc oriented
    end
    tau = tauNew;

    A.T_N(end+1)   = TN;
    A.tf_nd(end+1) = x(end)*Dx(end);          % t_f is the LAST unknown
    A.tf_days(end+1) = A.tf_nd(end)*tStar/86400;
    A.tauT(end+1) = tau(end);  A.sminX(end+1) = sminX;
    A.sminAug(end+1) = sminAug;  A.normR(end+1) = norm(R,inf);
    A.ds(end+1) = ds;  A.p{end+1} = pOf(x);
    lg(['  step %3d: T = %7.3f mN  tf = %7.3f d  tau_T = %+.4f  ' ...
        'smin(R_x) = %.2e  smin([R_x R_T]) = %.2e  ||R|| = %.1e'], ...
        step, TN*1000, A.tf_days(end), tau(end), sminX, sminAug, norm(R,inf));

    if numel(A.tauT) > 1 && sign(A.tauT(end)) ~= sign(A.tauT(end-1))
        A.fold(end+1) = numel(A.tauT);
        lg('  *** TANGENT THRUST COMPONENT CHANGED SIGN -> FOLD at T = %.3f mN ***', TN*1000);
    end
    if TN < Tstop, lg('  reached Tstop = %.1f mN', Tstop*1000); break, end
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
            F  = [Ri; tau'*(w - wp)];
            if norm(F, inf) < nTol, ok = true; break, end
            JF = [Jxi Rti; tau'];
            dw = -(JF \ F);
            if ~all(isfinite(dw)), break, end
            w = w + dw;
        end
        if ok
            x = w(1:end-1);  tpar = w(end);
            ds = min(dsMax, ds*1.3);
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
function h = resHandle(rv0, rvf, seed0, Tnd, cnd, muStar)
% RESHANDLE  The production ms residual for one thrust value.
% INPUTS: rv0; rvf; seed0; Tnd; cnd; muStar.  OUTPUTS: h (fhandle p -> [R,J]).
sd = seed0;  sd.Y(1:7,1) = [rv0(:); 1];
[~, inf_] = ms_tfmin(rv0(:), rvf(:), sd, Tnd, cnd, muStar, ...
    struct('assembleOnly', true, 'tfLo', 0.05, 'tfHi', 20, 'pMax', 1e7));
h = inf_.residual;
end

function Rt = fdRT(mkRes, p, TN, sT)
% FDRT  Central difference of the residual in the thrust parameter, at
% FIXED unknowns, returned in the SCALED parameter (R_T * sT).
% INPUTS: mkRes; p; TN; sT.  OUTPUTS: Rt [n x 1].
h = 1e-6*max(TN, 1e-3);
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
