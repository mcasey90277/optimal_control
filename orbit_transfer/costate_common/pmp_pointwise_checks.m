function P = pmp_pointwise_checks(t, Y, Tnd, cnd, mu, opts)
%% Purpose:
%
%   First-order (Pontryagin) checks evaluated ON A FLIGHT, at sampled times,
%   complementing the shooting residual: a small residual says the pieces
%   MATCH each other; these say the pieces are the right ones.
%
%     N2  Hamiltonian: autonomous, free final time  =>  H = 1 + lambda.f = 0.
%     N4  transversality on the free mass: lambda_m(t_f) = 0.
%     N5  adjoint equations lambda' = -dH/dx, with dH/dx by central
%         differences in the STATE at fixed costate, at two step sizes
%         (Richardson-combined; their disagreement is reported as fdAgree).
%         Differencing the propagator's OUTPUT instead measures its sample
%         spacing, not the equations.
%     N6  the minimum principle, EXACTLY, for the control the propagator
%         APPLIED: the thrust acceleration is recovered as the difference
%         between the powered and the coasting vector field, and the gap
%         H(applied) - min_alpha H = (T/m)(lambda_v . alpha + |lambda_v|)
%         is evaluated directly. Sampling the unit sphere against the
%         analytic minimiser tests nothing (the minimiser is constructed).
%         Also the applied throttle (must be 1) and the WEAK throttle
%         condition Q_mt = |lambda_v|/m + lambda_m/c >= 0, under which u = 1
%         minimises H; the STRICT version min Q_mt > 0 is the sufficiency
%         hypothesis S2.
%
%% Inputs:
%
%  t                        [N x 1]                 sample times [ND]
%  Y                        [N x 14]                [r v m lambda] samples
%  Tnd, cnd, mu             scalars                 thrust, exhaust speed,
%                                                   mass ratio [ND]
%  opts                     struct (optional)
%   .nSample [120] interior sample count, .hRel [1e-6] FD step relative to
%   max(1, |y_j|), .rhs [@mintime_rhs_point] the vector field under test
%   (signature rhs(y, Tmax, c, mu)); the test beside this file injects a
%   wrong-sign field here to prove the gap opens
%
%% Outputs:
%
%  P                        struct                  .Hmax .lamMf .adjErr
%                                                   .fdAgree .dirGap
%                                                   .throttleErr .minQmt
%                                                   .tSample [n x 1]
%                                                   .nSample .Hval [N x 1]
%                                                   .Qmt [N x 1]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
nS = min(d('nSample', 120), max(numel(t) - 2, 1));  hRel = d('hRel', 1e-6);
rhs = d('rhs', @mintime_rhs_point);
t = t(:);  N = numel(t);

% N2 on every sample
Hval = zeros(N, 1);
for k = 1:N
    F = rhs(Y(k,:).', Tnd, cnd, mu);
    Hval(k) = 1 + Y(k, 8:14)*F(1:7);
end
P.Hval = Hval;  P.Hmax = max(abs(Hval));

% N4
P.lamMf = abs(Y(end, 14));

% weak throttle condition along the whole flight
lamV = Y(:, 11:13);  lamVmag = sqrt(sum(lamV.^2, 2));
P.Qmt = lamVmag./Y(:, 7) + Y(:, 14)/cnd;
P.minQmt = min(P.Qmt);

% N5, N6 on interior samples
kk = unique(round(linspace(2, N - 1, nS)));
P.tSample = t(kk);  P.nSample = numel(kk);
adjErr = 0;  fdAgree = 0;  dirGap = 0;  throttleErr = 0;
for k = kk
    yk = Y(k, :).';
    F  = rhs(yk, Tnd, cnd, mu);
    F0 = rhs(yk, 0,   cnd, mu);                            % coasting field
    % adjoint: dH/dx at two step sizes, Richardson-combined
    h1 = hRel*max(1, abs(yk(1:7)));
    g1 = dHdx(rhs, yk, h1,   Tnd, cnd, mu);
    g2 = dHdx(rhs, yk, h1/2, Tnd, cnd, mu);
    gR = (4*g2 - g1)/3;
    adjErr  = max(adjErr,  norm(F(8:14) + gR)/max(norm(gR), 1));
    fdAgree = max(fdAgree, norm(g1 - g2)/max(norm(gR), 1));
    % minimum principle with the APPLIED control
    m = yk(7);  lv = yk(11:13);
    aT = F(4:6) - F0(4:6);                                 % (T/m) u alpha
    alphaApplied = aT*m/Tnd;                               % u alpha
    throttleErr = max(throttleErr, abs(norm(alphaApplied) - 1));
    gap = (Tnd/m)*(lv.'*alphaApplied + norm(lv));          % >= 0, 0 at the minimiser
    dirGap = max(dirGap, gap);
end
P.adjErr = adjErr;  P.fdAgree = fdAgree;  P.dirGap = dirGap;  P.throttleErr = throttleErr;
end

function g = dHdx(rhs, yk, h, Tnd, cnd, mu)
% DHDX  Central-difference dH/dx in the state at fixed costate.
% INPUTS: rhs; yk [14x1]; h [7x1]; Tnd; cnd; mu.  OUTPUTS: g [7x1].
g = zeros(7, 1);  lam = yk(8:14);
for jj = 1:7
    yp = yk;  yp(jj) = yp(jj) + h(jj);
    ym = yk;  ym(jj) = ym(jj) - h(jj);
    Fp = rhs(yp, Tnd, cnd, mu);
    Fm = rhs(ym, Tnd, cnd, mu);
    g(jj) = (lam.'*Fp(1:7) - lam.'*Fm(1:7))/(2*h(jj));
end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
