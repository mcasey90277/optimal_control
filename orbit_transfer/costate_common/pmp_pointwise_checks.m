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
%         APPLIED. The applied control is recovered from the field itself,
%         on BOTH rows it enters: the thrust acceleration as the difference
%         between the powered and the coasting vector field, b = u*alpha,
%         and the throttle again from the mass row, u_mass = -c F_m / T.
%         (Astra review 2026-09-11: the acceleration norm establishes only
%         the acceleration-side throttle; the field must consume mass at
%         the SAME throttle.) The gap against the FULL control minimum,
%           H(u, alpha) - min_{u', beta} H
%             = (T/m)(lambda_v . b + |b| |lambda_v|)      direction part
%             + T (max(Q, 0) - |b| Q),                     throttle part
%           Q = |lambda_v|/m + lambda_m/c,
%         is evaluated directly and is zero exactly at the minimiser;
%         both its MINIMUM and its MAXIMUM over the samples are kept, so
%         an over-unit thrust along the minimiser (a negative gap) cannot
%         hide behind a max that starts at zero. Sampling the unit sphere
%         against the analytic minimiser tests nothing (the minimiser is
%         constructed). The WEAK throttle condition Q >= 0 is reported; the
%         STRICT version min Q > 0 is the sufficiency hypothesis S2.
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
%                                                   .fdAgree
%                                                   .dirGap (max |direction
%                                                   part|) .fullGap (max
%                                                   |full gap|) .gapMin
%                                                   .gapMax (signed extremes
%                                                   of the full gap)
%                                                   .throttleErr (max of the
%                                                   two below)
%                                                   .throttleAccErr
%                                                   .throttleMassErr .minQmt
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
adjErr = 0;  fdAgree = 0;
dirGap = 0;  fullGap = 0;  gapMin = Inf;  gapMax = -Inf;
throttleAccErr = 0;  throttleMassErr = 0;
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
    % minimum principle with the APPLIED control, recovered on both rows
    m = yk(7);  lv = yk(11:13);  rho = norm(lv);  lamM = yk(14);
    aT = F(4:6) - F0(4:6);                                 % (T/m) u alpha
    b  = aT*m/Tnd;                                         % u alpha (not unit unless u = 1)
    u  = norm(b);
    uMass = -cnd*F(7)/Tnd;                                 % the throttle the MASS row burned at
    throttleAccErr  = max(throttleAccErr,  abs(u - 1));
    throttleMassErr = max(throttleMassErr, abs(uMass - 1));
    Q = rho/m + lamM/cnd;
    gDir = (Tnd/m)*(lv.'*b + u*rho);                       % direction part, 0 at -lv/rho
    gThr = Tnd*(max(Q, 0) - u*Q);                          % throttle part, 0 at u = 1 when Q >= 0
    gFull = gDir + gThr;
    dirGap  = max(dirGap,  abs(gDir));
    fullGap = max(fullGap, abs(gFull));
    gapMin  = min(gapMin, gFull);  gapMax = max(gapMax, gFull);
end
P.adjErr = adjErr;  P.fdAgree = fdAgree;
P.dirGap = dirGap;  P.fullGap = fullGap;  P.gapMin = gapMin;  P.gapMax = gapMax;
P.throttleAccErr = throttleAccErr;  P.throttleMassErr = throttleMassErr;
P.throttleErr = max(throttleAccErr, throttleMassErr);
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
