function lamSeed = costate_seed_from_nlp(tMesh, X, U, Tmax, c, muStar)
% COSTATE_SEED_FROM_NLP  Reconstruct initial costates from a direct solution.
%
% The proper covector mapping: along the (known) converged NLP trajectory
% the costate equation is LINEAR, lambda_dot = -A(x(t),u(t))' * lambda, so
% lambda(t) = Psi(t) * lambda(0) with Psi the 7x7 costate transition
% matrix (one 49-ODE integration). PMP then pins lambda(0) up to scale:
%   - on every burn arc the primer direction must match the NLP thrust
%     direction: lambda_v(t_k) x alpha_k = 0 (linear in lambda(0));
%   - free final mass: e7' * Psi(tf) * lambda(0) = 0.
% Solve the stacked homogeneous system by smallest singular vector; fix
% the SIGN so thrust opposes lambda_v, and the SCALE so the min-fuel
% switching function S = 1 - ||lambda_v||c/m - lambda_m vanishes at the
% NLP's first throttle switch. This sidesteps fmincon's multiplier
% estimates entirely (at an lbfgs step-tolerance stall they are
% feasibility-accurate but optimality-loose).
%
% INPUTS:
%   tMesh  - mesh times of the NLP solution [(N+1)x1] (ND)
%   X      - NLP states [7x(N+1)]
%   U      - NLP controls [w(3); s] [4x(N+1)]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   lamSeed - reconstructed initial costates [7x1] (Lagrange-form gauge)
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4 (covector mapping principle).

tMesh = tMesh(:);
tf    = tMesh(end);

% --- costate transition matrix along the NLP arc ---------------------------
xInterp = @(t) interp1(tMesh, X.', t, 'pchip').';
uInterp = @(t) interp1(tMesh, U.', t, 'pchip').';
    function dY = psiDot(t, Y)
        [~, A] = lt_dynamics_throttle(xInterp(t), uInterp(t), Tmax, c, muStar);
        dY = reshape(-A.' * reshape(Y, 7, 7), 49, 1);
    end

% sample times: burn nodes (s > 0.9), thinned to ~120, plus tf
s       = U(4, :);
burnIdx = find(s > 0.9);
if isempty(burnIdx)
    error('costate_seed_from_nlp:noBurn', ...
          'no node has throttle s > 0.9; the primer-direction rows need a burn arc');
end
burnIdx = burnIdx(round(linspace(1, numel(burnIdx), min(120, numel(burnIdx)))));
tSamp   = unique([tMesh(burnIdx); tf]);

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
tSpan = unique([0; tSamp]);
[tPsi, YPsi] = ode113(@psiDot, tSpan, reshape(eye(7), 49, 1), opts);
[tPsi, keep] = unique(tPsi, 'stable');
YPsi = YPsi(keep, :);

% --- stack the homogeneous conditions --------------------------------------
rows = zeros(3*numel(burnIdx) + 1, 7);
rPtr = 0;
for kB = burnIdx(:).'
    tk    = tMesh(kB);
    PsiK  = reshape(interp1(tPsi, YPsi, tk, 'pchip').', 7, 7);
    wK    = U(1:3, kB);
    alphK = wK./sqrt(sum(wK.^2));
    crossM = [0 -alphK(3) alphK(2); alphK(3) 0 -alphK(1); -alphK(2) alphK(1) 0];
    rows(rPtr+1:rPtr+3, :) = -crossM * PsiK(4:6, :);   % lambda_v x alpha = 0
    rPtr = rPtr + 3;
end
PsiF = reshape(interp1(tPsi, YPsi, tf, 'pchip').', 7, 7);
rows(rPtr+1, :) = PsiF(7, :);                          % lambda_m(tf) = 0
rows = rows(1:rPtr+1, :);

[~, ~, V] = svd(rows, 'econ');
lam0 = V(:, end);

% --- sign: thrust must OPPOSE lambda_v on the burn --------------------------
k1    = burnIdx(1);
Psi1  = reshape(interp1(tPsi, YPsi, tMesh(k1), 'pchip').', 7, 7);
lamV1 = Psi1(4:6, :)*lam0;
w1    = U(1:3, k1);
if (-lamV1).'*w1 < 0
    lam0 = -lam0;
end

% --- scale: S = 0 at the first throttle switch ------------------------------
% (node before the crossing; if the arc never switches, anchor at the last
% burn node instead -- S = 0 there is the weakest defensible gauge)
swIdx = find(abs(diff(s > 0.5)), 1, 'first');
if isempty(swIdx), swIdx = burnIdx(end); end
tSw   = tMesh(min(swIdx, numel(tMesh)));
PsiS  = reshape(interp1(tPsi, YPsi, tSw, 'pchip').', 7, 7);
lamS  = PsiS*lam0;
mSw   = interp1(tMesh, X(7,:).', tSw, 'pchip');
gauge = sqrt(sum(lamS(4:6).^2))*c/mSw + lamS(7);
lamSeed = lam0./gauge;
end
