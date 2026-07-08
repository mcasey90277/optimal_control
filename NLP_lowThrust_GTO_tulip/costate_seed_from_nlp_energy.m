function lamSeed = costate_seed_from_nlp_energy(tMesh, X, U, Tmax, c, muStar)
% COSTATE_SEED_FROM_NLP_ENERGY  Reconstruct min-ENERGY costates from a
% direct solution (covector mapping).
%
% Along the converged NLP arc the costate equation is LINEAR,
% lambda_dot = -A(x,u)' lambda, so lambda(t) = Psi(t) lambda(0) with Psi
% the 7x7 costate transition matrix (one 49-ODE integration). PMP pins
% lambda(0):
%   - DIRECTION: on every thrust-on node the primer must oppose the NLP
%     thrust direction, lambda_v(t_k) x alpha_k = 0 (homogeneous, linear in
%     lambda(0)); plus free final mass lambda_m(tf) = 0. Solve by smallest
%     singular vector -> lambda(0) up to sign and scale.
%   - SIGN: thrust must OPPOSE lambda_v.
%   - SCALE: min-energy has NO bang-bang switch to anchor against, but it
%     has something better -- the SMOOTH stationarity u = S_e holds at
%     every interior-throttle node:
%         s(t_a) = S_e = Tmax( ||lambda_v(t_a)||/m + lambda_m(t_a)/c ).
%     Pick an interior node (0.3 < s < 0.9) and scale lambda(0) so this
%     holds exactly. (For min-fuel the analogous anchor is S = 0 at a
%     switch; the energy version is better conditioned -- any interior
%     node works, and s there is O(1).)
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
%   lamSeed - reconstructed initial costates [7x1]
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4 (covector mapping principle).
%   [2] Caillau, Gergaud, Noailles, JOTA 2003 (min-energy transfer).

tMesh = tMesh(:);
tf    = tMesh(end);

xInterp = @(t) interp1(tMesh, X.', t, 'pchip').';
uInterp = @(t) interp1(tMesh, U.', t, 'pchip').';
    function dY = psiDot(t, Y)
        [~, A] = lt_dynamics_throttle(xInterp(t), uInterp(t), Tmax, c, muStar);
        dY = reshape(-A.' * reshape(Y, 7, 7), 49, 1);
    end

% thrust-on nodes (s > 0.3) for the direction rows, thinned to ~150
s      = U(4, :);
thrIdx = find(s > 0.3);
if isempty(thrIdx)
    error('costate_seed_from_nlp_energy:noThrust', ...
          'no node has throttle s > 0.3; need a thrust arc for primer rows');
end
thrIdx = thrIdx(round(linspace(1, numel(thrIdx), min(150, numel(thrIdx)))));
tSamp  = unique([tMesh(thrIdx); tf]);

opts   = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
tSpan  = unique([0; tSamp]);
[tPsi, YPsi] = ode113(@psiDot, tSpan, reshape(eye(7), 49, 1), opts);
[tPsi, keep] = unique(tPsi, 'stable');
YPsi = YPsi(keep, :);

% --- homogeneous system: direction rows + terminal mass row ----------------
rows = zeros(3*numel(thrIdx) + 1, 7);
rPtr = 0;
for kT = thrIdx(:).'
    tk    = tMesh(kT);
    PsiK  = reshape(interp1(tPsi, YPsi, tk, 'pchip').', 7, 7);
    wK    = U(1:3, kT);
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

% --- sign: thrust must OPPOSE lambda_v on a thrust node --------------------
k1    = thrIdx(1);
Psi1  = reshape(interp1(tPsi, YPsi, tMesh(k1), 'pchip').', 7, 7);
lamV1 = Psi1(4:6, :)*lam0;
w1    = U(1:3, k1);
if (-lamV1).'*w1 < 0
    lam0 = -lam0;
end

% --- scale: enforce u = S_e = s at an interior-throttle node ---------------
interior = find(s > 0.3 & s < 0.9);
if isempty(interior)                     % fully saturated arc: use max-s node
    [~, ia] = max(s);  ta = tMesh(ia);  sTgt = s(ia);
else
    ia = interior(round(numel(interior)/2));  ta = tMesh(ia);  sTgt = s(ia);
end
PsiA  = reshape(interp1(tPsi, YPsi, ta, 'pchip').', 7, 7);
lamA  = PsiA*lam0;
mA    = interp1(tMesh, X(7,:).', ta, 'pchip');
SeHat = Tmax*(sqrt(sum(lamA(4:6).^2))/mA + lamA(7)/c);
lamSeed = lam0*(sTgt/SeHat);
end
