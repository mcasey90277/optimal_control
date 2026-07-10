% DEMO_GAUSS_SUM_SPLIT  Curvature-split Gaussian-sum update vs single-Gaussian EKF.
%
% One nonlinear measurement update z = h(x) + noise on the saddle model, from an
% elongated prior that straddles the curved direction. We compare:
%   - TRUTH:  the exact Bayesian posterior on a dense grid (reference);
%   - EKF:    one Gaussian, linearize h at the prior mean, standard update;
%   - SPLIT:  split the prior along the max-curvature direction (demo_sr1_curvature
%             item 3) into K moment-matched children, EKF-update each, reweight by
%             measurement likelihood, recombine -> Gaussian mixture posterior.
% The split shrinks each child's spread in the curved direction, so its
% linearization error (~ curvature x spread^2) drops, and the mixture mean tracks
% the true posterior mean far better than the single Gaussian.
%
% Produces gauss_sum_split.png. REFERENCES: Alspach & Sorenson (1972), Gaussian sum.

clear; clc; rng(1);

% ---- problem setup (curvature-DOMINATED: weak direct channels, precise saddle) ----
% If the position channels (h_1=x1, h_2=x2) are precise the linear part pins the
% state and nonlinearity is irrelevant. Here we make them weak (large R) and the
% curved channel h_3=(x1^2-x2^2)/2 precise -- so the likelihood is a curved ridge
% that a single Gaussian cannot represent, and the split has something to do.
Ksplit = 5;                          % children along the max-curvature direction
x_true = [1.2; 0.7];
mu0    = [0.7; 0.7];                 % off-axis (x2~=0, non-degenerate linearization)
P0     = diag([0.7^2, 0.5^2]);
R      = diag([0.5^2, 0.5^2, 0.04^2]);   % moderately weak position, precise curvature
z      = saddle_h(x_true) + chol(R,'lower')*randn(3,1);

% ---- TRUTH: dense-grid posterior p(x|z) ~ N(x;mu0,P0) * exp(-1/2 ||z-h(x)||_R^-1) ----
g  = linspace(-1.5, 3.0, 401);
[X1,X2] = meshgrid(g,g);
logpost = zeros(size(X1));
iP0 = inv(P0); iR = inv(R);
for a = 1:numel(X1)
    xv = [X1(a); X2(a)];
    r  = z - saddle_h(xv);
    logpost(a) = -0.5*((xv-mu0)'*iP0*(xv-mu0)) - 0.5*(r'*iR*r);
end
post = exp(logpost - max(logpost(:)));  post = post / sum(post(:));
mean_true = [sum(X1(:).*post(:)); sum(X2(:).*post(:))];

% ---- EKF: single Gaussian, linearize at mu0 ----
[zhat, J] = saddle_h(mu0);
S  = J*P0*J' + R;
K  = P0*J'/S;
mu_ekf = mu0 + K*(z - zhat);
P_ekf  = (eye(2) - K*J)*P0;

% ---- SPLIT: along max-curvature direction, K=3 moment-matched children ----
H3 = [1 0; 0 -1];  L = chol(P0,'lower');
best_e=-inf; ud=[1;0];
for th = linspace(0,pi,1801)
    u=[cos(th);sin(th)]; v=L*u; e=(0.5*(v'*H3*v))^2;
    if e>best_e, best_e=e; ud=u; end
end
d = L*ud; d = d/norm(d);                      % split direction in state space
sig_d = sqrt(d'*P0*d);                        % prior std along d
% symmetric, equally-weighted K-split of the 1-D marginal N(0,sig_d^2): means at
% K offsets in [-a,a], per-component std gam*sig_d, with gam chosen to MATCH the
% parent marginal variance (mean 0 is automatic by symmetry).
K = Ksplit; a = 1.225; w = ones(1,K)/K;
u_off = linspace(-a, a, K);                   % offsets in units of sig_d
gam2 = 1 - mean(u_off.^2);                     % variance match: mean(o^2)+gam^2=1
assert(gam2 > 0, 'split too aggressive: gam^2<=0');
offs = u_off*sig_d;  dP = (gam2 - 1)*sig_d^2;  % variance change along d (reduce)
comp_mu = cell(1,K); comp_P = cell(1,K);
for c = 1:K
    comp_mu{c} = mu0 + offs(c)*d;
    comp_P{c}  = P0 + dP*(d*d');               % shrink covariance along d only
end
% sanity: the split must reproduce the prior (mean & covariance)
mmix = zeros(2,1); for c=1:K, mmix=mmix+w(c)*comp_mu{c}; end
Pmix = zeros(2); for c=1:K, dm=comp_mu{c}-mmix; Pmix=Pmix+w(c)*(comp_P{c}+dm*dm'); end
fprintf('split (K=%d) preserves prior:  ||mean-mu0||=%.2e  ||cov-P0||=%.2e\n', ...
        K, norm(mmix-mu0), norm(Pmix-P0,'fro'));

% EKF-update each child; reweight by measurement likelihood; recombine
wq = zeros(1,K); qmu = cell(1,K); qP = cell(1,K);
for c = 1:K
    [zc, Jc] = saddle_h(comp_mu{c});
    Sc = Jc*comp_P{c}*Jc' + R;
    Kc = comp_P{c}*Jc'/Sc;
    qmu{c} = comp_mu{c} + Kc*(z - zc);
    qP{c}  = (eye(2)-Kc*Jc)*comp_P{c};
    wq(c)  = w(c) * exp(-0.5*(z-zc)'*(Sc\(z-zc))) / sqrt(det(2*pi*Sc));
end
wq = wq / sum(wq);
mu_split = zeros(2,1); for c=1:K, mu_split=mu_split+wq(c)*qmu{c}; end

% ---- evaluate both approximations on the grid and score by KL(truth||approx) ----
% Mean-distance is a poor metric for the curved (non-Gaussian) posterior a nonlinear
% measurement produces; KL divergence from the true posterior is the right measure of
% how well each approximation represents it.
gauss2 = @(m,P) reshape(exp(-0.5*sum(([X1(:) X2(:)]'-m).*(P\([X1(:) X2(:)]'-m)),1))' ...
                        /(2*pi*sqrt(det(P))), size(X1));
g_ekf = gauss2(mu_ekf, P_ekf);
g_mix = zeros(size(X1));
for c=1:K, g_mix = g_mix + wq(c)*gauss2(qmu{c}, qP{c}); end
nrm = @(g) g/sum(g(:));
pe = nrm(g_ekf); pm = nrm(g_mix); pt = post;            % all normalized on the grid
fl = 1e-300;
KL = @(p,q) sum(p(:).*log((p(:)+fl)./(q(:)+fl)));        % KL(truth || approx)
KL_ekf = KL(pt, pe);  KL_mix = KL(pt, pm);

fprintf('\nhow well each represents the true posterior (KL(truth||approx), lower=better):\n');
fprintf('   single-Gaussian EKF : KL = %.4f\n', KL_ekf);
fprintf('   curvature-split (K=%d): KL = %.4f\n', K, KL_mix);
fprintf('   KL reduction factor  : %.2fx\n', KL_ekf/max(KL_mix,eps));
fprintf('   (means, secondary -- unreliable for curved posteriors: true [%+0.3f %+0.3f], EKF [%+0.3f %+0.3f], split [%+0.3f %+0.3f])\n', ...
        mean_true(1), mean_true(2), mu_ekf(1), mu_ekf(2), mu_split(1), mu_split(2));

% ---- figure ----
fig = figure('Visible','off','Position',[100 100 760 680]);
contour(X1,X2,post,12,'LineColor',[.6 .6 .6]); hold on; axis equal; grid on;
hT  = plot(x_true(1),x_true(2),'k*','MarkerSize',12,'LineWidth',1.5);
hMt = plot(mean_true(1),mean_true(2),'ks','MarkerSize',10,'LineWidth',1.5);
ell = @(m,P,col,lw) plot_ellipse_local(m,P,col,lw);
hE  = ell(mu_ekf, P_ekf, [0.85 0.2 0.2], 2);                 % EKF 1-sigma (red)
hC  = ell(qmu{1}, qP{1}, [0.1 0.45 0.85], 1.2);              % first child (blue, for legend)
for c=2:K, ell(qmu{c}, qP{c}, [0.1 0.45 0.85], 1.2); end     % remaining children
hSm = plot(mu_split(1),mu_split(2),'b^','MarkerSize',9,'LineWidth',1.5);
hEm = plot(mu_ekf(1),mu_ekf(2),'rv','MarkerSize',9,'LineWidth',1.5);
legend([hC hT hMt hE hSm hEm], ...
       {'split children','x_{true}','true mean','EKF 1\sigma','split mean','EKF mean'}, ...
       'Location','northwest','FontSize',9);
xlabel('x_1'); ylabel('x_2');
title(sprintf('Curvature-split (KL=%.2f) vs single-Gaussian EKF (KL=%.2f), saddle measurement', KL_mix, KL_ekf));
xlim([-1.5 3.0]); ylim([-1.5 2.0]);
exportgraphics(fig, 'gauss_sum_split.png', 'Resolution', 130);
fprintf('\nwrote gauss_sum_split.png\n');

% ---- local: 1-sigma ellipse (returns the line handle) ----
function hh = plot_ellipse_local(m, P, col, lw)
    t = linspace(0,2*pi,80); c = [cos(t); sin(t)];
    L = chol(P,'lower'); e = L*c + m;
    hh = plot(e(1,:), e(2,:), '-', 'Color', col, 'LineWidth', lw);
end
