% P0H_DIAG_FLOOR  Diagnose the P0h 1.6e-2 floor at 75 mN: residual breakdown,
% CS-vs-central-difference Jacobian cross-check, weakest singular direction.
%
% Decides between "CS-through-ode113 Jacobian quality is the rate limiter"
% (cure = Z0's exact variational STM) and "the floor lives in a specific
% residual component / near-null direction" (cure = formulation-level).

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

A = load(fullfile(here, 'results', 'p0h_gn_finish.mat'));
a = A.anchor;
Tmax = (a.Tmax_mN/25)*a.P.Tmax25;
resFun = @(lam) shoot_residual_energy(lam, a.tf, a.rv0, 1, a.rvf, Tmax, a.P.c, a.P.muStar);

[R, Jcs] = resFun(a.lam0);
lbl = {'rx','ry','rz','vx','vy','vz','lam_m'};
fprintf('R breakdown at the 75 mN floor (||R|| = %.3e):\n', norm(R));
for k = 1:7, fprintf('  %-5s %+.4e\n', lbl{k}, R(k)); end

Jfd = zeros(7);
for k = 1:7
    h = max(1e-6*abs(a.lam0(k)), 1e-8);
    ep = zeros(7,1);  ep(k) = h;
    Jfd(:,k) = (resFun(a.lam0+ep) - resFun(a.lam0-ep))/(2*h);
end
fprintf('J: CS vs central-diff rel error = %.3e   cond(Jcs) = %.2e  cond(Jfd) = %.2e\n', ...
        norm(Jcs-Jfd)/norm(Jfd), cond(Jcs), cond(Jfd));

[~, S, V] = svd(Jcs);  sv = diag(S);
fprintf('singular values: %.3e .. %.3e\n', sv(1), sv(end));
vn = abs(V(:,end));  [~, ix] = sort(vn, 'descend');
fprintf('weakest-direction mass: %s %.2f | %s %.2f | %s %.2f\n', ...
        lbl{ix(1)}, vn(ix(1))^2, lbl{ix(2)}, vn(ix(2))^2, lbl{ix(3)}, vn(ix(3))^2);

% Newton step size along the weakest direction vs the floor residual
dzFull = -Jcs\R;
fprintf('full Newton step norm: %.3e  (per-component max %.3e)\n', norm(dzFull), max(abs(dzFull)));
Rt = resFun(a.lam0 + dzFull);
fprintf('full Newton step: ||R|| %.3e -> %.3e\n', norm(R), norm(Rt));
