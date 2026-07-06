% VERIFY_CHECKPOINTS  Headless capture of every tutorial checkpoint number.
%
%   Runs each checkpoint snippet exactly as the tutorial will present it and
%   prints the true values (scientific notation), so the tutorial's expected
%   outputs are verified rather than plausible. Safe under matlab -batch.
%
% INPUTS:  (none -- script)
% OUTPUTS: (none -- prints checkpoint values)

here = fileparts(mfilename('fullpath'));
cd(here);
rng(0);

mu = 1;
fprintf('===== CHECKPOINT A1: two_body_accel =====\n');
g = two_body_accel([1;0], mu);
fprintf('g([1;0])       = [% .6f % .6f]  (expect [-1 0])\n', g);
g2 = two_body_accel([0;2], mu);
fprintf('g([0;2])       = [% .6f % .6f]  (expect [0 -0.25])\n', g2);

fprintf('\n===== CHECKPOINT A2: gravity_gradient =====\n');
r_test = [0.8; -0.6];                          % ||r|| = 1
G = gravity_gradient(r_test, mu);
fprintf('G([0.8;-0.6])  = [% .6f % .6f; % .6f % .6f]\n', G(1,1), G(1,2), G(2,1), G(2,2));
fprintf('symmetry ||G-G''|| = %.3e\n', norm(G - G'));
fprintf('trace(G)*r^3/mu  = %.10f   (expect 1: trace G = mu/r^3)\n', trace(G)*1^3/mu);
% finite-difference check of G = dg/dr
h = 1e-7; Gfd = zeros(2);
for ci = 1:2
    e = zeros(2,1); e(ci) = h;
    Gfd(:,ci) = (two_body_accel(r_test + e, mu) - two_body_accel(r_test - e, mu)) / (2*h);
end
fprintf('FD Jacobian match ||G - G_fd|| = %.3e\n', norm(G - Gfd));

fprintf('\n===== CHECKPOINT B1: Hamiltonian conservation =====\n');
x0  = [1; 0; 0; 1];
lam_test = [0.01; -0.02; 0.03; 0.01];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
solB = ode45(@(t,z) ocp_dynamics(t,z,mu), [0 4.4], [x0; lam_test], opts);
tq = linspace(0, 4.4, 401);
zq = deval(solB, tq);
Hq = zeros(1, numel(tq));
for k = 1:numel(tq)
    z = zq(:,k); u = -z(7:8);
    f = [z(3:4); two_body_accel(z(1:2), mu) + u];
    Hq(k) = 0.5*(u'*u) + z(5:8)'*f;
end
fprintf('H(0) = %.10f   max|H(t)-H(0)| = %.3e\n', Hq(1), max(abs(Hq - Hq(1))));

fprintf('\n===== CHECKPOINT B2: coast residual (lam0 = 0) =====\n');
rf = 1.5; xf = [-rf; 0; 0; -sqrt(mu/rf)]; tf = 4.4;
psi0 = shoot_residual(zeros(4,1), x0, xf, tf, mu);
fprintf('psi(0) = [% .6f % .6f % .6f % .6f]   ||psi|| = %.6f\n', psi0, norm(psi0));

fprintf('\n===== CHECKPOINT C: indirect solve =====\n');
tic;
[lam0, sol, J_ind, psi_norm] = solve_indirect(x0, xf, tf, mu, zeros(4,1));
tC = toc;
fprintf('lam0   = [% .10f % .10f % .10f % .10f]\n', lam0);
fprintf('||psi|| = %.3e   J_ind = %.10f   time %.2f s\n', psi_norm, J_ind, tC);
umax = max(vecnorm(-deval(sol, linspace(0,tf,1001), 7:8)));
fprintf('max ||u(t)|| = %.6f  (sanity: small, near-Hohmann)\n', umax);

fprintf('\n===== CHECKPOINT D: collocation solve =====\n');
N = 40;
tic;
[Xc, Uc, tc, J_col, out] = collocation_transfer(x0, xf, tf, mu, N);
tD = toc;
fprintf('N = %d  exitflag = %d  max|defect| = %.3e  J_col = %.10f  time %.1f s\n', ...
        N, out.exitflag, out.max_defect, J_col, tD);
fprintf('J_col - J_ind = %+.3e\n', J_col - J_ind);

fprintf('\n===== CHECKPOINT E: primer + cross-method =====\n');
[err_primer, ~] = primer_check(sol, tf, mu);
fprintf('max ||p_ode - (-lam_v)|| = %.3e\n', err_primer);
z_nodes = deval(sol, tc);
fprintf('max ||u_dir - u_ind|| at nodes = %.3e\n', max(vecnorm(Uc + z_nodes(7:8,:))));
fprintf('max ||x_dir - x_ind|| at nodes = %.3e\n', max(vecnorm(Xc - z_nodes(1:4,:))));

% N-sweep for the convergence report-back
fprintf('\n===== CHECKPOINT D2: collocation convergence in N =====\n');
for Ns = [10 20 40 80]
    [~, ~, ~, Jn, on] = collocation_transfer(x0, xf, tf, mu, Ns);
    fprintf('N = %3d   J = %.10f   J - J_ind = %+.3e   exitflag %d\n', ...
            Ns, Jn, Jn - J_ind, on.exitflag);
end
