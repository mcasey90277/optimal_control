% RUN_ORBIT_TRANSFER  Demo: planar min-energy transfer, three lenses.
%
%   Solves the canonical tutorial instance (circular r=1 to circular r=1.5,
%   near-Hohmann transfer time) by:
%     (1) the INDIRECT method  -- shooting on the costate BVP (fsolve), and
%     (2) the DIRECT method    -- trapezoidal collocation NLP (fmincon),
%   then (3) verifies both against Lawden's PRIMER equation p'' = Gp and
%   against each other. Prints every tutorial checkpoint quantity and draws
%   the transfer + control/primer profiles.
%
% INPUTS:  (none -- script)
% OUTPUTS: (none -- prints checkpoints, opens one figure)
%
% REFERENCES:
%   [1] myLatex/notes/cov_pmp_orbit_transfer.tex (companion theory note).

clear; clc;

% --- the tutorial instance (canonical units: mu = 1, r0 = 1) -------------
mu  = 1;
r0  = 1;            v0 = sqrt(mu/r0);            % circular speed at r0
rf  = 1.5;          vf = sqrt(mu/rf);            % circular speed at rf
x0  = [ r0; 0; 0;  v0];                          % depart (r0, 0), CCW
xf  = [-rf; 0; 0; -vf];                          % arrive (-rf, 0), CCW
tf  = 4.4;                                       % ~ Hohmann half-period 4.39

fprintf('=== Min-energy planar transfer: r=1 -> r=1.5, tf = %.2f TU ===\n\n', tf);

% --- (1) indirect: shooting ----------------------------------------------
tic;
[lam0, sol, J_ind, psi_norm] = solve_indirect(x0, xf, tf, mu, zeros(4,1));
t_ind = toc;
fprintf('[indirect]  lam0 = [% .10f  % .10f  % .10f  % .10f]\n', lam0);
fprintf('[indirect]  ||psi|| = %.3e   J = %.10f   (%.2f s)\n\n', psi_norm, J_ind, t_ind);

% --- (2) direct: trapezoidal collocation ---------------------------------
N = 40;
tic;
[Xc, Uc, tc, J_col, out] = collocation_transfer(x0, xf, tf, mu, N);
t_col = toc;
fprintf('[direct]    N = %d   exitflag = %d   max|defect| = %.3e\n', ...
        N, out.exitflag, out.max_defect);
fprintf('[direct]    J = %.10f   (J_dir - J_ind = %+.3e)   (%.2f s)\n\n', ...
        J_col, J_col - J_ind, t_col);

% --- (3) primer verification ---------------------------------------------
[err_primer, ~] = primer_check(sol, tf, mu);
fprintf('[primer]    max ||p_ode(t) + lam_v(t)|| = %.3e   (p''''=Gp verified)\n', err_primer);

% control agreement between the two methods, sampled at the nodes
z_nodes  = deval(sol, tc);
U_ind    = -z_nodes(7:8, :);
err_ctrl = max(vecnorm(Uc - U_ind));
fprintf('[cross]     max ||u_direct - u_indirect|| at nodes = %.3e\n', err_ctrl);

% state agreement at the nodes
err_state = max(vecnorm(Xc - z_nodes(1:4, :)));
fprintf('[cross]     max ||x_direct - x_indirect|| at nodes = %.3e\n\n', err_state);

% --- figure ---------------------------------------------------------------
tq = linspace(0, tf, 1001);
zq = deval(sol, tq);
th = linspace(0, 2*pi, 361);

figure('Name', 'Min-energy orbit transfer', 'Position', [80 80 1100 450]);

subplot(1,2,1); hold on; axis equal; grid on;
plot(r0*cos(th), r0*sin(th), 'k:');
plot(rf*cos(th), rf*sin(th), 'k:');
plot(zq(1,:), zq(2,:), 'b-', 'LineWidth', 1.5);
plot(Xc(1,:), Xc(2,:), 'ro', 'MarkerSize', 4);
plot(x0(1), x0(2), 'ks', 'MarkerFaceColor', 'g');
plot(xf(1), xf(2), 'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
xlabel('x_1'); ylabel('x_2');
title('Transfer: indirect (line) vs collocation (nodes)');
legend('r = 1', 'r = 1.5', 'indirect', 'collocation', 'depart', 'arrive', ...
       'Location', 'southoutside', 'NumColumns', 3);

subplot(1,2,2); hold on; grid on;
plot(tq, -zq(7,:), 'b-', 'LineWidth', 1.5);
plot(tq, -zq(8,:), 'r-', 'LineWidth', 1.5);
plot(tc, Uc(1,:), 'bo', 'MarkerSize', 4);
plot(tc, Uc(2,:), 'ro', 'MarkerSize', 4);
xlabel('t [TU]'); ylabel('control = primer');
title('u^*(t) = p(t): indirect (lines) vs collocation (markers)');
legend('p_1 = u_1 (ind)', 'p_2 = u_2 (ind)', 'u_1 (dir)', 'u_2 (dir)', ...
       'Location', 'best');
