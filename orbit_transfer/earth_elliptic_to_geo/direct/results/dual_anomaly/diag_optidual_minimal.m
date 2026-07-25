function diag_optidual_minimal()
% DIAG_OPTIDUAL_MINIMAL  Minimal reproduction of the opti.dual() sign behavior
% that caused the Campaign-B PMP primer anomaly.
%
% A 3-variable NLP with one VECTOR equality constraint, solved three ways that
% differ only in how the same mathematical constraint is WRITTEN. For each we
% compare, entry by entry:
%     opti.dual(con)            (the convenience accessor -- what the solvers
%                                used to call, once per constraint)
%     opti.lam_g(rows)          (the raw low-level multiplier vector -- what
%                                they call now)
% and the ground-truth multiplier from hand-solved KKT stationarity
%     grad_f + A' * lam = 0.
%
% The full-scale finding this isolates (results/dual_anomaly/diag_rawdual.m, on
% the certified 10 N MEE row): identical magnitudes, entry-wise SIGN
% differences, and the raw lam_g is the one that satisfies stationarity
% (||grad_x L||_inf = 1.5e-14) and yields a 0.000 deg primer.
%
% INPUTS:  none
% OUTPUTS: none (prints a comparison table per constraint-writing style)
%
% REFERENCES:
%   [1] diag_t1_beta.m, diag_rawdual.m (the full-scale evidence).
%   [2] process/DESIGN_dual_map.md sec 3 T1 (H2-new: extraction convention).
import casadi.*
c     = [0.5; 0.5; 0.5];
tgt   = [1; 2; 3];
forms = {'x - c == 0', 'c - x == 0', 'x == c'};

fprintf('\nMinimal opti.dual vs opti.lam_g probe\n');
fprintf('objective  min sum((x-[1;2;3]).^2),  constraint pins x = [0.5;0.5;0.5]\n');
fprintf('hand KKT   grad_f = 2*(x-tgt) = [%g %g %g];  lam s.t. grad_f + A''lam = 0\n\n', ...
    2*(c - tgt));

for q = 1:numel(forms)
    opti = casadi.Opti();
    x = opti.variable(3,1);
    opti.minimize(sum((x - tgt).^2));
    r0 = size(opti.g,1) + 1;
    switch q
        case 1, con = (x - c) == 0;
        case 2, con = (c - x) == 0;
        case 3, con = x == c;
    end
    opti.subject_to(con);
    rows = r0:size(opti.g,1);
    opti.solver('ipopt', struct('print_time', false), struct('print_level', 0));
    sol = opti.solve();

    lamRaw  = full(sol.value(opti.lam_g));
    lamRaw  = lamRaw(rows);
    lamDual = full(sol.value(opti.dual(con)));
    xv      = full(sol.value(x));

    % ground truth from stationarity: grad_f + A'*lam = 0, A = d(con)/dx
    Fa   = Function('Fa', {opti.x}, {jacobian(opti.g, opti.x)});
    Aall = full(Fa(full(sol.value(opti.x))));
    A    = Aall(rows, :);
    gradf = 2*(xv - tgt);
    lamTrue = -(A.') \ gradf;

    fprintf('--- constraint written as:  %s\n', forms{q});
    fprintf('    %-12s %-12s %-12s %-12s\n', 'lam_g(rows)', 'opti.dual', 'KKT truth', 'dual/lam_g');
    for r = 1:3
        fprintf('    %-12.6g %-12.6g %-12.6g %-+12.3f\n', ...
            lamRaw(r), lamDual(r), lamTrue(r), lamDual(r)/lamRaw(r));
    end
    fprintf('    stationarity resid: with lam_g %.2e | with opti.dual %.2e\n\n', ...
        norm(gradf + A.'*lamRaw, inf), norm(gradf + A.'*lamDual, inf));
end
end
