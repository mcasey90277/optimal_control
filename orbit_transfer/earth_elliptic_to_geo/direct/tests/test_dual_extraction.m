% TEST_DUAL_EXTRACTION  Regression guard for the Campaign-B dual-extraction bug
% (resolved 2026-07-25): defect duals must come from opti.lam_g by constraint-row
% range, NEVER from opti.dual().
%
% WHAT BROKE. opti.dual(con) returns the multiplier for CasADi's CANONICALIZED
% constraint orientation, not the orientation of the opti.g row it pairs with in
% grad_f + A'*lam = 0. The collocation defect canonicalizes row by row, so the
% duals came back with exact magnitudes and ~44-60% of signs inverted -- which
% read as a structured, eccentricity-correlated primer rotation and kept the
% first-order PMP gates failing (32.4 deg / 78.4%) across two campaigns.
%
% This file locks BOTH halves of the fix:
%   Test 1 -- the CasADi behavior itself, on a 3-variable NLP. If a future
%     CasADi release changes opti.dual's convention, this test tells us
%     (it asserts the CURRENT behavior, including the failure), so the
%     lesson documents reality rather than folklore.
%   Test 2 -- the property that actually matters, on a real (small) MEE solve:
%     the duals casadi_lt_mee hands back must satisfy Lagrangian stationarity,
%     measured by the tangential residual (I - beta*beta') dL/dbeta at burn
%     nodes. Tangential because the ||beta||^2 == 1 multiplier is purely radial
%     and legitimately absorbs the radial component.
%
% Test 2 is the one to trust: it would fail for ANY dual corruption, not just
% the sign-flip class this episode happened to produce.
%
% REFERENCES:
%   [1] process/DESIGN_dual_map.md (status banner: root cause + evidence).
%   [2] results/dual_anomaly/diag_optidual_minimal.m (the interactive version
%       of Test 1, with the printed comparison table).
%   [3] results/dual_anomaly/diag_t1_beta.m (the full-scale version of Test 2).
%   [4] process/LESSONS_DUAL_EXTRACTION.md (why it took so long to find).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
import casadi.*

%% Test 1 -- opti.dual follows canonicalized orientation; opti.lam_g does not.
c = [0.5; 0.5; 0.5];  tgt = [1; 2; 3];
resid = struct('lamg', zeros(1,3), 'dual', zeros(1,3));
for q = 1:3
    opti = casadi.Opti();
    x = opti.variable(3,1);
    opti.minimize(sum((x - tgt).^2));
    r0 = size(opti.g,1) + 1;
    switch q
        case 1, con = (x - c) == 0;      % variable-first
        case 2, con = (c - x) == 0;      % NEGATED form -- the trap
        case 3, con = x == c;            % lhs/rhs form
    end
    opti.subject_to(con);
    rows = r0:size(opti.g,1);
    opti.solver('ipopt', struct('print_time', false), struct('print_level', 0));
    sol = opti.solve();

    xv    = full(sol.value(x));
    Fa    = Function('Fa', {opti.x}, {jacobian(opti.g, opti.x)});
    A     = full(Fa(full(sol.value(opti.x))));
    A     = A(rows, :);
    gradf = 2*(xv - tgt);

    lamRaw  = full(sol.value(opti.lam_g));  lamRaw = lamRaw(rows);
    lamDual = full(sol.value(opti.dual(con)));
    resid.lamg(q) = norm(gradf + A.'*lamRaw,  inf);
    resid.dual(q) = norm(gradf + A.'*lamDual, inf);
end
assert(all(resid.lamg < 1e-8), ...
    ['Test1: opti.lam_g must satisfy grad_f + A''*lam = 0 for EVERY way of ' ...
     'writing the constraint (got %s) -- if this fails the fix''s premise is ' ...
     'wrong, not just its convention'], mat2str(resid.lamg, 3));
assert(resid.dual(1) < 1e-8 && resid.dual(3) < 1e-8, ...
    'Test1: opti.dual should still be correct for the variable-first forms');
assert(resid.dual(2) > 1e-3, ...
    ['Test1: opti.dual is expected to FAIL stationarity for the negated form ' ...
     '(c - x == 0); it did not (residual %.3e). CasADi may have changed its ' ...
     'dual convention -- re-read process/LESSONS_DUAL_EXTRACTION.md before ' ...
     'assuming the old advice still applies.'], resid.dual(2));
fprintf(['Test 1 (opti.dual canonicalization trap) PASSED -- lam_g resid %s | ' ...
         'dual resid %s\n'], mat2str(resid.lamg,3), mat2str(resid.dual,3));

%% Test 2 -- casadi_lt_mee's duals satisfy beta-stationarity at a CONVERGED point.
% Uses the certified 10 N row (same dependency as test_sosc_recover_10N.m): a
% warm re-solve there lands machine-tight in ~1 min, whereas a cold small solve
% does NOT reach a KKT point cheaply (a 75-node 3-rev eps=1 build stops at
% Maximum_Iterations_Exceeded / Infeasible_Problem_Detected), and stationarity
% residuals are meaningless at a non-converged iterate. Skips if the cache is
% absent -- results/*.mat are gitignored campaign caches.
row = fullfile(module_root(), 'results', 'MEE_M2_10N.mat');
if ~isfile(row)
    fprintf(['Test 2 SKIPPED -- %s absent (gitignored cache). Re-run after ' ...
             'the 10 N row is built.\n'], row);
    fprintf('test_dual_extraction: Test 1 PASS, Test 2 skipped\n');
    return
end
saved = sosc_load_row(row);
p     = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
out   = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, ...
    struct('par', p, 'mode', 'fixedtf', 'eps', 0, 'tfTarget', saved.tfTarget, ...
           'x0', saved.X(:,1), 'xf', saved.xf, 'maxIter', saved.maxIter, ...
           'warmTight', true, 'printLevel', 0, 'returnModel', true));
assert(out.success && out.maxDefect < 1e-8, ...
    ['Test2: warm re-solve of the 10 N row did not land machine-tight ' ...
     '(%s, defect %.3e) -- cannot judge stationarity'], ...
    out.ipoptStatus, out.maxDefect);

opti = out.model.opti;  sol = opti.debug;
xv   = full(sol.value(opti.x));
lamv = full(sol.value(opti.lam_g));
N1   = size(out.U, 2);
nX   = 7 * N1;
assert(max(abs(xv(nX + (1:4*N1)) - out.U(:))) < 1e-10, ...
    'Test2: Opti variable-layout assumption (X then U) failed');

Fk = Function('Fk', {opti.x, opti.lam_g}, ...
    {gradient(opti.f, opti.x), jacobian(opti.g, opti.x)});
[gfD, AD_] = Fk(xv, lamv);
gL = full(gfD) + sparse(AD_).' * lamv;          % s=+1 convention, asserted below
assert(norm(gL, inf) < 1e-6, ...
    ['Test2: full Lagrangian stationarity ||grad_x L||_inf = %.3e is not ' ...
     'small -- the recovered multipliers are not a KKT certificate'], norm(gL, inf));

burn = out.U(4,:) > 0.5;
assert(any(burn), 'Test2: no burn nodes to test beta-stationarity on');
tanRes = zeros(1, N1);
for k = 1:N1
    b = out.U(1:3,k) / norm(out.U(1:3,k));
    v = gL(nX + (k-1)*4 + (1:3));
    tanRes(k) = norm(v - (v.'*b)*b);
end
tanMed = median(tanRes(burn));
assert(tanMed < 1e-8, ...
    ['Test2: tangential dL/dbeta on burn nodes has median %.3e -- beta ' ...
     'stationarity is violated, which is the exact signature of corrupted ' ...
     'defect duals (see process/LESSONS_DUAL_EXTRACTION.md)'], tanMed);

% and the delivered out.lamDef must BE those raw multipliers, not a re-derivation
di = find(strcmp({out.model.creg.label}, 'defect'), 1);
lamDefRaw = reshape(lamv(out.model.creg(di).rows), 7, N1-1);
assert(max(abs(lamDefRaw(:) - out.lamDef(:))) < 1e-12, ...
    ['Test2: out.lamDef does not match opti.lam_g on the defect rows -- ' ...
     'someone reintroduced an opti.dual() extraction path']);
fprintf(['Test 2 (casadi_lt_mee duals satisfy beta-stationarity) PASSED -- ' ...
         '||grad_x L||_inf %.2e, tangential median %.2e\n'], norm(gL,inf), tanMed);

fprintf('test_dual_extraction: ALL PASS\n');
