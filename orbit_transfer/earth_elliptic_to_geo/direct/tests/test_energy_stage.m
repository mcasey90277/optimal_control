% TEST_ENERGY_STAGE  Coplanar eps=1 (energy) solve at tf = 1.5*tfmin, fixed L_f (gate G2).
%
% Single-path solve: loose call (warmTight=false, maxIter=3000) on eps=1 energy
% problem, gate assertions unchanged.
% If this solve ever fails on a rerun (max-iter), apply the warm-continuation
% pattern from run_mintime.m — deliberately not coded here to keep the gate test single-path.
%
% Gate G2: success && maxDefect<1e-8 && |tf - tfTarget| < 1e-6.
%
% REFERENCES: [1] task-9-brief.md Step 2-3. [2] run_mintime.m. [3] process/DESIGN.md sec 4.

root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;

p   = kepler_lt_params(10, 1500, 2000);
P0  = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0, 0, pi, p.mu);
rv0 = [r0; v0];
mt  = run_mintime(10, 0, 600);
tf  = 1.5 * mt.tfmin;
Lf  = pi + (1.12*1.5 + 0.09) * mt.dL_mt;      % paper law R2 (c_Lf ~ 1.77)
[sg, X0, U0, tauf0, si] = seed_2body(p, rv0, ...
      struct('sbar', 1/1.5, 'tDur', tf, 'N', 600, 'targetLf', Lf));
term = geo_terminal('fixed', p, Lf);

isGood = @(o) o.success && o.maxDefect < 1e-8;

fprintf('ENERGY STAGE: loose solve @ eps=1, tf=%.4f, Lf=%.4f...\n', tf, Lf);
o = casadi_lt_2body(sg, X0, U0, tauf0, term, struct('par',p,'mode','fixedtf', ...
      'eps',1,'tfTarget',tf,'rv0',rv0,'maxIter',3000,'printLevel',3));
fprintf('  status=%s defect=%.3e tf=%.6f (target %.6f)\n', ...
        o.ipoptStatus, o.maxDefect, o.tf, tf);

assert(o.success, 'energy solve failed: %s', o.ipoptStatus);
assert(o.maxDefect < 1e-8, 'energy defect %.2e', o.maxDefect);
assert(abs(o.tf - tf) < 1e-6, 'tf pin violated');
save(fullfile(module_root(),'results','energy_M0_coplanar.mat'), 'o', 'sg', 'tauf0', 'tf', 'Lf', 'rv0', 'si');
fprintf('test_energy_stage: ALL PASS (mf=%.2f kg, edge=%.1f%%)\n', o.m_f_kg, 100*o.edge);
