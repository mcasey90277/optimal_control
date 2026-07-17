% TEST_ENERGY_STAGE  Coplanar eps=1 (energy) solve at tf = 1.5*tfmin, fixed L_f (gate G2).
%
% First call is LOOSE (warmTight defaults false -- genuine move off the
% propagated seed), maxIter=3000. Task-8 lesson: cold/loose solves on this
% problem routinely hit Maximum_Iterations_Exceeded at maxIter=3000 and then
% converge with one or two warm-start continuation rounds (warmTight=true,
% re-solving from the last iterate's X,U). If the single call does not reach
% the gate but the defect has come down substantially, apply that pattern
% here (up to 3 rounds; stall rule: <1 decade defect improvement -> stop).
% The eps=1 energy problem is convex-ish in the control so it is expected to
% behave at least as well as run_mintime's mintime stage did.
%
% CONTINGENCY (brief step 3, used only if continuation also fails): re-target
% the terminal at the seed's own arrival longitude si.Larr (zero topology
% gap), then walk L_f to the law value in 2-3 fixed-longitude warm-started
% continuation steps (the paper's own c_Lf device).
%
% Gate G2 stays fixed regardless of path taken: success && maxDefect<1e-8 &&
% |tf - tfTarget| < 1e-6.
%
% REFERENCES: [1] task-9-brief.md Step 2-3. [2] run_mintime.m (continuation
%   pattern this test reuses). [3] DESIGN.md sec 4.
CONTINGENCY_MAX_ROUNDS = 3;

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

fprintf('ENERGY STAGE: first (loose) solve @ eps=1, tf=%.4f, Lf=%.4f...\n', tf, Lf);
o = casadi_lt_2body(sg, X0, U0, tauf0, term, struct('par',p,'mode','fixedtf', ...
      'eps',1,'tfTarget',tf,'rv0',rv0,'maxIter',3000,'printLevel',3));
fprintf('  first: status=%s defect=%.3e tf=%.6f (target %.6f)\n', ...
        o.ipoptStatus, o.maxDefect, o.tf, tf);

round_ = 0;
usedContinuation = false;
if ~isGood(o) && o.maxDefect < 1e-2
    usedContinuation = true;
    while ~isGood(o) && round_ < CONTINGENCY_MAX_ROUNDS
        round_ = round_ + 1;
        prevDefect = o.maxDefect;
        oNew = casadi_lt_2body(sg, o.X, o.U, tauf0, term, struct('par',p,'mode','fixedtf', ...
              'eps',1,'tfTarget',tf,'rv0',rv0,'maxIter',3000,'warmTight',true,'printLevel',3));
        fprintf('  continuation round %d: defect %.3e -> %.3e, status=%s\n', ...
                round_, prevDefect, oNew.maxDefect, oNew.ipoptStatus);
        if ~isGood(oNew)
            decadeImprove = log10(max(prevDefect, realmin)) - log10(max(oNew.maxDefect, realmin));
            if decadeImprove < 1
                fprintf('  continuation stalled at round %d: %.3e -> %.3e (%.2f decades)\n', ...
                        round_, prevDefect, oNew.maxDefect, decadeImprove);
                o = oNew;
                break;
            end
        end
        o = oNew;
    end
end

usedContingency = false;
if ~isGood(o)
    % CONTINGENCY: retarget terminal at seed's own arrival longitude (zero
    % topology gap), then walk L_f to the law value in fixed-longitude steps.
    usedContingency = true;
    fprintf('ENERGY STAGE: CONTINGENCY -- retargeting terminal at seed arrival Larr=%.4f\n', si.Larr);
    termSeed = geo_terminal('fixed', p, si.Larr);
    oc = casadi_lt_2body(sg, X0, U0, tauf0, termSeed, struct('par',p,'mode','fixedtf', ...
          'eps',1,'tfTarget',tf,'rv0',rv0,'maxIter',3000,'printLevel',3));
    fprintf('  contingency seed-Larr solve: status=%s defect=%.3e\n', oc.ipoptStatus, oc.maxDefect);
    LfSteps = si.Larr + (0:1/3:1) * (Lf - si.Larr);   % 4 steps incl. endpoints
    LfSteps = LfSteps(2:end);                         % skip the seed-Larr point already solved
    for kL = 1:numel(LfSteps)
        termK = geo_terminal('fixed', p, LfSteps(kL));
        oc = casadi_lt_2body(sg, oc.X, oc.U, tauf0, termK, struct('par',p,'mode','fixedtf', ...
              'eps',1,'tfTarget',tf,'rv0',rv0,'maxIter',3000,'warmTight',true,'printLevel',3));
        fprintf('  contingency Lf-walk step %d/%d (Lf=%.4f): status=%s defect=%.3e\n', ...
                kL, numel(LfSteps), LfSteps(kL), oc.ipoptStatus, oc.maxDefect);
    end
    o = oc;
end

assert(o.success, 'energy solve failed: %s', o.ipoptStatus);
assert(o.maxDefect < 1e-8, 'energy defect %.2e', o.maxDefect);
assert(abs(o.tf - tf) < 1e-6, 'tf pin violated');
save(fullfile('results','energy_M0_coplanar.mat'), 'o', 'sg', 'tauf0', 'tf', 'Lf', 'rv0', 'si', ...
     'round_', 'usedContinuation', 'usedContingency');
fprintf('test_energy_stage: ALL PASS (mf=%.2f kg, edge=%.1f%%, continuationRounds=%d, contingency=%d)\n', ...
        o.m_f_kg, 100*o.edge, round_, usedContingency);
