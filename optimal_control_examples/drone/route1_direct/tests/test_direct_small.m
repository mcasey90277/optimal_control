function ok = test_direct_small()
%% Purpose:
%
%   The M = 4 paper circle through the solver and the gates:
%     1. hard mode: verify_drones_direct says ok, every gate's number is in
%        range, J is above the rigorous two-pair bound 0.2003202 (relax to the
%        two antipodal head-on pairs, each the two-agent closed form) and
%        within 0.5 % of J_LB = 0.200. Finding: the hard optimum 0.2006465 is
%        above the paper's OCFormer 0.2003 (which sits at the pair bound with
%        no room for cross-pair avoidance) and SympOCNet's 0.2001 is infeasible
%        under the stated clearance; the paper numbers are reported, not gated;
%     2. the mirror solve (opposite swirl) agrees in cost to 1e-6;
%     3. barrier mode: solves; out.Unodes equals barrier_terms at the
%        solution nodes to 1e-12 (pins the CasADi copy of the barrier);
%        its fine-grid clearance is above l; its cost is above the hard
%        optimum (the barrier can only add);
%     4. a failed solve (max_iter = 2) is NOT verified, and .why names the
%        IPOPT status;
%     5. the sphere layout (not mirror-symmetric) skips G5 with a note.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));

[x0, xT, geom] = drone_layouts('paper', 4);
oH = solve_drones_direct(x0, xT, geom, struct('mode', 'hard'));
rep = verify_drones_direct(oH, struct('mirror', true, 'mesh', false));
ok = chk(ok, rep.ok, sprintf('hard M = 4: verify ok (%s)', tern(isempty(rep.why), 'no reason recorded', rep.why)));
% rigorous lower bound: relax to the two antipodal head-on pairs, each the two-agent closed form (test_direct_oracle)
D = 1;  R = 2*geom.r;  T = 10;
ell = 2*sqrt(D^2 - R^2) + R*(pi - 2*acos(R/D));
pairBound = 2*ell^2/(4*T);
ok = chk(ok, oH.J >= pairBound - 1e-9, sprintf('hard M = 4: J = %.7f >= two-pair rigorous bound %.7f', oH.J, pairBound));
ok = chk(ok, (oH.J - oH.JLB)/oH.JLB <= 5e-3, sprintf('hard M = 4: gap to J_LB %.3f %% (<= 0.5 %%)', 100*(oH.J - oH.JLB)/oH.JLB));
fprintf('  NOTE  paper OCFormer 0.2003 (SympOCNet 0.2001): ours - OCFormer = %+.2e; 0.2001 lies below the rigorous pair bound %.7f\n', oH.J - 0.2003, pairBound);
ok = chk(ok, rep.G1.minHfine >= -1e-2*2*geom.r, sprintf('G1 fine-grid clearance %.2e (bar %.1e)', rep.G1.minHfine, -1e-2*2*geom.r));
ok = chk(ok, rep.G3.maxDefect < 1e-9, sprintf('G3 defects %.1e', rep.G3.maxDefect));
ok = chk(ok, rep.G5.ran && rep.G5.relDiff < 1e-6, sprintf('G5 mirror: J = %.6f vs %.6f (%.1e rel)', oH.J, rep.G5.Jmirror, rep.G5.relDiff));

[h0, hT, gH] = drone_layouts('headon');
oM = solve_drones_direct(h0, hT, gH, struct('mode', 'hard', 'N', 50));
repM = verify_drones_direct(oM, struct('mirror', false, 'mesh', true));
ok = chk(ok, repM.G4.ran && repM.G4.relChange <= 1e-3, ...
         sprintf('G4 mesh (headon, N = 50 -> 100): ran = %d, J2N = %.7f, relChange = %.2e (<= 1e-3)', repM.G4.ran, repM.G4.J2N, repM.G4.relChange));

oB = solve_drones_direct(x0, xT, geom, struct('mode', 'barrier'));
ok = chk(ok, oB.success, sprintf('barrier M = 4: %s', oB.ipoptStatus));
Ucheck = zeros(1, oB.N + 1);
for k = 1:oB.N + 1
    Ucheck(k) = barrier_terms(oB.X(:, k), oB.pairs, geom);
end
ok = chk(ok, max(abs(Ucheck - oB.Unodes)) < 1e-12, sprintf('CasADi barrier = barrier_terms at the nodes (%.1e)', max(abs(Ucheck - oB.Unodes))));
repB = verify_drones_direct(oB, struct('mirror', false));
ok = chk(ok, repB.ok, sprintf('barrier M = 4: verify ok (%s)', tern(isempty(repB.why), '', repB.why)));
ok = chk(ok, repB.G1.minHfine > geom.l, sprintf('barrier: fine-grid clearance %.2e > l = %.0e', repB.G1.minHfine, geom.l));
ok = chk(ok, oB.J > oH.J, sprintf('barrier cost %.6f > hard cost %.6f', oB.J, oH.J));

oF = solve_drones_direct(x0, xT, geom, struct('mode', 'hard', 'maxIter', 2));
repF = verify_drones_direct(oF, struct('mirror', false));
ok = chk(ok, ~oF.success && ~repF.ok && contains(repF.why, oF.ipoptStatus), ...
         sprintf('a failed solve is not verified: %s', repF.why));

[s0, sT, gS] = drone_layouts('sphere', 6);
oS = solve_drones_direct(s0, sT, gS, struct('mode', 'hard', 'N', 60));
repS = verify_drones_direct(oS, struct('mirror', true));
ok = chk(ok, oS.success && repS.ok && ~repS.G5.ran && ~isempty(repS.G5.note), ...
         sprintf('sphere M = 6: ok, mirror skipped (%s), ball clearance %.2e', repS.G5.note, repS.G1.minHballFine));

if ok, fprintf('TEST_DIRECT_SMALL: ALL PASS\n'); else, fprintf('TEST_DIRECT_SMALL: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
