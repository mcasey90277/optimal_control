function ok = test_polish_rendezvous()
%% Purpose:
%
%   Verifies polish_rendezvous on the real schedule: (1) spacecraft 3's
%   chosen candidate polished without gates -- root found, pinned arrival
%   epoch, departure and slot states exact, flown miss inside the inner
%   rule, a cell-sized t_f change; (2) the same with gates on -- every
%   measured gate finite, flight / pointwise / hypothesis gates pass, the
%   conjugate verdict printed; (3a) the preflight departure-before-epoch
%   refusal; (3b) the AT-ROOT refusal on spacecraft 1's real candidate
%   (its cell departs at clock 0 and the root needs t_f 1.45 h longer:
%   measured tDep = -0.0137 ND); (4) a deliberate gate failure through
%   opts.gateOverride (fly, gates off); (5) a gate error through the
%   forceError seam (hyp); (6) a deliberate non-fly gate failure (dimS
%   requirement, gates on).
%
%   Spacecraft 3, not 1, carries checks 1, 2 and 4 because spacecraft 1
%   is refused by the tDep < 0 rule (measured 2026-11-10); the assertions
%   are the brief's, unchanged.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));

M = mission_inputs();
S = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, ...
                        M.slotPhase0, struct('propKgPerNd', M.propKgPerNd));
k = 3;
j = S.assign(k);
cand = struct('iD', S.iD(k), 'iA', S.iA(k), 'nRev', S.nRev(k), ...
              'tDep', S.tDep(k), 'tArr', S.tArr(k));
dayNd = M.tStar/86400;

% ---- 1. polish without gates ----
P = polish_rendezvous(M, k, j, cand, struct('gates', false));
fprintf(['  sc %d slot %d: ok %d (%s)  tf cell %.6f d -> %.6f d  nEval %d  ' ...
         'g %.2e  fly %.3g km %.3g m/s  %.1f s\n'], k, j, P.ok, P.reason, ...
        P.tfCell*dayNd, P.tf*dayNd, P.nEval, P.gResid, P.flight.flyKm, ...
        P.flight.flyVms, P.wallS);
assert(P.ok, 'polish not ok: %s', P.reason);
assert(abs(P.gResid) < 1e-9, 'g residual %.3g', P.gResid);
assert(P.flight.flyKm < 100 && P.flight.flyVms < 10, 'flown miss');
assert(abs((P.tDep + P.tf) - cand.tArr) < 1e-12, 'arrival epoch not pinned');
rv0Want = M.depState(mod(M.droPhase0(k) + P.tDep/M.tauD, 1));
assert(max(abs(P.rv0 - rv0Want)) < 1e-12, 'departure state');
rvfWant = M.arrState(mod(M.slotPhase0(j) + cand.tArr/M.tauT, 1));
assert(max(abs(P.rvf - rvfWant)) < 1e-12, 'slot state');
relChange = abs(P.tf - P.tfCell)/P.tfCell;
assert(relChange < 0.05, 'tf moved %.3g of the cell value', relChange);
assert(P.nEval <= 25, 'nEval %d', P.nEval);

% ---- 2. with gates ----
G = polish_rendezvous(M, k, j, cand, struct('gates', true));
g = G.gates;
pw = g.pw;  hy = g.hyp;
fprintf('  gates: resid %.2e  fly %.3g km %.3g m/s  (passFly %d)\n', ...
        g.resid, g.flyKm, g.flyVms, g.passFly);
fprintf(['  pointwise: Hmax %.2e  lamMf %.2e  adjErr %.2e  fullGap %.2e  ' ...
         'fieldGap %.2e  throttleErr %.2e  minQmt %.3g  (passPW %d)\n'], ...
        pw.Hmax, pw.lamMf, pw.adjErr, pw.fullGap, pw.fieldGap, ...
        pw.throttleErr, pw.minQmt, g.passPW);
fprintf('  hypothesis: minLamV %.3g  minQmt %.3g  dimS %d  (passHyp %d)\n', ...
        hy.minLamV, hy.minQmt, hy.dimS, g.passHyp);
fprintf('  conjugate: conjPass %g  covered %g  (%s)  passConj %d\n', ...
        g.conjPass, g.conjCovered, g.conjReason, g.passConj);
fprintf('  tfMin acceptance: accDz %.2e  (passAcc %d)\n', g.accDz, g.passAcc);
fprintf('  ok %d (%s)  wall %.1f s (root %.1f, cert %.1f)\n', G.ok, G.reason, ...
        G.wallS, G.wallRootS, G.wallCertS);
vals = [g.resid g.flyKm g.flyVms pw.Hmax pw.lamMf pw.adjErr pw.fullGap ...
        pw.fieldGap pw.throttleErr pw.minQmt hy.minLamV hy.minQmt hy.dimS g.accDz];
assert(all(isfinite(vals)), 'a measured gate value is not finite');
assert(g.passFly, 'flight gate');
assert(g.passPW, 'pointwise gate');
assert(g.passHyp, 'hypothesis gate');
assert(g.passAcc, 'tfMin acceptance gate');
assert(ismember(g.conjPass, [0 1]) || isnan(g.conjPass), 'conjPass value');
assert(ischar(g.conjReason) && ~isempty(g.conjReason), 'conjugate reason');

% ---- 3a. departure before epoch 0 (preflight) ----
bad = struct('iD', 1, 'iA', 1, 'nRev', 0, 'tDep', 0, 'tArr', 0.5*M.tf(1, 1));
R = polish_rendezvous(M, 1, 1, bad, struct('gates', false));
fprintf('  refusal: ok %d (%s)\n', R.ok, R.reason);
assert(~R.ok && startsWith(R.reason, 'departure before epoch 0'), ...
       'refusal: %s', R.reason);

% ---- 3b. departure before epoch 0 at the root (spacecraft 1) ----
k1 = 1;
c1 = struct('iD', S.iD(k1), 'iA', S.iA(k1), 'nRev', S.nRev(k1), ...
            'tDep', S.tDep(k1), 'tArr', S.tArr(k1));
R1 = polish_rendezvous(M, k1, S.assign(k1), c1, struct('gates', true));
fprintf('  sc 1: ok %d (%s)  tDep %.6f ND  g %.2e  nEval %d\n', R1.ok, ...
        R1.reason, R1.tDep, R1.gResid, R1.nEval);
assert(~R1.ok && startsWith(R1.reason, 'departure before epoch 0'), ...
       'sc 1: %s', R1.reason);
assert(R1.tDep < 0 && abs(R1.gResid) < 1e-9, 'sc 1: root found with tDep < 0');
assert(isnan(R1.gates.passPW) && isnan(R1.gates.passConj) && ...
       strcmp(R1.gates.conjReason, 'skipped'), 'sc 1: no gates run');

% ---- 4. deliberate gate failure ----
F = polish_rendezvous(M, k, j, cand, struct('gates', false, ...
                      'gateOverride', struct('flyKmMax', 1e-9)));
fprintf('  override: ok %d (%s)\n', F.ok, F.reason);
assert(~F.ok && startsWith(F.reason, 'gate failed'), 'override: %s', F.reason);

% ---- 5. gate error through the forceError seam ----
E = polish_rendezvous(M, k, j, cand, struct('gates', true, ...
                      'gateOverride', struct('forceError', 'hyp')));
fprintf('  forced error: ok %d (%s)\n', E.ok, E.reason);
assert(~E.ok && startsWith(E.reason, 'gate error: hyp'), 'forced: %s', E.reason);
assert(isequal(fieldnames(E.gates), fieldnames(G.gates)) && ...
       isequal(fieldnames(E), fieldnames(G)), 'forced: struct shape');

% ---- 6. deliberate non-fly gate failure (dimS requirement) ----
D = polish_rendezvous(M, k, j, cand, struct('gates', true, ...
                      'gateOverride', struct('dimS', 2)));
fprintf('  dimS override: ok %d (%s)\n', D.ok, D.reason);
assert(~D.ok && startsWith(D.reason, 'gate failed') && ...
       contains(D.reason, 'hyp') && ~contains(D.reason, 'fly'), 'dimS: %s', D.reason);

fprintf('  PASS\n');
ok = true;
end
