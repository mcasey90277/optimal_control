function ok = test_schedule_gto_dro()
%% Purpose:
%
%   Verifies schedule_gto_dro on synthetic tables small enough to check by
%   hand: (1) an analytic 2-departure-phase, 6-arrival-phase case (T0 and
%   every tDep computed here from the formulas), (2) the orientation
%   bound (the chosen orient0 lies in the only admissible window), (3)
%   relabelling (leg-2 satellite 1 is at DRO phase 0 at T0, recomputed
%   from tArr and the residual), (4) refusal with
%   schedule_gto_dro:noAssignment and the six candidate counts in the
%   message, (5) opts.exclude (the previous choice is removed, the next
%   candidate is returned), (6) the propellant tie-break between triples
%   with equal T0, (7) agreement with a brute-force triple-by-triple
%   oracle written here with plain loops, (8) a malformed exclude list is
%   refused.
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
%  M. Casey                                                   (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));

sA6   = (0:5)/6;
tauG  = 1/6;  tauD = 1;  sgn = -1;
base  = struct('orient0Grid', 0, 'gtoPhase0Grid', 0, 'deltaGrid', 0.05, ...
               'nMaxRevs', 3, 'tolPhase', 0.06, 'orientBound', 25, ...
               'propKgPerNd', 1);

% ---- 1. analytic case ----
% tauG = 1/6 puts the departure instants on the slot grid; cell (iD=1)
% has tf = 2, cell (iD=2) tf = 2.04. For n = 0 slot j takes iA = j at
% tDep = 0, tArr = 2, resid = +-0.05 (delta = 0.05, tol 0.06).
tf1 = cat(1, 2*ones(1, 6), 2.04*ones(1, 6));            % 2 x 6 x 1
S = schedule_gto_dro(tf1, [0 0.5], sA6, 0, tauG, tauD, sgn, base);
T0exp = 2 + 1/6 - 0.05;                                 % mod(0.05 + t, 1/6) = 0
assert(abs(S.lastArr - 2) < 1e-12, 'analytic lastArr %.12g', S.lastArr);
assert(abs(S.T0 - T0exp) < 1e-12, 'analytic T0 %.12g vs %.12g', S.T0, T0exp);
assert(all(S.iD == 1) && all(S.n == 0) && all(S.iS == 1), 'analytic iD, n, iS');
assert(all(abs(S.tDep) < 1e-12) && all(abs(S.tArr - 2) < 1e-12), 'analytic tDep/tArr');
assert(all(abs(S.tf - 2) < 1e-12) && all(abs(S.orientDeg) < 1e-12), 'analytic tf, orient');
assert(all(abs(abs(S.resid) - 0.05) < 1e-12), 'analytic resid');
assert(isequal(sort(S.iA), 1:6), 'analytic iA is a permutation');
assert(all(S.nCand >= 1) && numel(S.cands) == 6, 'analytic cands');
assert(S.orient0 == 0 && S.gtoPhase0 == 0 && S.delta == 0.05, 'analytic triple');
fprintf('  analytic: T0 %.6f (expect %.6f)  lastArr %g  max|tDep| %.2g\n', ...
        S.T0, T0exp, S.lastArr, max(abs(S.tDep)));

% T0 guard: a last arrival that is on the slot grid up to floating noise
% above must not cost a whole slot spacing
o1g = base;  o1g.deltaGrid = 0;
Sg = schedule_gto_dro(tf1 + 2^-51, [0 0.5], sA6, 0, tauG, tauD, sgn, o1g);
assert(6*Sg.lastArr > 12 && abs(Sg.T0 - Sg.lastArr) < 1e-12, ...
       'T0 guard: lastArr %.17g  T0 %.17g', Sg.lastArr, Sg.T0);
fprintf('  guard: lastArr %.17g  T0 %.17g\n', Sg.lastArr, Sg.T0);

% ... and a last arrival clearly above a grid instant (1e-9) jumps to the next
Sg2 = schedule_gto_dro(tf1 + 1e-9, [0 0.5], sA6, 0, tauG, tauD, sgn, o1g);
assert(abs(Sg2.T0 - (Sg2.lastArr - 1e-9 + 1/6)) < 1e-12, ...
       'T0 above guard: lastArr %.12g  T0 %.12g', Sg2.lastArr, Sg2.T0);

% ---- 1b. centred triple wins a T0 + propellant tie ----
% one departure at tDep = 0, so orientDeg = orient0; orient0 in {-20, 0, 15}
% all admissible (bound 25), identical T0 and propellant: the smallest
% departure offset from the sheet orientation (0 deg) must win, not the
% lexicographically smallest orient0 = -20
oC = base;  oC.orient0Grid = [-20 0 15];  oC.deltaGrid = 0;  oC.tolPhase = 0.02;
SC = schedule_gto_dro(2*ones(1, 6), 0, sA6, 0, tauG, tauD, sgn, oC);
assert(SC.orient0 == 0 && abs(SC.maxOrientOffsetDeg) < 1e-9, ...
       'centred: orient0 %g offset %g', SC.orient0, SC.maxOrientOffsetDeg);
oC.orient0Grid = [-20 15];
SC2 = schedule_gto_dro(2*ones(1, 6), 0, sA6, 0, tauG, tauD, sgn, oC);
assert(SC2.orient0 == 15 && abs(SC2.maxOrientOffsetDeg - 15) < 1e-9, 'centred 2: %g', SC2.orient0);
fprintf('  centred: orient0 %g (offset %g deg); of {-20, 15}: %g (offset %g)\n', ...
        SC.orient0, SC.maxOrientOffsetDeg, SC2.orient0, SC2.maxOrientOffsetDeg);

% ---- 2. orientation bound ----
% sheet 1 (orient 0) unusable, sheet 2 (orient 90) usable: tDep = 0 so
% orientDeg = orient0, admissible iff |wrap(orient0 - 90)| <= 25
tf2 = cat(3, nan(2, 6), tf1(:, :, 1));
o2  = base;  o2.orient0Grid = -180:5:175;
S2 = schedule_gto_dro(tf2, [0 0.5], sA6, [0 90], tauG, tauD, sgn, o2);
assert(S2.orient0 >= 65 && S2.orient0 <= 115, 'orient0 %g outside [65, 115]', S2.orient0);
assert(all(S2.iS == 2), 'orientation: sheet 2 only');
assert(S2.orient0 == 90 && abs(S2.maxOrientOffsetDeg) < 1e-9, 'centred orient0 (%g)', S2.orient0);
fprintf('  orientation: chosen orient0 = %g (window [65, 115])\n', S2.orient0);

% ---- 3. relabelling ----
for kk = 1:6
    phT0 = mod(sA6(S.iA(kk)) - S.resid(kk) + (S.T0 - S.tArr(kk))/tauD, 1);
    d = abs(phT0 - (kk-1)/6);
    d = min(d, 1 - d);
    assert(d < 1e-9, 'satellite %d phase at T0 %.9f vs %.9f', kk, phT0, (kk-1)/6);
end
assert(S.slot(1) == 6, 'slot at phase 0 at T0 is slot 6 here (slot %d)', S.slot(1));
assert(isequal(sort(S.slot), 1:6), 'relabel is a permutation');
fprintf('  relabel: satellite 1 on slot %d, phase at T0 = %.3g\n', S.slot(1), ...
        mod(sA6(S.iA(1)) - S.resid(1) + (S.T0 - S.tArr(1))/tauD, 1));

% ---- 4. refusal ----
tf4 = tf1;  tf4(:, 3) = NaN;  tf4(2, :) = NaN;          % slot 3 needs iA = 3
o4 = base;  o4.nMaxRevs = 0;  o4.tolPhase = 0.06;
msgId = '';  msg4 = '';
try
    schedule_gto_dro(tf4, [0 0.5], sA6, 0, tauG, tauD, sgn, o4);
catch err
    msgId = err.identifier;  msg4 = err.message;
end
assert(strcmp(msgId, 'schedule_gto_dro:noAssignment'), 'refusal id "%s"', msgId);
assert(contains(msg4, '1 1 0 1 1 1'), 'message lacks the six counts: %s', msg4);
fprintf('  refusal: %s\n', msg4);

% ---- 5. exclude (hand-derived values) ----
% slot 1 at n = 0, iD = 1: tDep = 0, tArr = 2, iA = 1. The next candidate
% of slot 1 is iD = 2 (tDep = tauG/2 = 1/12, tf = 2.04), n = 0:
% tArr = 1/12 + 2.04 = 2.123333, DRO phase mod(0.05 + tArr, 1) = 0.173333,
% nearest arrival phase sA(2) = 1/6.
tArrNext = tauG*0.5 + 2.04;
iANext   = round(6*mod(0.05 + tArrNext, 1)) + 1;
assert(iANext == 2, 'hand derivation of the next iA');
o5 = base;  o5.exclude = [1 1 1 0];
S5 = schedule_gto_dro(tf1, [0 0.5], sA6, 0, tauG, tauD, sgn, o5);
k5 = find(S5.slot == 1);
assert(S5.iD(k5) == 2 && S5.iA(k5) == iANext && S5.n(k5) == 0 && ...
       abs(S5.tArr(k5) - tArrNext) < 1e-12, 'next candidate not returned');
assert(abs(S5.tDep(k5) - tauG/2) < 1e-12, 'next candidate tDep');
c5 = S5.cands{k5};
assert(c5(1).iD == 1 && c5(1).iA == 1 && c5(1).n == 0 && c5(1).excluded && ...
       abs(c5(1).tArr - 2) < 1e-12 && ~c5(2).excluded, 'cands list / excluded flag');
fprintf('  exclude: slot 1 tArr 2 -> %.4f (iD %d iA %d n %d)\n', ...
        S5.tArr(k5), S5.iD(k5), S5.iA(k5), S5.n(k5));

% ---- 5b. lockTriple: re-scheduling one satellite moves nobody else ----
% All six tie at tArr = 2 on the base table. Bar slot 1's choice: its
% replacement arrives at 2.123333, now the bottleneck, so
% T0 = smallest t >= 2.123333 with mod(0.05 + t, 1/6) = 0
%    = 2 + 2/6 - 0.05, one slot spacing later than before; every phase at
% T0 advances by 1/6, so slot j becomes satellite k = mod(j + 1, 6) + 1
% (it was mod(j, 6) + 1).
oL = base;  oL.exclude = [1 1 1 0];
oL.orient0Grid = 90;  oL.gtoPhase0Grid = 0.5;  oL.deltaGrid = 0.3;   % must be ignored
oL.lockTriple = struct('orient0', 0, 'gtoPhase0', 0, 'delta', 0.05);
SL = schedule_gto_dro(tf1, [0 0.5], sA6, 0, tauG, tauD, sgn, oL);
assert(SL.orient0 == 0 && SL.gtoPhase0 == 0 && SL.delta == 0.05, 'locked triple changed');
T0new = 2 + 2/6 - 0.05;
assert(abs(SL.T0 - T0new) < 1e-12, 'locked T0 %.12g vs %.12g', SL.T0, T0new);
assert(abs(SL.lastArr - tArrNext) < 1e-12, 'locked lastArr');
for js = 1:6
    kOld = find(S.slot == js);  kNew = find(SL.slot == js);
    assert(kNew == mod(js + 1, 6) + 1 && kOld == mod(js, 6) + 1, 'labels slot %d', js);
    if js == 1, continue, end
    same = S.iD(kOld) == SL.iD(kNew) && S.iA(kOld) == SL.iA(kNew) && ...
           S.iS(kOld) == SL.iS(kNew) && S.n(kOld) == SL.n(kNew) && ...
           S.tDep(kOld) == SL.tDep(kNew) && S.tArr(kOld) == SL.tArr(kNew);
    assert(same, 'slot %d transfer moved under lockTriple', js);
end
assert(~isequal(S.slot, SL.slot), 'labels should change');
% locked triple infeasible under exclude: throws noAssignment
oLx = oL;
oLx.exclude = [ones(6, 1), (1:6).', ones(6, 1), zeros(6, 1)];   % bars every iA of iD 1 at n = 0
oLx.nMaxRevs = 0;  oLx.tolPhase = 0.01;  oLx.lockTriple.delta = 0;                          % ... and no other candidate
msgId = '';
try
    schedule_gto_dro(tf1(1, :), 0, sA6, 0, tauG, tauD, sgn, oLx);
catch err
    msgId = err.identifier;
end
assert(strcmp(msgId, 'schedule_gto_dro:noAssignment'), 'locked infeasible id "%s"', msgId);
% (b) without the lock the same exclude is simply re-optimised over a grid
oF = base;  oF.exclude = [1 1 1 0];  oF.deltaGrid = [0 0.05 0.1];
SF = schedule_gto_dro(tf1, [0 0.5], sA6, 0, tauG, tauD, sgn, oF);
assert(isfinite(SF.T0) && SF.T0 >= SF.lastArr - 1e-12 && ...
       all(arrayfun(@(kk) ~SF.cands{kk}(find(abs([SF.cands{kk}.tArr] - SF.tArr(kk)) < 1e-12, 1)).excluded, 1:6)), ...
       'unlocked re-schedule infeasible');
fprintf('  lockTriple: T0 %.6f -> %.6f (expect %.6f), slot 1 label %d -> %d; unlocked: delta %g T0 %.6f\n', ...
        S.T0, SL.T0, T0new, find(S.slot == 1), find(SL.slot == 1), SF.delta, SF.T0);

% ---- 6. tie to less propellant ----
% one departure phase, two sheets (orient 0, 90). Triple (orient0 0,
% gtoPhase0 0): sheet 1, tf 2, tDep 0. Triple (90, 0.5): sheet 2,
% tDep = tauG/2 = 1/12, tf = 2 - 1/12; same tArr = 2 and same T0, less
% propellant. (0, 0.5) and (90, 0) miss the arrival phase by 1/12.
tfA = 2*ones(1, 6);  tfB = (2 - 1/12)*ones(1, 6);
tf6 = cat(3, tfA, tfB);                                  % 1 x 6 x 2
o6 = base;  o6.orient0Grid = [0 90];  o6.gtoPhase0Grid = [0 0.5];
o6.deltaGrid = 0;  o6.tolPhase = 0.02;  o6.propKgPerNd = 3;
S6 = schedule_gto_dro(tf6, 0, sA6, [0 90], tauG, tauD, sgn, o6);
assert(S6.orient0 == 90 && S6.gtoPhase0 == 0.5, 'tie-break chose (%g, %g)', ...
       S6.orient0, S6.gtoPhase0);
assert(abs(S6.totalPropKg - 6*3*(2 - 1/12)) < 1e-12, 'propellant %.9g', S6.totalPropKg);
assert(abs(S6.T0 - 2) < 1e-9, 'tie T0 %.12g', S6.T0);
o6b = o6;  o6b.orient0Grid = 0;  o6b.gtoPhase0Grid = 0;
S6b = schedule_gto_dro(tf6, 0, sA6, [0 90], tauG, tauD, sgn, o6b);
assert(abs(S6b.T0 - S6.T0) < 1e-9 && S6b.totalPropKg > S6.totalPropKg + 1, 'tie set-up');
fprintf('  tie: chose (%g, %g), %.4f kg vs %.4f kg for the lexicographic one\n', ...
        S6.orient0, S6.gtoPhase0, S6.totalPropKg, S6b.totalPropKg);

% ---- 7. brute-force oracle on a random table ----
rng(11);
nD = 4;  nA = 6;  nS = 2;  nMax = 4;
tf7 = 2 + 3*rand(nD, nA, nS);  tf7(rand(nD, nA, nS) < 0.15) = NaN;
sD7 = (0:nD-1)/nD;  sA7 = (0:nA-1)/nA;  sO7 = [-20 30];
tG = 0.37;  tD = 1.13;  sg = -1;
g7 = struct('orient0Grid', [-30 0 25 60], 'gtoPhase0Grid', [0 0.25 0.5], ...
            'deltaGrid', [0 0.04 0.1], 'nMaxRevs', nMax, 'tolPhase', 0.08, ...
            'orientBound', 25, 'propKgPerNd', 1.9);
S7 = schedule_gto_dro(tf7, sD7, sA7, sO7, tG, tD, sg, g7);
bestT0 = inf;  bestProp = inf;  bestTrip = [];
for o0 = g7.orient0Grid
  for g0 = g7.gtoPhase0Grid
    for dl = g7.deltaGrid
      tarrSlot = nan(1, 6);  propSlot = nan(1, 6);  full = true;
      for js = 1:6
        bt = inf;  bp = inf;
        for id = 1:nD
          for ia = 1:nA
            for is = 1:nS
              for nr = 0:nMax
                dp = mod(sD7(id) - g0, 1);
                if dp > 1 - 1e-9, dp = 0; end
                td = tG*(dp + nr);
                od = o0 + sg*(180/pi)*td;
                wr = mod(od - sO7(is) + 180, 360) - 180;
                if abs(wr) > 25 || isnan(tf7(id, ia, is)), continue, end
                ta = td + tf7(id, ia, is);
                r = mod(sA7(ia) - mod(dl + (js-1)/6 + ta/tD, 1), 1);
                if r > 0.5, r = r - 1; end
                if abs(r) > 0.08, continue, end
                if ta < bt - 1e-12 || (abs(ta - bt) <= 1e-12 && tf7(id, ia, is) < bp)
                    bt = ta;  bp = tf7(id, ia, is);
                end
              end
            end
          end
        end
        if isinf(bt), full = false;  break, end
        tarrSlot(js) = bt;  propSlot(js) = bp;
      end
      if ~full, continue, end
      la = max(tarrSlot);
      kq = ceil(6*(dl + la/tD) - 1e-12);
      t0 = tD*(kq/6 - dl);
      pr = 1.9*sum(propSlot);
      if t0 < bestT0 - 1e-12 || (abs(t0 - bestT0) <= 1e-12 && pr < bestProp - 1e-12)
          bestT0 = t0;  bestProp = pr;  bestTrip = [o0 g0 dl];
      end
    end
  end
end
assert(~isempty(bestTrip), 'oracle found no feasible triple (bad seed)');
assert(abs(S7.T0 - bestT0) < 1e-10, 'oracle T0 %.10g vs %.10g', S7.T0, bestT0);
assert(isequal([S7.orient0 S7.gtoPhase0 S7.delta], bestTrip), 'oracle triple [%s]', ...
       sprintf('%g ', [S7.orient0 S7.gtoPhase0 S7.delta]));
assert(abs(S7.totalPropKg - bestProp) < 1e-9, 'oracle propellant');
fprintf('  oracle: T0 %.6f  triple [%g %g %g]  prop %.4f kg\n', S7.T0, bestTrip, bestProp);

% ---- 8. malformed exclude ----
bad = {[1 2 3], [1.5 1 1 0], [0 1 1 0], [3 1 1 0], [1 7 1 0], [1 1 2 0], [1 1 1 -1], ...
       [1 1 1 4], [1 1 1 NaN], 'abc'};
for kb = 1:numel(bad)
    msgId = '';
    o8 = base;  o8.exclude = bad{kb};
    try
        schedule_gto_dro(tf1, [0 0.5], sA6, 0, tauG, tauD, sgn, o8);
    catch err
        msgId = err.identifier;
    end
    assert(strcmp(msgId, 'schedule_gto_dro:exclude'), 'malformed exclude %d id "%s"', kb, msgId);
end
fprintf('  malformed exclude: %d cases refused\n', numel(bad));

ok = true;
fprintf('  PASS\n');
end
