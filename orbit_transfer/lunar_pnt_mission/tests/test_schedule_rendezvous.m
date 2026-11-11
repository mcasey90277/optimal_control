function ok = test_schedule_rendezvous()
%% Purpose:
%
%   Verifies schedule_rendezvous: (1) a hand-checkable analytic case,
%   (2) agreement with a brute-force oracle written here with plain nested
%   loops (no shared helper), (3) tolerance respected and admissible counts
%   per pair, (4) refusal with schedule_rendezvous:noAssignment, (5) the
%   six-column ranked matrices and the tDep matrix, (8) opts.exclude:
%   empty reproduces the result bit for bit (random case and real
%   catalog), the brute-force oracle with the same exclusions, a malformed
%   list refused, and on the real catalog excluding spacecraft 1's chosen
%   candidate moves that pair's best arrival to its next-ranked one;
%   (9) HEADLINE PINS on the real catalog: last arrival 18.537 d and
%   assignment [2 4 6 1 5 3] unpolished; with the study's three refused
%   candidates excluded, 18.799 d and [2 5 6 1 3 4].
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
%  M. Casey  exclusion checks (8)                             (c) 11/10/2026
%  M. Casey  headline pins (9)                                (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));

% ---- 1. analytic case ----
o1 = struct('tol', 0.01, 'nMaxRevs', 3, 'propKgPerNd', 1);
sA = [0 0.25 0.5 0.75];
S = schedule_rendezvous([2 2 2 2; 3 3 3 3], [0 0.5], sA, 1, 4, [0 0.5], [0 0.5], o1);
assert(isequal(S.assign, [1 2]), 'analytic assign');
assert(isequal(S.best.tArr, [2 2; 3 3]), 'analytic best tArr');
assert(isequal(S.best.iD, [1 1; 2 2]) && isequal(S.best.nRev, zeros(2)), 'analytic iD, n');
assert(isequal(S.best.iA, [3 1; 4 2]), 'analytic iA');
assert(S.tLast == 3, 'analytic tLast');
assert(all(abs(S.resid) < 1e-12) && all(abs(S.best.resid(:)) < 1e-12), 'analytic resid');
assert(isequal(S.tDep, [0 0]) && isequal(S.best.tDep, zeros(2)), 'analytic tDep');
assert(all(size(S.ranked{1, 1}) == [size(S.ranked{1, 1}, 1) 6]), 'ranked has six columns');
assert(S.ranked{1, 1}(1, 1) == 2 && S.ranked{2, 2}(1, 1) == 3, 'ranked sorted');
fprintf('  analytic: tLast %g  assign [%d %d]  max|resid| %.2g\n', ...
        S.tLast, S.assign, max(abs(S.resid)));

% ---- 2. brute force oracle ----
rng(7);
nD = 6; nA = 12; K = 3;
tf = 2 + rand(nD, nA);
sD2 = (0:nD-1)/nD;  sA2 = (0:nA-1)/nA;
tauD = 1; tauT = 5.2; dp0 = (0:2)/3; sp0 = [0.11 0.47 0.82];
tol = 0.02; nMax = 6;
o2 = struct('tol', tol, 'nMaxRevs', nMax, 'propKgPerNd', 1.7);
S2 = schedule_rendezvous(tf, sD2, sA2, tauD, tauT, dp0, sp0, o2);
bfArr = inf(K, K);  bfCnt = zeros(K, K);  bfProp = inf(K, K);
for kk = 1:K
    for jj = 1:K
        for id = 1:nD
            for ia = 1:nA
                for nn = 0:nMax
                    d = mod(sD2(id) - dp0(kk), 1);
                    if d > 1 - 1e-9, d = 0; end
                    tdep = tauD*d + nn*tauD;
                    tarr = tdep + tf(id, ia);
                    r = mod(sA2(ia) - mod(sp0(jj) + tarr/tauT, 1), 1);
                    if r > 0.5, r = r - 1; end
                    if abs(r) <= tol
                        bfCnt(kk, jj) = bfCnt(kk, jj) + 1;
                        pk = 1.7*tf(id, ia);
                        if tarr < bfArr(kk, jj) - 1e-12 || ...
                           (abs(tarr - bfArr(kk, jj)) <= 1e-12 && pk < bfProp(kk, jj))
                            bfArr(kk, jj) = tarr;  bfProp(kk, jj) = pk;
                        end
                    end
                end
            end
        end
    end
end
allP = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];   % lexicographic
bfLast = inf;  bfAssign = [];  bfTot = inf;
for kp = 1:6
    p = allP(kp, :);
    tl = 0;  tp = 0;
    for kk = 1:K
        tl = max(tl, bfArr(kk, p(kk)));
        tp = tp + bfProp(kk, p(kk));
    end
    if tl < bfLast - 1e-12 || (abs(tl - bfLast) <= 1e-12 && tp < bfTot - 1e-12)
        bfLast = tl;  bfAssign = p;  bfTot = tp;
    end
end
assert(isfinite(bfLast), 'oracle found no assignment (bad test seed)');
assert(abs(S2.tLast - bfLast) <= 1e-12, 'brute tLast %g vs %g', S2.tLast, bfLast);
assert(isequal(S2.assign, bfAssign), 'brute assign');
bt = S2.best.tArr;  bt(isnan(bt)) = inf;
assert(max(abs(bt(isfinite(bfArr)) - bfArr(isfinite(bfArr)))) <= 1e-12 && ...
       isequal(isfinite(bt), isfinite(bfArr)), 'brute best.tArr');
assert(abs(S2.totalPropKg - bfTot) <= 1e-9, 'brute total propellant');
fprintf('  brute: tLast %.6f  assign [%d %d %d]  total prop %.4f kg\n', ...
        S2.tLast, S2.assign, S2.totalPropKg);

% ---- 3. tolerance respected, counts per pair ----
assert(all(abs(S2.resid) <= tol), 'chosen resid within tol');
assert(isequal(S2.nAdmissible, bfCnt), 'nAdmissible vs brute');
for kk = 1:K
    for jj = 1:K
        R = S2.ranked{kk, jj};
        assert(size(R, 1) == bfCnt(kk, jj) && all(abs(R(:, 5)) <= tol) && ...
               issorted(R(:, 1)), 'ranked pair (%d,%d)', kk, jj);
    end
end
fprintf('  counts per pair:%s\n', sprintf(' %d', S2.nAdmissible'));

% ---- 4. refusal ----
o4 = struct('tol', 1e-6, 'nMaxRevs', 3, 'propKgPerNd', 1);
msgId = '';
try
    schedule_rendezvous([2 2 2 2; 3 3 3 3], [0 0.5], [0.1 0.35 0.6 0.85], ...
                        1, 4, [0 0.5], [0 0.5], o4);
catch err
    msgId = err.identifier;
    fprintf('  refusal message: %s\n', strrep(err.message, newline, ' | '));
end
assert(strcmp(msgId, 'schedule_rendezvous:noAssignment'), 'refusal id "%s"', msgId);

% ---- 5. an empty pair with a complete permutation through the others ----
sA5 = [0.5 0.25];  tf5 = [2 2; 3 3];
o5 = struct('tol', 0.01, 'nMaxRevs', 0, 'propKgPerNd', 1);
S5 = schedule_rendezvous(tf5, [0 0.5], sA5, 1, 4, [0 0.5], [0 0.5], o5);
assert(isnan(S5.best.tArr(1, 2)) && S5.nAdmissible(1, 2) == 0, 'empty pair NaN');
assert(isequal(S5.assign, [1 2]) && S5.tLast == 3, 'assignment around the empty pair');
assert(isempty(S5.ranked{1, 2}), 'empty ranked');

% ---- 6. Hall violation with finite entries ----
sA6 = [0.5 0.75];
msgId = '';  msg6 = '';
try
    schedule_rendezvous([2 2; 3 3], [0 0.5], sA6, 1, 4, [0 0.5], [0 0.5], o5);
catch err
    msgId = err.identifier;  msg6 = err.message;
end
assert(strcmp(msgId, 'schedule_rendezvous:noAssignment'), 'Hall violation id "%s"', msgId);
% finite entries exist (slot 1 reachable by both): the message names only
% the empty pairs and carries finite numbers in the table
assert(contains(msg6, '(sc 1, slot 2)') && ~contains(msg6, '(sc 1, slot 1)') && ...
       ~contains(msg6, '(sc 2, slot 1)'), 'Hall message pairs');

% ---- 7. snap: a departure phase equal to the DRO phase up to noise ----
a7 = 3/5;  b7 = 3*(1/5);
assert(a7 - b7 < 0, 'snap test needs negative floating noise');
S7 = schedule_rendezvous(2, a7, 0.5, 1, 4, b7, 0, ...
                         struct('tol', 0.01, 'nMaxRevs', 0, 'propKgPerNd', 1));
assert(S7.tDep == 0 && S7.tArr == 2, 'snap: tDep %g', S7.tDep);
fprintf('  snap: sD - droPhase0 = %.3g, tDep = %g\n', a7 - b7, S7.tDep);

% ---- 8. exclusion list (opts.exclude rows [k j iD iA nRev]) ----
% 8a. an empty exclude reproduces the result without the field, bit for bit
o8 = o2;  o8.exclude = zeros(0, 5);
S8 = schedule_rendezvous(tf, sD2, sA2, tauD, tauT, dp0, sp0, o8);
assert(isequal(rmfield(S8, 'opts'), rmfield(S2, 'opts')), 'empty exclude changes the result');
% 8b. brute-force oracle with exclusions: every chosen candidate of S2 plus
% the second-ranked candidate of pair (1,1) are excluded
ex = [(1:K).', S2.assign.', S2.iD.', S2.iA.', S2.nRev.'];
if size(S2.ranked{1, 1}, 1) >= 2
    ex = [ex; 1 1 S2.ranked{1, 1}(2, 2:4)];
end
o8b = o2;  o8b.exclude = ex;
S8b = schedule_rendezvous(tf, sD2, sA2, tauD, tauT, dp0, sp0, o8b);
bfArr = inf(K, K);  bfCnt = zeros(K, K);  bfProp = inf(K, K);
for kk = 1:K
    for jj = 1:K
        for id = 1:nD
            for ia = 1:nA
                for nn = 0:nMax
                    if any(ex(:, 1) == kk & ex(:, 2) == jj & ex(:, 3) == id & ...
                           ex(:, 4) == ia & ex(:, 5) == nn)
                        continue
                    end
                    d = mod(sD2(id) - dp0(kk), 1);
                    if d > 1 - 1e-9, d = 0; end
                    tdep = tauD*d + nn*tauD;
                    tarr = tdep + tf(id, ia);
                    r = mod(sA2(ia) - mod(sp0(jj) + tarr/tauT, 1), 1);
                    if r > 0.5, r = r - 1; end
                    if abs(r) <= tol
                        bfCnt(kk, jj) = bfCnt(kk, jj) + 1;
                        pk = 1.7*tf(id, ia);
                        if tarr < bfArr(kk, jj) - 1e-12 || ...
                           (abs(tarr - bfArr(kk, jj)) <= 1e-12 && pk < bfProp(kk, jj))
                            bfArr(kk, jj) = tarr;  bfProp(kk, jj) = pk;
                        end
                    end
                end
            end
        end
    end
end
bfLast = inf;  bfAssign = [];  bfTot = inf;
for kp = 1:6
    p = allP(kp, :);
    tl = 0;  tp = 0;
    for kk = 1:K
        tl = max(tl, bfArr(kk, p(kk)));
        tp = tp + bfProp(kk, p(kk));
    end
    if tl < bfLast - 1e-12 || (abs(tl - bfLast) <= 1e-12 && tp < bfTot - 1e-12)
        bfLast = tl;  bfAssign = p;  bfTot = tp;
    end
end
assert(isfinite(bfLast), 'oracle with exclusions found no assignment (bad test seed)');
assert(abs(S8b.tLast - bfLast) <= 1e-12, 'excluded brute tLast %g vs %g', S8b.tLast, bfLast);
assert(isequal(S8b.assign, bfAssign), 'excluded brute assign');
assert(isequal(S8b.nAdmissible, bfCnt), 'excluded nAdmissible vs brute');
bt = S8b.best.tArr;  bt(isnan(bt)) = inf;
assert(isequal(isfinite(bt), isfinite(bfArr)) && ...
       max(abs(bt(isfinite(bfArr)) - bfArr(isfinite(bfArr)))) <= 1e-12, 'excluded best.tArr');
for kk = 1:K
    for jj = 1:K
        R = S8b.ranked{kk, jj};
        hit = ex(ex(:, 1) == kk & ex(:, 2) == jj, 3:5);
        assert(~any(ismember(R(:, 2:4), hit, 'rows')), 'excluded row survives in pair (%d,%d)', kk, jj);
    end
end
fprintf('  exclude (oracle): %d rows, tLast %.6f -> %.6f, assign [%d %d %d]\n', ...
        size(ex, 1), S2.tLast, S8b.tLast, S8b.assign);
% 8c. malformed exclude is refused
msgId = '';
o8c = o2;  o8c.exclude = [1 2 3];
try
    schedule_rendezvous(tf, sD2, sA2, tauD, tauT, dp0, sp0, o8c);
catch err
    msgId = err.identifier;
end
assert(strcmp(msgId, 'schedule_rendezvous:exclude'), 'malformed exclude id "%s"', msgId);
% 8d. REAL catalog: excluding spacecraft 1's chosen candidate for its slot
% moves that pair's best arrival to its next-ranked candidate; an empty
% exclude reproduces the committed result bit for bit
M = mission_inputs();
oR = struct('tol', 0.005, 'nMaxRevs', 10, 'propKgPerNd', M.propKgPerNd);
SR = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, oR);
oRe = oR;  oRe.exclude = zeros(0, 5);
SRe = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, oRe);
assert(isequal(rmfield(SRe, 'opts'), rmfield(SR, 'opts')), 'real catalog: empty exclude changes the result');
j1 = SR.assign(1);
R1 = SR.ranked{1, j1};
oRx = oR;  oRx.exclude = [1 j1 SR.iD(1) SR.iA(1) SR.nRev(1)];
SRx = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, oRx);
assert(size(R1, 1) >= 2, 'real catalog: pair (1,%d) has a single candidate', j1);
assert(isequal(SRx.ranked{1, j1}, R1(2:end, :)), 'real catalog: ranked list not shifted by one');
assert(SRx.best.tArr(1, j1) == R1(2, 1) && SRx.best.iD(1, j1) == R1(2, 2) && ...
       SRx.best.iA(1, j1) == R1(2, 3) && SRx.best.nRev(1, j1) == R1(2, 4), ...
       'real catalog: best of pair (1,%d) is not the next-ranked candidate', j1);
assert(SRx.nAdmissible(1, j1) == SR.nAdmissible(1, j1) - 1, 'real catalog: count');
d = M.tStar/86400;
fprintf('  exclude (real): pair (1,%d) best %.4f d -> %.4f d; last arrival %.4f d -> %.4f d\n', ...
        j1, SR.best.tArr(1, j1)*d, SRx.best.tArr(1, j1)*d, SR.tLast*d, SRx.tLast*d);

% ---- 9. HEADLINE PINS on the real catalog (README, Measured results) ----
% 9a. the unpolished optimum
assert(abs(SR.tLast*d - 18.537) < 1e-3, 'real catalog: last arrival %.6f d (pin 18.537)', SR.tLast*d);
assert(isequal(SR.assign, [2 4 6 1 5 3]), 'real catalog: assignment [%s] (pin [2 4 6 1 5 3])', ...
       sprintf('%d ', SR.assign));
% 9b. with the study's three refused candidates excluded (rows [k j iD iA
% nRev] read from the saved full run's rounds: sc 1 slot 2 cell (1,1)
% and sc 2 slot 4 cell (5,15) refused in round 1, sc 6 slot 4 cell
% (21,16) in round 2; all zero-wait)
exStudy = [1 2 1 1 0; 2 4 5 15 0; 6 4 21 16 0];
oRs = oR;  oRs.exclude = exStudy;
SRs = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, oRs);
for kx = 1:2                                   % the round-1 refusals are round-1 choices
    kc = exStudy(kx, 1);
    assert(isequal(exStudy(kx, :), [kc SR.assign(kc) SR.iD(kc) SR.iA(kc) SR.nRev(kc)]), ...
           'exclusion row %d is not spacecraft %d''s round-1 choice', kx, kc);
end
assert(abs(SRs.tLast*d - 18.799) < 1e-3, 'excluded: last arrival %.6f d (pin 18.799)', SRs.tLast*d);
assert(isequal(SRs.assign, [2 5 6 1 3 4]), 'excluded: assignment [%s] (pin [2 5 6 1 3 4])', ...
       sprintf('%d ', SRs.assign));
fprintf('  pins: unpolished %.4f d [%s]; three refusals excluded %.4f d [%s]\n', SR.tLast*d, ...
        sprintf('%d ', SR.assign), SRs.tLast*d, sprintf('%d ', SRs.assign));
ok = true;
fprintf('  PASS\n');
end
