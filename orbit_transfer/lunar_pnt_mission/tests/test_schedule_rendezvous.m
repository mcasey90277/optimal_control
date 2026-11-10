function ok = test_schedule_rendezvous()
%% Purpose:
%
%   Verifies schedule_rendezvous: (1) a hand-checkable analytic case,
%   (2) agreement with a brute-force oracle written here with plain nested
%   loops (no shared helper), (3) tolerance respected and admissible counts
%   per pair, (4) refusal with schedule_rendezvous:noAssignment, (5) the
%   six-column ranked matrices and the tDep matrix.
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
ok = true;
fprintf('  PASS\n');
end
