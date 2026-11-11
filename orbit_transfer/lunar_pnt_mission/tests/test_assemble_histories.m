function ok = test_assemble_histories()
%% Purpose:
%
%   Verifies assemble_histories: (1) on-slots mode over one tulip period
%   (shape, phase, slot positions at t = 0, Moon distance band); (2) one
%   real polished transfer (spacecraft 3, gates off; the other five marked
%   onStation) -- phase sequence 1 -> 2 -> 3 with the switches at tDep and
%   tArr to one sample, hand-over jumps, no teleport; (3) deliberate
%   failure: tArr shifted by +0.05 ND breaks the arrival continuity number;
%   (4) refusals: an unaccepted P that is not onStation, and the periodic
%   interpolants against the phase_state handles; wall time of a 50-day,
%   300 s, six-spacecraft history is printed.
%
%   Spacecraft 3, not 1: spacecraft 1 is refused by the polish (its exact
%   slot needs a departure before clock 0), measured in Task 3.
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
kmPerNd = M.lStar;
moon = [1 - M.muStar, 0, 0];

% ---- 1. on-slots over one tulip period ----
H = assemble_histories(M, [], 1:6, struct('onSlots', true, 'tEnd', M.tauT));
nn = numel(H.t);
assert(isequal(size(H.r), [nn 3 6]), 'size of H.r');
assert(all(H.phase(:) == 3), 'phase must be 3 everywhere');
assert(H.t(1) == 0 && abs(H.t(2) - 300/M.tStar) < 1e-15, 'time axis');
assert(all(isnan(H.tDep)) && all(isnan(H.tArr)), 'tDep/tArr in onSlots mode');
assert(H.tFirst == 0 && H.tLast == 0, 'onSlots window');
for kk = 1:6
    want = M.arrState(M.slotPhase0(kk));
    assert(max(abs(squeeze(H.r(1, :, kk)).' - want(1:3))) < 1e-12, ...
           'slot %d position at t = 0', kk);
end
rEnd = squeeze(H.r(end, :, :)).';                       % 6 x 3
dMoon = sqrt(sum((rEnd - moon).^2, 2))*kmPerNd;
fprintf('  Moon distance at tEnd: %.0f .. %.0f km\n', min(dMoon), max(dMoon));
assert(all(dMoon > 1900 & dMoon < 120000), 'Moon distance band');
% permuted assignment: spacecraft m sits on slot assign(m), not slot m
perm = [3 1 2 6 4 5];
Hp = assemble_histories(M, [], perm, struct('onSlots', true, 'tEnd', M.tauT));
for kk = 1:6
    for ii = [1 400]
        tt = Hp.t(ii);
        want  = M.arrState(mod(M.slotPhase0(perm(kk)) + tt/M.tauT, 1));
        other = M.arrState(mod(M.slotPhase0(kk) + tt/M.tauT, 1));
        got = squeeze(Hp.r(ii, :, kk)).';
        assert(max(abs(got - want(1:3))) < 1e-12, 'permuted slot, sc %d', kk);
        assert(max(abs(got - other(1:3))) > 1e-3, 'sc %d sits on slot kk', kk);
    end
end
assert(isequal(Hp.assign, perm), 'H.assign');
% same set of epochs reached through a later slot phase: drive of the
% periodic interpolant against the handle at random phases
sTest = rand(1, 200);
rvH   = M.arrState(sTest);
rvS   = zeros(6, 200);
for kk = 1:200, rvS(:, kk) = M.arrState(sTest(kk)); end
assert(isequal(size(rvH), [6 200]) && max(abs(rvH(:) - rvS(:))) < 1e-14, ...
       'phase_state handle must accept a vector of phases');

% ---- 2. one real transfer ----
S = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, ...
                        M.slotPhase0, struct('propKgPerNd', M.propKgPerNd));
k = 3;
j = S.assign(k);
cand = struct('iD', S.iD(k), 'iA', S.iA(k), 'nRev', S.nRev(k), ...
              'tDep', S.tDep(k), 'tArr', S.tArr(k));
P3 = polish_rendezvous(M, k, j, cand, struct('gates', false));
assert(P3.ok, 'polish not ok: %s', P3.reason);
Pc = cell(1, 6);  Pc{k} = P3;
onSt = true(1, 6);  onSt(k) = false;
asg = S.assign(:).';
opts = struct('onStation', onSt, 'serviceDays', 10);
t0 = tic;
H = assemble_histories(M, Pc, asg, opts);
fprintf('  history: N = %d, wall %.2f s\n', numel(H.t), toc(t0));
ph = H.phase(:, k);
assert(all(diff(double(ph)) >= 0) && ph(1) == 1 && ph(end) == 3, 'phase order');
iDep = find(ph >= 2, 1);  iArr = find(ph == 3, 1);
dtNd = 300/M.tStar;
assert(abs(H.t(iDep) - P3.tDep) <= dtNd && H.t(iDep) >= P3.tDep - 1e-12, ...
       'departure switch');
assert(abs(H.t(iArr) - P3.tArr) <= dtNd && H.t(iArr) >= P3.tArr - 1e-12, ...
       'arrival switch');
assert(all(H.phase(:, [1 2 4 5 6]) == 3, 'all'), 'onStation phases');
fprintf('  sc %d jump: dep %.3e km  arr %.3e km\n', k, H.jumpKm.dep(k), ...
        H.jumpKm.arr(k));
assert(H.jumpKm.dep(k) < 1, 'departure jump');
assert(H.jumpKm.arr(k) < 100, 'arrival jump');
step = sqrt(sum(diff(squeeze(H.r(:, :, k))).^2, 2))*kmPerNd;
fprintf('  largest step sc %d: %.1f km per %g s (limit %.0f)\n', k, max(step), ...
        300, 3*300);
assert(max(step) < 3*300, 'teleport at a hand-over');
assert(isnan(H.jumpKm.dep(1)) && isnan(H.jumpKm.arr(1)), 'onStation jumps are NaN');
assert(abs(H.tEnd - (max(H.tArr) + 10*86400/M.tStar)) < dtNd, 'default tEnd');

assert(H.tFirst == P3.tArr && H.tLast == P3.tArr, 'assembly window');
Hd = assemble_histories(M, Pc, asg, struct('onStation', onSt));
assert(Hd.t(end) >= Hd.tLast + 30*86400/M.tStar - dtNd, 'default 30 d service window');
ia = iArr + 100;                                   % a sample after arrival
wantA = M.arrState(mod(M.slotPhase0(S.assign(k)) + H.t(ia)/M.tauT, 1));
assert(max(abs(squeeze(H.r(ia, :, k)).' - wantA(1:3))) < 1e-12, ...
       'post-arrival slot indirection');
fprintf('  window: tFirst = tLast = %.6f ND; default t(end) - tLast = %.3f d\n', ...
        Hd.tLast, (Hd.t(end) - Hd.tLast)*M.tStar/86400);

% ---- 3. deliberate failure: shifted tArr ----
Pb = P3;  Pb.tArr = P3.tArr + 0.05;
Pcb = cell(1, 6);  Pcb{k} = Pb;
try
    assemble_histories(M, Pcb, asg, opts);
    error('test:noRefusal', 'a shifted tArr was accepted by default');
catch err
    assert(strcmp(err.identifier, 'assemble_histories:epochs'), ...
           'shifted tArr id: %s', err.identifier);
end
optsInf = opts;  optsInf.epochTol = inf;
Hb = assemble_histories(M, Pcb, asg, optsInf);
fprintf('  shifted tArr: arr jump %.1f km\n', Hb.jumpKm.arr(k));
assert(Hb.jumpKm.arr(k) > 100, 'the continuity number must be able to fail');

% ---- 4. refusals and the struct-array form ----
Pbad = P3;  Pbad.ok = false;
Pcr = cell(1, 6);  Pcr{k} = Pbad;
try
    assemble_histories(M, Pcr, asg, opts);
    error('test:noRefusal', 'an unaccepted P was flown');
catch err
    assert(strcmp(err.identifier, 'assemble_histories:notAccepted'), ...
           'refusal id: %s', err.identifier);
end
% tEnd guards
for bad = {-1, [1 2], H.tLast - 1, NaN}
    optsT = opts;  optsT.tEnd = bad{1};
    try
        assemble_histories(M, Pc, asg, optsT);
        error('test:noRefusal', 'bad tEnd accepted');
    catch err
        assert(strcmp(err.identifier, 'assemble_histories:tEnd'), ...
               'tEnd id: %s', err.identifier);
    end
end
% duplicated flight sample: same history; reversed time: flightGrid
Pd = P3;
Pd.flight.t = [P3.flight.t(1:5); P3.flight.t(5); P3.flight.t(6:end)];
Pd.flight.Y = [P3.flight.Y(1:5, :); P3.flight.Y(5, :); P3.flight.Y(6:end, :)];
Pcd = cell(1, 6);  Pcd{k} = Pd;
Hdup = assemble_histories(M, Pcd, asg, opts);
assert(isequal(Hdup.r, H.r), 'duplicate flight sample changed the history');
Pr = P3;
Pr.flight.t([5 6]) = Pr.flight.t([6 5]);
Pcr2 = cell(1, 6);  Pcr2{k} = Pr;
try
    assemble_histories(M, Pcr2, asg, opts);
    error('test:noRefusal', 'non-monotone flight accepted');
catch err
    assert(strcmp(err.identifier, 'assemble_histories:flightGrid'), ...
           'flightGrid id: %s', err.identifier);
end

% ---- 4b. phase boundaries exactly on grid nodes (synthetic straight flight)
dt0 = 300/M.tStar;
tD0 = 10*dt0;  tA0 = 40*dt0;
rA = M.depState(mod(M.droPhase0(1) + tD0/M.tauD, 1));
rB = M.arrState(mod(M.slotPhase0(asg(1)) + tA0/M.tauT, 1));
fl = struct('t', (0:30).'*dt0, 'Y', zeros(31, 14));
fl.Y(:, 1:3) = rA(1:3).'.*(1 - (0:30).'/30) + rB(1:3).'.*((0:30).'/30);
Psyn = struct('ok', true, 'tDep', tD0, 'tArr', tA0, 'flight', fl);
Pcs = cell(1, 6);  Pcs{1} = Psyn;
oS = true(1, 6);  oS(1) = false;
Hs2 = assemble_histories(M, Pcs, asg, struct('onStation', oS, 'serviceDays', 1));
iD0 = find(abs(Hs2.t - tD0) < 1e-12);  iA0 = find(abs(Hs2.t - tA0) < 1e-12);
assert(iD0 == 11 && iA0 == 41, 'grid nodes of the events');
assert(Hs2.phase(iD0 - 1, 1) == 1 && Hs2.phase(iD0, 1) == 2 && ...
       Hs2.phase(iA0 - 1, 1) == 2 && Hs2.phase(iA0, 1) == 3, ...
       'phase at an event node');
fprintf('  event-node phases: %d %d | %d %d\n', Hs2.phase(iD0 - 1, 1), ...
        Hs2.phase(iD0, 1), Hs2.phase(iA0 - 1, 1), Hs2.phase(iA0, 1));

Ps = repmat(P3, 1, 6);
Hs = assemble_histories(M, Ps, asg, struct('onStation', onSt, 'serviceDays', 10));
assert(isequal(Hs.r, H.r), 'struct array form must equal the cell form');

% ---- 5. wall time: 50 days, 300 s, six spacecraft ----
t0 = tic;
Hw = assemble_histories(M, [], 1:6, struct('onSlots', true, ...
                        'tEnd', 50*86400/M.tStar));
fprintf('  50-day six-spacecraft onSlots: N = %d, wall %.2f s\n', ...
        numel(Hw.t), toc(t0));

fprintf('  PASS\n');
ok = true;
end
