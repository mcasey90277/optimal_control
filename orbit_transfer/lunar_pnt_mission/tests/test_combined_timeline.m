function ok = test_combined_timeline()
%% Purpose:
%
%   Verifies combined_timeline on SYNTHETIC leg-1 and leg-2 results files
%   (temporary files in a temporary folder, never in the repository):
%   (1) the bar ends of every satellite against hand-worked values (tStar
%   = 172800 s, so 1 ND = 2 d; leg-1 T0 = 1.5 ND; leg-1 satellite k departs
%   at 0.1k ND and arrives 0.8 ND later; leg-2 spacecraft k departs at 0.2k
%   ND and arrives 2 ND later on the leg-2 clock; the leg-2 service window
%   ends at 10 ND); (2) a figure is written when outPng is given and the
%   function closes the figure it opened; (3) Review Focus 4, DELIBERATE
%   FAILURES, each refused with its identifier: the leg-1 package identity
%   differs from the current package file (identity; accepted with
%   allowStalePackage), the leg-2 DRO differs (dro), the leg-2 droPhase0
%   differs (droPhase0), a leg-1 DRO period other than 1 (dro), the time
%   units differ (units), fewer than six accepted leg-1 transfers
%   (incomplete), T0 before the last leg-1 arrival (epoch), a leg-1
%   verdict other than PASS (verdict), a leg-1 transfer accepted with
%   G.pass false (incomplete), a leg-2 gate table with G0 FAIL or without
%   G0 (leg2); (4) then,
%   when both real results files exist, one run on them (compatibility
%   passes, six rows, leg-2 bars shifted by T0).
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
%  M. Casey  fix round 1: verdict, G.pass, leg-2 gate refusals  (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_lunar_pnt_paths.m'));
scratch = '/private/tmp/claude-501/-Users-msc-Desktop-optimal-control/70650963-d331-4c2e-9b00-6156c2bd5033/scratchpad/leg1';
if ~isfolder(scratch), mkdir(scratch); end
tmp = tempname(scratch);
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
nFig0 = numel(findall(groot, 'Type', 'figure'));

% ---- synthetic files ----
pkg = fullfile(tmp, 'pkg.mat');
pkgVar = magic(4);
save(pkg, 'pkgVar');
sha = file_sha(pkg);
md5Dro = 'aaaabbbbccccddddeeeeffff00001111';
tStar = 172800;
[l1, l2] = synthetic(tStar, sha, md5Dro);
f1 = fullfile(tmp, 'leg1.mat');  f2 = fullfile(tmp, 'leg2.mat');
save_leg1(f1, l1);  save_leg2(f2, l2);

% ---- 1. bar ends, hand-worked (days) ----
png = fullfile(tmp, 'timeline.png');
[T, info] = combined_timeline(f1, f2, struct('catalogMat', pkg, 'outPng', png));
assert(height(T) == 6, 'six rows expected, got %d', height(T));
for k = 1:6
    gw = [0, 0.2*k];
    l1b = [0.2*k, 0.2*k + 1.6];
    dc = [0.2*k + 1.6, 3];
    dw = [3, 3 + 0.4*k];
    l2b = [3 + 0.4*k, 3 + 0.4*k + 4];
    os = [3 + 0.4*k + 4, 23];
    got = [T.gtoWait(k, :); T.leg1(k, :); T.droCoast(k, :); T.droWait(k, :); T.leg2(k, :); T.onStation(k, :)];
    want = [gw; l1b; dc; dw; l2b; os];
    err = max(abs(got(:) - want(:)));
    assert(err < 1e-12, 'satellite %d: bar ends off by %.3g d', k, err);
    assert(T.sat(k) == k, 'row %d is satellite %d', k, T.sat(k));
end
fprintf('  synthetic: bar ends match the hand-worked values (T0 %.3f d, end %.3f d)\n', info.T0Days, info.endDays);
fprintf('    sat 3: GTO [0 0.6] leg1 [0.6 2.2] coast [2.2 3] wait [3 4.2] leg2 [4.2 8.2] on [8.2 23] d\n');
assert(abs(info.T0Days - 3) < 1e-12, 'T0 %.15g d, expected 3', info.T0Days);
assert(isfile(png), 'no figure written');
assert(numel(findall(groot, 'Type', 'figure')) == nFig0, 'combined_timeline left a figure open');
fprintf('  figure written and closed\n');

% ---- 3. refusals (Review Focus 4) ----
cases = {};
x1 = l1;  x1.run_meta.identity.sha256 = repmat('0', 1, 64);
cases(end+1, :) = {'package replaced', x1, l2, struct(), 'combined_timeline:identity'};
x2 = l2;  x2.fingerprint.droMd5 = 'ffffffffffffffffffffffffffffffff';
cases(end+1, :) = {'leg-2 DRO differs', l1, x2, struct(), 'combined_timeline:dro'};
x2 = l2;  x2.fingerprint.droPhase0 = (0:5)/6 + 0.01;
cases(end+1, :) = {'leg-2 droPhase0 differs', l1, x2, struct(), 'combined_timeline:droPhase0'};
x1 = l1;  x1.run_meta.tauD = 0.5;
cases(end+1, :) = {'leg-1 DRO period 0.5', x1, l2, struct(), 'combined_timeline:dro'};
x2 = l2;  x2.fingerprint.tStar = tStar + 1;
cases(end+1, :) = {'time units differ', l1, x2, struct(), 'combined_timeline:units'};
x1 = l1;  x1.P1(4).accepted = false;
cases(end+1, :) = {'satellite 4 not accepted', x1, l2, struct(), 'combined_timeline:incomplete'};
x1 = l1;  x1.T0 = 1.35;
cases(end+1, :) = {'T0 before the last leg-1 arrival', x1, l2, struct(), 'combined_timeline:epoch'};
x1 = l1;  x1.run_meta.verdict = 'not passed: GH';
cases(end+1, :) = {'leg-1 verdict not PASS', x1, l2, struct(), 'combined_timeline:verdict'};
x1 = l1;  x1.P1(2).G.pass = false;
cases(end+1, :) = {'satellite 2 accepted but G.pass false', x1, l2, struct(), 'combined_timeline:incomplete'};
x2 = l2;  x2.gates(2, 2) = {'FAIL'};
cases(end+1, :) = {'leg-2 validation G0 FAIL', l1, x2, struct(), 'combined_timeline:leg2'};
x2 = l2;  x2.gates = x2.gates(1, :);
cases(end+1, :) = {'leg-2 gate table without G0', l1, x2, struct(), 'combined_timeline:leg2'};
for kc = 1:size(cases, 1)
    g1 = fullfile(tmp, sprintf('c%d_leg1.mat', kc));  g2 = fullfile(tmp, sprintf('c%d_leg2.mat', kc));
    save_leg1(g1, cases{kc, 2});  save_leg2(g2, cases{kc, 3});
    o = cases{kc, 4};  o.catalogMat = pkg;
    id = '';
    try
        combined_timeline(g1, g2, o);
    catch err
        id = err.identifier;
        msg = err.message;
    end
    assert(strcmp(id, cases{kc, 5}), '%s: expected %s, got "%s"', cases{kc, 1}, cases{kc, 5}, id);
    fprintf('  refused (%s): %s\n', cases{kc, 1}, msg);
end
x1 = l1;  x1.run_meta.identity.sha256 = repmat('0', 1, 64);
g1 = fullfile(tmp, 'stale_leg1.mat');  save_leg1(g1, x1);
[Ts, infoS] = combined_timeline(g1, f2, struct('catalogMat', pkg, 'allowStalePackage', true));
assert(height(Ts) == 6 && ~infoS.packageMatch, 'allowStalePackage should build and report the mismatch');
fprintf('  allowStalePackage: built, packageMatch = %d\n', infoS.packageMatch);

% ---- 4. the real files, when present ----
r1 = fullfile(root, 'results', 'lunar_pnt_leg1_study.mat');
r2 = fullfile(root, 'results', 'lunar_pnt_mission_study.mat');
if isfile(r1) && isfile(r2)
    [Tr, infoR] = combined_timeline(r1, r2, struct());
    R1 = load(r1, 'T0', 'P1');
    R2 = load(r2, 'P');
    dayNd = infoR.tStar/86400;
    e1 = max(abs(Tr.leg2(:, 1) - (R1.T0 + [R2.P.tDep].')*dayNd));
    e2 = max(abs(Tr.leg1(:, 2) - [R1.P1.tArr].'*dayNd));
    assert(height(Tr) == 6 && e1 < 1e-12 && e2 < 1e-12, 'real files: shift error %.3g / %.3g d', e1, e2);
    fprintf('  real files: compatible; T0 %.3f d; leg-2 departures %.3f .. %.3f d; on station from %.3f .. %.3f d\n', ...
            infoR.T0Days, min(Tr.leg2(:, 1)), max(Tr.leg2(:, 1)), min(Tr.onStation(:, 1)), max(Tr.onStation(:, 1)));
else
    fprintf('  real files: not both present, skipped\n');
end
fprintf('  PASS\n');
ok = true;
end

% ------------------------------------------------------------------------
function [l1, l2] = synthetic(tStar, sha, md5Dro)
% SYNTHETIC  Minimal leg-1 and leg-2 results structs (see the header).
%   INPUTS: tStar [s]; sha char [1 x 64]; md5Dro char [1 x 32].
%   OUTPUTS: l1, l2 structs (variables of the two files).
P1 = struct('tDep', num2cell(0.1*(1:6)), 'tArr', num2cell(0.1*(1:6) + 0.8), 'accepted', true, ...
            'G', struct('pass', true));
l1.P1 = P1;
l1.T0 = 1.5;
l1.run_meta = struct('identity', struct('sha256', sha, 'file', 'pkg.mat'), 'droMd5', md5Dro, ...
                     'droPhase0', (0:5)/6, 'tauD', 1, 'regime', struct('tStar', tStar), 'verdict', 'PASS');
l2.P = struct('tDep', num2cell(0.2*(1:6)), 'tArr', num2cell(0.2*(1:6) + 2));
l2.fingerprint = struct('droMd5', md5Dro, 'droPhase0', (0:5)/6, 'tStar', tStar);
l2.H = struct('tEnd', 10);
l2.gates = {'I1', 'PASS', 1, 'synthetic'; 'G0', 'PASS', 0, 'synthetic'};
end

% ------------------------------------------------------------------------
function save_leg1(f, s)
% SAVE_LEG1  Write the leg-1 variables.  INPUTS: f char; s struct.
P1 = s.P1;  T0 = s.T0;  run_meta = s.run_meta;
save(f, 'P1', 'T0', 'run_meta');
end

% ------------------------------------------------------------------------
function save_leg2(f, s)
% SAVE_LEG2  Write the leg-2 variables.  INPUTS: f char; s struct.
P = s.P;  fingerprint = s.fingerprint;  H = s.H;  gates = s.gates;
save(f, 'P', 'fingerprint', 'H', 'gates');
end

% ------------------------------------------------------------------------
function h = file_sha(f)
% FILE_SHA  Hex SHA-256 of a file's bytes.  INPUTS: f char.  OUTPUTS: h char.
md  = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(f, 'r');
raw = fread(fid, Inf, '*uint8');
fclose(fid);
md.update(raw);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
