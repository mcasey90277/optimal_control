function ok = test_run_costate_library_seams()
% TEST_RUN_COSTATE_LIBRARY_SEAMS  Two P0 repairs to the front door (FINDINGS
% 78): (1) the sheet builder receives the campaign's departure phase and its
% anchor, not the shipped defaults; (2) 'packaged' is no longer the only
% success signal -- package, audit and sweep are reported separately in
% out.stages so a caller can require all three.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'), ...
        fullfile(fileparts(fileparts(here)), 'campaign_common'));
tmp = fullfile(tempdir, sprintf('rcl_seams_%s', char(java.util.UUID.randomUUID())));
mkdir(tmp);  cleanup = onCleanup(@() rmdir(tmp, 's'));

% ---- 1. the sheet builder sees sD(1), the anchor file and its phase -------
sD = [0.2 0.5];  sA = [0.1 0.6];
recF = fullfile(tmp, 'received.mat');
engine = struct('thrustN', 0.07, 'ispS', 900, 'm0kg', 150);
orbits = struct('tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1);
fake = @(o) fakeSheet(o, recF, sA, engine, orbits);
o = run_costate_library(struct('outDir', tmp, 'sD', sD, 'sA', sA, 'launch', false, ...
        'anchorMat', '/x/anchor.mat', 'anchorSA', 0.61, 'sheetBuilder', fake, 'librarySeeds', false, ...
        'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false)));
ok = chk(ok, isfile(recF), 'the injected sheet builder was called');
if isfile(recF)
    R = load(recF);  r = R.received;
    ok = chk(ok, isfield(r, 'sD') && r.sD == 0.2, 'sheet builder receives sD = sD(1), not the default 0');
    ok = chk(ok, isfield(r, 'anchorMat') && strcmp(r.anchorMat, '/x/anchor.mat'), 'sheet builder receives the campaign anchor file');
    ok = chk(ok, isfield(r, 'sA0') && r.sA0 == 0.61, 'sheet builder receives the anchor''s arrival phase');
end
ok = chk(ok, strcmp(o.state, 'blocked'), 'an uncertified fake sheet leaves the call blocked (nothing else ran)');

% ---- 2. stage outcomes are reported separately -----------------------------
H = run_costate_library('localfunctions');
on = @(f) true;
s = H.stageOutcomes({'audit: 3 of 576 entries bad', '2 column(s) are shorter than 23 points: ...'}, on);
ok = chk(ok, ~s.audit && s.sweep && s.package, 'stageOutcomes: an audit blocker fails audit only');
s = H.stageOutcomes({'sweep incomplete: 4 entries still to measure'}, on);
ok = chk(ok, ~s.sweep && s.audit, 'stageOutcomes: a sweep blocker fails sweep only');
s = H.stageOutcomes({'sweep threw: boom'}, on);
ok = chk(ok, ~s.sweep, 'stageOutcomes: a sweep that threw fails sweep');
s = H.stageOutcomes({}, on);
ok = chk(ok, s.package && s.audit && s.sweep, 'stageOutcomes: no blockers, every stage ok');
s = H.stageOutcomes({'audit never run'}, @(f) ~strcmp(f, 'audit'));
ok = chk(ok, isnan(s.audit) || islogical(s.audit), 'stageOutcomes: a stage that was off is reported, not failed');

if ok, fprintf('test_run_costate_library_seams: ALL PASS\n'); else, fprintf('test_run_costate_library_seams: FAIL\n'); end
end

function fakeSheet(o, recF, sA, engine, orbits)
% FAKESHEET  Stand-in for build_arrival_sheet: record what it was given and
% write a sheet with no certified column, carrying the problem identity the
% front door checks.  INPUTS: o (the builder's options); recF; sA; engine;
% orbits.  OUTPUTS: none (files).
received = o; %#ok<NASGU>
save(recF, 'received');
nA = numel(sA);
problem = struct('thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
                 'tauDRO', orbits.tauDRO, 'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, 'sD', o.sD);
S = struct('sA', sA, 'TF', nan(1, nA), 'Z8', nan(8, nA), 'cand', {cell(1, nA)}, 'problem', problem); %#ok<NASGU>
save(o.out, 'S');
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
