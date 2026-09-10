function ok = test_deliverable_audit_gate()
%% Purpose:
%
%   Tests the AUDIT GATE on build_dro_deliverable: a catalog does not ship
%   until an audit that re-derived it from its own keys comes back clean.
%
%   Fixing a builder says nothing about a file already shipped, which is why
%   this is enforced in code rather than written in a checklist. The gate
%   must refuse three ways: no audit at all, an audit of a DIFFERENT
%   catalog, and an audit with findings. It may be waived only by an
%   explicit, named opt-out.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip'), fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
catMat = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', ...
                  'costate_catalog_dro_tulip_70mN.mat');
if ~isfile(catMat)
    fprintf('  SKIP  the 70 mN catalog is not on disk\n');  return
end
outDir = fullfile(tempdir, 'deliv_gate');

% (1) no audit -> refused
ok = chk(ok, refuses(catMat, outDir, struct()), 'no audit: refused');

% (2) an audit of a DIFFERENT catalog -> refused
A = struct('nBad', 0, 'nOk', 1, 'catMat', 'some/other/catalog.mat', 'rows', struct('ok', true));
m1 = fullfile(tempdir, 'audit_wrong.mat');  save(m1, 'A');
ok = chk(ok, refuses(catMat, outDir, struct('auditMat', m1)), 'audit of another catalog: refused');

% (3) an audit WITH findings -> refused
A = struct('nBad', 2, 'nOk', 51, 'catMat', catMat, 'rows', struct('ok', false));
m2 = fullfile(tempdir, 'audit_bad.mat');  save(m2, 'A');
ok = chk(ok, refuses(catMat, outDir, struct('auditMat', m2)), 'audit with findings: refused');

% (4) the real, clean audit -> ships
real_ = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', 'audit_70mN.mat');
if isfile(real_)
    try
        out = build_dro_deliverable(struct('catMat', catMat, 'outDir', outDir, ...
                                           'zip', false, 'auditMat', real_));
        ok = chk(ok, ~isempty(out.dir) && any(contains(out.files, 'README.md')), ...
                 sprintf('clean audit: ships (%d files)', numel(out.files)));
    catch ME
        ok = chk(ok, false, ['clean audit should ship but threw: ' ME.message]);
    end
else
    fprintf('  SKIP  no real audit on disk\n');
end

% (5) explicit named opt-out is honoured
try
    build_dro_deliverable(struct('catMat', catMat, 'outDir', outDir, 'zip', false, ...
                                 'skipAuditBecause', 'test'));
    ok = chk(ok, true, 'explicit named opt-out is honoured');
catch ME
    ok = chk(ok, false, ['opt-out should work but threw: ' ME.message]);
end

if ok, fprintf('TEST_DELIVERABLE_AUDIT_GATE: ALL PASS\n'); else, fprintf('TEST_DELIVERABLE_AUDIT_GATE: FAIL\n'); end
end

function tf = refuses(catMat, outDir, extra)
% REFUSES  True when the builder refuses to package.  INPUTS: catMat;
% outDir; extra.  OUTPUTS: tf.
cfg = extra;  cfg.catMat = catMat;  cfg.outDir = outDir;  cfg.zip = false;
tf = false;
try
    build_dro_deliverable(cfg);
catch
    tf = true;
end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
