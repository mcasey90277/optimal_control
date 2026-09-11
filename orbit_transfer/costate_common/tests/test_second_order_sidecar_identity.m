function ok = test_second_order_sidecar_identity()
%% Purpose:
%
%   Tests that second_order_pass binds every sidecar record to the catalog
%   entry it measured. The sidecar used to be POSITIONAL: record q was
%   whatever entry sat q-th in find(has_solution) at write time, so a
%   re-packaged catalog with a different entry set would have silently
%   inherited measurements of other entries on resume, and the chain
%   script pointed at the FIRST sweep's sidecar, whose lift margins differ
%   from the catalog's by up to 1.35e4 (found 2026-09-11).
%
%   Checks, all with maxEntries = 0 so no instrument runs:
%     1. an unkeyed (legacy) sidecar is refused by default;
%     2. adoptLegacy accepts one only when every written-back value in the
%        catalog equals the sidecar's, and then stamps the keys;
%     3. the first-sweep sidecar, which does not match, is refused even
%        with adoptLegacy;
%     4. a keyed record that points at a different entry is refused;
%     5. a keyed record whose stored z8 differs from the catalog's is
%        refused (same cell, different solution).
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
res = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results');
catMat = fullfile(res, 'costate_catalog_dro_tulip_70mN.mat');
v2 = fullfile(res, 'second_order_progress_v2.mat');
v1 = fullfile(res, 'second_order_progress.mat');
if ~(isfile(catMat) && isfile(v2) && isfile(v1))
    fprintf('  SKIP  70 mN catalog or its sidecars not on disk\n');  return
end
tmp = tempname;  mkdir(tmp);
q0 = struct('maxEntries', 0, 'logFile', '');

% 1. unkeyed sidecar refused by default
f = fullfile(tmp, 'legacy.mat');  copyfile(v2, f);
ok = chk(ok, refuses(@() second_order_pass(catMat, setfield(q0, 'sideMat', f)), ...
         'second_order_pass:unkeyedSidecar'), 'an unkeyed sidecar is refused by default');

% 2. adoptLegacy on the matching sidecar: accepted, keys stamped
o = q0;  o.sideMat = f;  o.adoptLegacy = true;
S = second_order_pass(catMat, o);
P = load(f);
ok = chk(ok, S.done && isfield(P.R, 'key') && isfield(P.R, 'z8') && ...
         all(arrayfun(@(r) numel(r.key) == 3 && numel(r.z8) == 8, P.R)), ...
         'the matching sidecar is adopted and every record is keyed');

% 3. the first-sweep sidecar does not match the catalog: refused
g = fullfile(tmp, 'v1.mat');  copyfile(v1, g);
o = q0;  o.sideMat = g;  o.adoptLegacy = true;
ok = chk(ok, refuses(@() second_order_pass(catMat, o), 'second_order_pass:legacyMismatch'), ...
         'a legacy sidecar that disagrees with the catalog is refused');

% 4. a keyed record pointing at another entry: refused
R = P.R;  R([1 2]) = R([2 1]);  h = fullfile(tmp, 'swapped.mat');  save(h, 'R');
ok = chk(ok, refuses(@() second_order_pass(catMat, setfield(q0, 'sideMat', h)), ...
         'second_order_pass:staleSidecar'), 'a record keyed to a different entry is refused');

% 5. same key, different solution: refused
R = P.R;  R(1).z8 = R(1).z8*(1 + 1e-6);  h = fullfile(tmp, 'moved.mat');  save(h, 'R');
ok = chk(ok, refuses(@() second_order_pass(catMat, setfield(q0, 'sideMat', h)), ...
         'second_order_pass:staleSidecar'), 'a record whose z8 differs from the catalog is refused');

rmdir(tmp, 's');
if ok, fprintf('TEST_SECOND_ORDER_SIDECAR_IDENTITY: ALL PASS\n');
else,  fprintf('TEST_SECOND_ORDER_SIDECAR_IDENTITY: FAIL\n'); end
end

function r = refuses(fh, id)
% REFUSES  True when fh throws the named identifier.  INPUTS: fh; id.
% OUTPUTS: r [logical].
try
    fh();  r = false;
catch ME
    r = strcmp(ME.identifier, id);
    if ~r, fprintf('        (threw %s: %s)\n', ME.identifier, ME.message); end
end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
