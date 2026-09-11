function ok = test_guard_catalog_overwrite()
%% Purpose:
%
%   Tests guard_catalog_overwrite, the check build_70mN_library runs before
%   its PACKAGE stage. Packaging rewrites the catalog file, and the
%   second-order sweep's measurements live only in that file (catalogs are
%   gitignored; the .bak_2nd backup predates the writeback). With the
%   shipped switches -- package on, sweep off -- a plain run would have
%   erased them with no way back (found 2026-09-11).
%
%   Checks:
%     1. no catalog on disk: allowed, no backup;
%     2. a catalog WITHOUT second-order fields: allowed, backed up;
%     3. a catalog WITH them and the sweep stage off: refused, file intact;
%     4. the same with the sweep stage on: allowed, backed up byte-equal.
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
addpath(fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
tmp = tempname;  mkdir(tmp);
f = fullfile(tmp, 'costate_catalog_test.mat');

% 1. nothing to overwrite
bak = guard_catalog_overwrite(f, false);
ok = chk(ok, isempty(bak), 'no catalog on disk: allowed, nothing backed up');

% 2. a plain catalog
costate_catalog_test = struct('n_entries', 3, 'sheets', struct('z8', zeros(8,3)));
save(f, 'costate_catalog_test');
bak = guard_catalog_overwrite(f, false);
ok = chk(ok, ~isempty(bak) && isfile(bak), 'a catalog without second-order fields: allowed and backed up');

% 3. carrying the sweep's writeback, sweep off
costate_catalog_test.second_order = struct('date', '2026-09-11');
save(f, 'costate_catalog_test');
before = fileread_bytes(f);
ok = chk(ok, refuses(@() guard_catalog_overwrite(f, false), ...
         'guard_catalog_overwrite:wouldStripSecondOrder'), ...
         'second-order fields with the sweep off: refused');
ok = chk(ok, isequal(before, fileread_bytes(f)), 'and the catalog is untouched');

% 4. sweep on: allowed, backup byte-equal to the original
bak = guard_catalog_overwrite(f, true);
ok = chk(ok, ~isempty(bak) && isequal(fileread_bytes(bak), before), ...
         'with the sweep on: allowed, backup byte-equal');

rmdir(tmp, 's');
if ok, fprintf('TEST_GUARD_CATALOG_OVERWRITE: ALL PASS\n');
else,  fprintf('TEST_GUARD_CATALOG_OVERWRITE: FAIL\n'); end
end

function b = fileread_bytes(f)
% FILEREAD_BYTES  Raw bytes of a file.  INPUTS: f.  OUTPUTS: b [uint8].
fid = fopen(f, 'r');  b = fread(fid, inf, '*uint8');  fclose(fid);
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
