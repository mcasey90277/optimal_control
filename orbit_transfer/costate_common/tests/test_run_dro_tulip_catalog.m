function ok = test_run_dro_tulip_catalog()
%% Purpose:
%
%   Tests that run_dro_tulip is served by the CERTIFIED 70 mN CATALOG, not
%   only by the ten solutions dro_tulip_library gathered from the anchor and
%   sweep files. Before 2026-09-11 the front door knew 10 of the catalog's
%   115 entries, so asking for any of the other 105 phase pairs started a
%   continuation walk -- minutes to hours -- to re-derive a transfer that was
%   already certified, audited and on disk.
%
%   The catalog stores z8 but no junction states, so the front door builds
%   the multiple-shooting seed from z8 on demand (seed_from_z8) and sends it
%   through the same gate stack as every other route.
%
%   Checks:
%     1. dro_tulip_library lists the catalog's entries when asked
%        (includeCatalog), and still lists only the file-backed ten by
%        default -- build_arrival_sheet seeds from the default list, and
%        seeding a sheet with its own previous output would be circular;
%     2. a catalog-only pair is served from the library, not walked, and
%        certifies at the catalog's own t_f.
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
ind = fullfile(fileparts(here), 'DRO_tulip', 'indirect');
addpath(here, ind);
catMat = fullfile(ind, 'results', 'costate_catalog_dro_tulip_70mN.mat');
if ~isfile(catMat), fprintf('  SKIP  no 70 mN catalog on disk\n'); return, end
L = load(catMat);  fn = fieldnames(L);  s = L.(fn{1}).sheets(1);
tS = 382981.289129055/86400;

% the test pair: row 2 (sD = 1/12), column 3 -- certified in the catalog,
% absent from the file-backed library
iD = 2;  iA = 3;
assert(s.has_solution(iD, iA, 1), 'test pair (%d,%d) is not in the catalog', iD, iA);
sD = s.sD_frac(iD);  sA = s.sA_frac(iA);
z8cat = s.z8(:, s.entry_index(iD, iA, 1));

% 1. listing
try
    libDef = dro_tulip_library(ind);
    libCat = dro_tulip_library(ind, struct('includeCatalog', true));
catch ME
    fprintf('  FAIL  dro_tulip_library does not take includeCatalog (%s)\n', ME.message);
    fprintf('TEST_RUN_DRO_TULIP_CATALOG: FAIL\n');  ok = false;  return
end
hit = @(lib) find(abs([lib.sD] - sD) < 1e-9 & abs([lib.sA] - sA) < 1e-9, 1);
ok = chk(ok, isempty(hit(libDef)), 'the default listing does not include catalog entries');
k = hit(libCat);
ok = chk(ok, ~isempty(k) && strcmp(libCat(k).src, 'catalog') && isequal(libCat(k).z(:), z8cat), ...
         'includeCatalog lists the pair, sourced from the catalog, with its z8');
ok = chk(ok, numel(libCat) >= nnz(s.has_solution), ...
         sprintf('includeCatalog lists at least the catalog''s %d entries (got %d)', nnz(s.has_solution), numel(libCat)));
if ~ok, fprintf('TEST_RUN_DRO_TULIP_CATALOG: FAIL\n'); return, end

% 2. the front door serves it from the library and certifies it
t0 = tic;
T = run_dro_tulip(sD, sA, struct('quiet', true, 'wallSec', 120));
ok = chk(ok, strcmp(T.source, 'library'), sprintf('served from the library (source: %s, %.0f s)', T.source, toc(t0)));
ok = chk(ok, T.ok, sprintf('certified (%s)', T.reason));
ok = chk(ok, abs(T.tfDays - z8cat(8)*tS) < 1e-6, ...
         sprintf('t_f = %.6f d against the catalog''s %.6f d', T.tfDays, z8cat(8)*tS));

if ok, fprintf('TEST_RUN_DRO_TULIP_CATALOG: ALL PASS\n'); else, fprintf('TEST_RUN_DRO_TULIP_CATALOG: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
