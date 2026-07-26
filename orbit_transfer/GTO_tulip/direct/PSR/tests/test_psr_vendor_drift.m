% TEST_PSR_VENDOR_DRIFT  Guard PSR/lib's vendored copies against silent drift.
%
% WHY THIS EXISTS. PSR is self-contained by design: PSR/lib/README.md records
% that run_psr reaches only PSR/, PSR/lib/ and pumpkyn, verified with
% requiredFilesAndProducts under restoredefaultpath. That isolation works
% exactly as intended -- which is the problem it creates: upstream fixes never
% arrive, and nothing notices.
%
% It has already cost us once. casadi_minfuel_sundman drifted 4.7 kB apart from
% its origin (16.6 vs 11.9 kB) as returnModel/creg, vBox/rBox and boundSat were
% added upstream. run_psr's stage 5c then called a capability the vendored copy
% did not have, threw, and was swallowed by a catch -- an advisory warning on
% every run, with the gate never once executing, undetected for a day.
%
% The defect was never the vendoring. It was that no test compared the copies.
% This test is that comparison. It consolidates nothing; the drift half reads
% files by absolute path, and the resolution half saves and restores the path
% exactly, so PSR's isolation is unaffected either way.
%
% All 20 .m files in PSR/lib are copies of something. This test asserts:
%   * the 16 that are byte-identical to their origin STAY byte-identical;
%   * the 4 known-divergent ones are still divergent for the recorded reason,
%     and it names them, so "expected drift" can never quietly become
%     "unnoticed drift";
%   * the PATH RESOLUTION of the shadowed pair does not silently flip -- see
%     the second half. Byte-identity alone is not enough: two of these files
%     are dead code, and a one-line addpath change would bring them to life.
%
% If this test fails, do NOT just re-copy the file. Decide deliberately:
% either port the upstream change into PSR (and re-run the certified PSR
% result), or record the divergence as intentional in PSR/lib/README.md and
% move the file into the KNOWN_DIVERGENT list below with its reason.
%
% INPUTS:  none
% OUTPUTS: none (prints a per-file table; throws on unexpected drift)
%
% REFERENCES:
%   [1] GTO_tulip/direct/PSR/lib/README.md (the vendoring contract + manifest).
%   [2] GTO_tulip/doc/EXECUTION_PATHS.md sec 1 (the measured duplicate table).
%   [3] orbit_transfer/CODE_STRUCTURE.md (why a drift test rather than a
%       consolidation: the vendoring is deliberate; the silence was the bug).
here    = fileparts(mfilename('fullpath'));          % .../PSR/tests
psrLib  = fullfile(here, '..', 'lib');
tulipD  = fullfile(here, '..', '..');                % .../direct
common  = fullfile(here, '..', '..', '..', '..', 'cr3bp_common');
msBand  = fullfile(here, '..', '..', '..', 'indirect', 'ms_band');
refine  = fullfile(tulipD, 'sundman_minfuel', 'refine');

% name -> origin folder. Provenance, not just a checksum list: if an origin
% moves, this test fails loudly rather than silently skipping the file.
IDENTICAL = {
    'refine_loop',           refine
    'refine_sigma',          refine
    'prep_refine_seed',      refine
    'pmp_refine_indicator',  refine
    'warmstart_on_mesh',     refine
    'sms_eom',               msBand
    'sms_jacobian_cs',       msBand
    'sms_pack',              msBand
    'sms_problem',           msBand
    'sms_residual',          msBand
    'sms_seed_duals',        msBand
    'sms_unpack',            msBand
    'verify_direct_pmp',     msBand
    'beta_from_duals',       msBand
    'gto_tulip_endpoints',   common
    'insertion_states',      common
};

% name -> {origin, why it is allowed to differ}
KNOWN_DIVERGENT = {
    'casadi_minfuel_sundman', fullfile(tulipD,'sundman_minfuel'), ...
        ['DEAD CODE -- shadowed by the entry-point addpath(../sundman_minfuel), ' ...
         'so PSR runs the UPSTREAM copy. Divergence is therefore harmless today ' ...
         'AND un-shadowing it would change PSR''s solver. See PSR/lib/README.md']
    'minfuel_at_tf',          fullfile(tulipD,'sundman_minfuel'), ...
        'DEAD CODE -- shadowed by the same addpath; PSR runs the upstream copy'
    'cr3bp_lt_params',        common, ...
        'PSR-owned variant'
    'minfuel_config',         common, ...
        'documented in PSR/lib/README.md: dirs repointed to ../../sundman_minfuel/results'
};

% Direct content comparison -- no hashing dependency, and an exact answer.
sameFile = @(a,b) isequal(fileread(a), fileread(b));

fprintf('\n=== PSR/lib vendor-drift check ===\n');
nBad = 0;  msgs = {};

for k = 1:size(IDENTICAL,1)
    nm = IDENTICAL{k,1};  org = IDENTICAL{k,2};
    pv = fullfile(psrLib, [nm '.m']);   po = fullfile(org, [nm '.m']);
    assert(isfile(pv), 'test_psr_vendor_drift: vendored file missing: %s', pv);
    if ~isfile(po)
        nBad = nBad + 1;
        msgs{end+1} = sprintf(['%s: ORIGIN MISSING at %s -- the source moved or was ' ...
            'deleted. Update this test''s provenance map deliberately.'], nm, po); %#ok<AGROW>
        fprintf('  %-24s ORIGIN MISSING\n', nm);  continue
    end
    same = sameFile(pv, po);
    if same
        fprintf('  %-24s identical\n', nm);
    else
        fprintf('  %-24s DRIFTED  <--\n', nm);
    end
    if ~same
        nBad = nBad + 1;
        msgs{end+1} = sprintf(['%s DRIFTED from its origin (%s). Do NOT simply re-copy: ' ...
            'either port the upstream change into PSR and re-run the certified PSR ' ...
            'result, or move this entry to KNOWN_DIVERGENT with the reason.'], ...
            nm, po); %#ok<AGROW>
    end
end

fprintf('\n  known-divergent (asserted STILL divergent, so expected drift\n');
fprintf('  can never quietly become unnoticed drift):\n');
for k = 1:size(KNOWN_DIVERGENT,1)
    nm = KNOWN_DIVERGENT{k,1};  org = KNOWN_DIVERGENT{k,2};  why = KNOWN_DIVERGENT{k,3};
    pv = fullfile(psrLib, [nm '.m']);  po = fullfile(org, [nm '.m']);
    assert(isfile(pv) && isfile(po), 'test_psr_vendor_drift: missing %s or %s', pv, po);
    if sameFile(pv, po)
        nBad = nBad + 1;
        msgs{end+1} = sprintf(['%s is now IDENTICAL to its origin but is listed as ' ...
            'known-divergent. The divergence was resolved -- move it to the ' ...
            'IDENTICAL list so it is guarded from now on.'], nm); %#ok<AGROW>
        fprintf('  %-24s converged (was divergent) <-- update the list\n', nm);
    else
        fprintf('  %-24s divergent as expected  (%s)\n', nm, why);
    end
end

%% ------------------------------------------------------------------------
%% Path resolution: which copy actually EXECUTES
%% ------------------------------------------------------------------------
% Byte-identity is only half the guarantee. The other half is which copy the
% path picks -- and on 2026-07-15 that silently flipped for two files when the
% insertion-points feature added addpath(../sundman_minfuel) to every PSR entry
% point. addpath PREPENDS, so the upstream folder went in front of PSR/lib.
%
% This half pins the resolution that real runs get, so a future addpath edit
% cannot quietly swap PSR's solver in either direction. Both outcomes are
% failures worth catching:
%   * upstream -> PSR/lib would swap the certified pipeline onto the frozen
%     2026-07-12 solver snapshot;
%   * PSR/lib -> upstream for a file expected local would bypass a PSR-owned
%     variant (refine_loop, verify_direct_pmp, minfuel_config).
% The path/cwd juggling lives in a local FUNCTION, not inline: its onCleanup
% then fires on return -- including on error -- whereas an onCleanup created in
% a script's workspace survives until the variable is cleared, and would leave
% the caller in PSR/ with a mutated path.
upstream = fullfile(tulipD, 'sundman_minfuel');
% name -> folder the ENTRY-POINT path state must resolve it to
EXPECT_RESOLUTION = {
    'casadi_minfuel_sundman', upstream,  'shadowed: PSR runs the upstream solver'
    'minfuel_at_tf',          upstream,  'shadowed: PSR runs the upstream driver'
    'insertion_states',       psrLib,    'vendored here; upstream no longer defines it'
    'refine_loop',            psrLib,    'PSR-owned variant must win'
    'verify_direct_pmp',      psrLib,    'PSR-owned variant must win'
    'minfuel_config',         psrLib,    'PSR-owned variant (repointed dirs) must win'
};

fprintf('\n  path resolution under a real entry-point path state:\n');
resolved = resolveUnderEntryPointPath(fullfile(here, '..'), EXPECT_RESOLUTION(:,1));

for k = 1:size(EXPECT_RESOLUTION,1)
    nm = EXPECT_RESOLUTION{k,1};  wantDir = EXPECT_RESOLUTION{k,2};  why = EXPECT_RESOLUTION{k,3};
    got = resolved{k};
    if isempty(got)
        nBad = nBad + 1;
        msgs{end+1} = sprintf('%s does not resolve at all under PSR''s entry-point path.', nm); %#ok<AGROW>
        fprintf('    %-24s NOT FOUND <--\n', nm);  continue
    end
    % Compare canonicalised folders, not strings: '..' segments differ by caller.
    gotDir = canon(fileparts(got));
    if strcmp(gotDir, canon(wantDir))
        fprintf('    %-24s %-9s (%s)\n', nm, folderTag(gotDir, psrLib), why);
    else
        nBad = nBad + 1;
        fprintf('    %-24s RESOLUTION CHANGED <--\n', nm);
        msgs{end+1} = sprintf(['%s now resolves to %s but must resolve to %s (%s). ' ...
            'An addpath change has swapped which copy executes -- this changes ' ...
            'behaviour on a certified pipeline. Do not "fix" the test; decide ' ...
            'deliberately and re-run the certified PSR result.'], ...
            nm, gotDir, canon(wantDir), why); %#ok<AGROW>
    end
end

if nBad > 0
    error('test_psr_vendor_drift:drift', ...
        'PSR/lib vendor drift, %d issue(s):\n  - %s\n', nBad, strjoin(msgs, '\n  - '));
end
fprintf('\ntest_psr_vendor_drift: ALL PASS (%d identical, %d known-divergent, %d resolutions pinned)\n', ...
        size(IDENTICAL,1), size(KNOWN_DIVERGENT,1), size(EXPECT_RESOLUTION,1));

function c = canon(p)
% CANON  Resolve '..' segments so folder comparisons are textual-safe.
% INPUTS:  p - folder path [char]
% OUTPUTS: c - canonical absolute path [char]
d = dir(p);
if isempty(d), c = p; else, c = d(1).folder; end
end

function t = folderTag(d, psrLib)
% FOLDERTAG  Short label for which side a resolution landed on.
% INPUTS:  d - canonical folder [char]; psrLib - PSR/lib folder [char]
% OUTPUTS: t - 'PSR/lib' or 'upstream' [char]
if strcmp(d, canon(psrLib)), t = 'PSR/lib'; else, t = 'upstream'; end
end

function locs = resolveUnderEntryPointPath(psrRoot, names)
% RESOLVEUNDERENTRYPOINTPATH  which() each name under PSR's real entry-point
% path state, then restore the caller's path and cwd exactly.
%
% Reproduces what run_psr / psr_run_one / gen_energy_seed do: setup_paths
% followed by addpath(../sundman_minfuel). Being a function, its onCleanup
% fires on return or on error, so the caller never inherits the mutated state.
%
% INPUTS:
%   psrRoot - the PSR/ directory [char]
%   names   - function names to resolve [Nx1 cell of char]
%
% OUTPUTS:
%   locs    - full path of each resolved file, '' if unresolved [Nx1 cell]
oldPath = path();  oldDir = pwd;
cleanup = onCleanup(@() restoreState(oldPath, oldDir)); %#ok<NASGU>

cd(psrRoot);
addpath(psrRoot);  setup_paths();
addpath(fullfile(psrRoot, '..', 'sundman_minfuel'));   % as every entry point does

locs = cell(numel(names), 1);
for k = 1:numel(names), locs{k} = which(names{k}); end
end

function restoreState(oldPath, oldDir)
% RESTORESTATE  Restore the caller's path and cwd exactly.
% INPUTS:  oldPath - path string to restore [char]; oldDir - cwd [char]
% OUTPUTS: none
path(oldPath);  cd(oldDir);
end
