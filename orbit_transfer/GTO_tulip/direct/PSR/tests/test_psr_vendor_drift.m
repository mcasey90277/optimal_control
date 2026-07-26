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
% This test is that comparison. It does NOT consolidate anything and does not
% touch the path: it reads files by absolute path, so PSR's self-containment
% guarantee is untouched.
%
% All 19 .m files in PSR/lib are copies of something. This test asserts:
%   * the 15 that are byte-identical to their origin STAY byte-identical;
%   * the 4 known-divergent ones are still divergent for the recorded reason,
%     and it names them, so "expected drift" can never quietly become
%     "unnoticed drift".
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
};

% name -> {origin, why it is allowed to differ}
KNOWN_DIVERGENT = {
    'casadi_minfuel_sundman', fullfile(tulipD,'sundman_minfuel'), ...
        'upstream gained returnModel/creg, vBox/rBox, boundSat; PSR frozen at the 2026-07-12 snapshot'
    'minfuel_at_tf',          fullfile(tulipD,'sundman_minfuel'), ...
        'PSR-owned variant'
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

if nBad > 0
    error('test_psr_vendor_drift:drift', ...
        'PSR/lib vendor drift, %d issue(s):\n  - %s\n', nBad, strjoin(msgs, '\n  - '));
end
fprintf('\ntest_psr_vendor_drift: ALL PASS (%d identical, %d known-divergent)\n', ...
        size(IDENTICAL,1), size(KNOWN_DIVERGENT,1));
