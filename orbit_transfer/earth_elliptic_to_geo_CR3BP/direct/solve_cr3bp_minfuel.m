function best = solve_cr3bp_minfuel(opts)
% SOLVE_CR3BP_MINFUEL  Gate-3 eps-sharpen: certified CR3BP ENERGY (gain=1)
% -> certified CR3BP MIN-FUEL, at a fixed t_f, via the shared homotopy_mee
% eps: 1 -> 0 sweep.
%
% Loads the Task-4 artifact (bridge_mu_continuation.m's certified, full-Moon
% gain=1 ENERGY solution at tfTarget = 1.5*tfmin) and sharpens it to the
% fuel problem (eps=0) with par.pert.gain pinned at 1 throughout (spec D5:
% the gain knob is a Phase-1 STAGE-1 bridge device only; by gate 3 the walk
% is already complete and this driver never touches gain again). The
% eps-schedule itself, and all per-step resume caching, belong to
% homotopy_mee.m (E2B/core/homotopy_mee.m) -- this driver does not
% reimplement any of that machinery, only assembles its inputs and enforces
% the post-hoc certification block.
%
% CERTIFICATION (review amendment E) -- ALL of:
%   best.certified && best.epsReached==0 && strcmp(best.ipoptStatus,
%   'Solve_Succeeded') && best.maxDefect<1e-6 && best.maxUnit<1e-8 &&
%   best.termErr<1e-8, AND no casadi_lt_mee:boundSaturation warning fired
%   during any eps-step solve in this invocation (enforced live via a
%   warning->error escalation around the homotopy_mee call, not a
%   post-hoc lastwarn() check -- lastwarn() would only see the LAST of
%   ~14 sequential solves and silently miss an earlier occurrence).
%
% CAVEAT (recorded, not silently skipped): the 2-body PMP/primer verifier
% mee_primer_switch.m (E2B/verify/mee_primer_switch.m) is NOT valid under
% lunar gravity as-is -- its Hamiltonian/primer derivation assumes
% lt_mee_rhs's 2-body dXdt; under par.pert (gain=1) dXdt carries an
% additional lunar term that is NOT proportional to thrust, so the
% zero-throttle ballistic dXdt would have to be subtracted out of the
% costate/primer bracket before any switching-function argument would
% still hold. A CR3BP-aware primer check is a recorded TODO (not
% implemented here); certification below rests on the four numeric NLP
% metrics (defect/unit-norm/terminal-error/IPOPT status) plus the
% bound-saturation check, NOT on primer/switching-function agreement.
%
% INPUTS:
%   opts - (optional) struct, all fields optional:
%     .thrustN - max thrust level [N], must match an existing Task-4
%                artifact (default 10)                              [scalar]
%     .phi0    - lunar phase at t=0 [rad] (default 0)                [scalar]
%     .maxIter - IPOPT max iterations per eps-step (default 1500)    [scalar]
%
% OUTPUTS:
%   best - homotopy_mee's eps=0 solver output (casadi_lt_mee out-struct:
%          .X .U .dL .success .ipoptStatus .maxDefect .maxUnit .termErr
%          .mf .m_f_kg .dV_kms .tf .switches .edge ...), PLUS
%          .epsReached (0 iff certified) .certified
%
% Saves E3B/results/minfuel_cr3bp_T<thrustTag>N_phi<phiTag>.mat holding
% best, fp (thrustN, m0kg, ispS, tfTarget, muM, DM, nM, phi0, gain=1),
% provenance (source artifact, cert2body row, dm_f, eps-table), tbl.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md
%       sec 2 D5, sec 8 gate 3, sec 9 (honesty/STOP rules).
%   [2] docs/superpowers/plans/2026-07-22-elliptic-geo-cr3bp-phase1.md Task 5.
%   [3] E3B/bridge_mu_continuation.m (Task 4; produces this driver's input
%       artifact; check_cache_fp pattern replicated locally below).
%   [4] earth_elliptic_to_geo/core/homotopy_mee.m (shared eps-sharpen
%       engine; resume-safe per-step caching, NOT reimplemented here).
%   [5] earth_elliptic_to_geo/reproduce/table3_certified.m (2-body
%       certified m_f_kg comparator).
%   [6] earth_elliptic_to_geo/verify/mee_primer_switch.m (2-body-only
%       primer verifier; see CAVEAT above).

if nargin < 1 || isempty(opts), opts = struct(); end
setup_paths();   % adds E2B core/lib/reproduce/verify -- must run before optdef is used

d = @(f,v) optdef(opts, f, v);
thrustN = d('thrustN', 10);
phi0    = d('phi0', 0);
maxIter = d('maxIter', 1500);

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

thrTag  = num_tag(thrustN);
phiTag  = num_tag(phi0);
rungTag = sprintf('T%sN', thrTag);
tag     = sprintf('cr3bp_%s_phi%s', rungTag, phiTag);

% --- Stage 1: load the Task-4 certified full-Moon (gain=1) ENERGY artifact -
artifactFile = fullfile(resDir, sprintf('energy_cr3bp_%s_phi%s.mat', rungTag, phiTag));
assert(exist(artifactFile, 'file') == 2, 'solve_cr3bp_minfuel:missingArtifact', ...
    ['Task-4 artifact not found: %s -- run bridge_mu_continuation(struct(' ...
     '''thrustN'',%g,''phi0'',%g)) first'], artifactFile, thrustN, phi0);
S = load(artifactFile);
sigma    = S.sigma;
X        = S.X;
U        = S.U;
dL       = S.dL;
tfTarget = S.tfTarget;
fpArtifact = S.fp;

par      = kepler_lt_params(thrustN, 1500, 2000);
par.pert = lunar_params(par, phi0, 1);   % gain=1: full physical Moon (Task 5
                                          % sharpens the gate-3 certified gain=1
                                          % energy solution; gain is not
                                          % re-walked here)

fpExpected = struct('thrustN', thrustN, 'm0kg', par.m0kg, 'ispS', par.ispS, ...
    'tfTarget', tfTarget, 'muM', par.pert.muM, 'DM', par.pert.DM, ...
    'nM', par.pert.nM, 'phi0', phi0, 'gain', 1);
check_cache_fp(struct('fp', fpArtifact), fpExpected, artifactFile, tag);
assert(isfield(fpArtifact, 'gain') && fpArtifact.gain == 1, ...
    'solve_cr3bp_minfuel:notFullMoon', ['Task-4 artifact %s has gain=%.4f ' ...
    '!= 1 -- this driver only sharpens the full-Moon (gain=1) certified ' ...
    'energy solution (gate 3 precondition)'], artifactFile, fpArtifact.gain);

fprintf('SOLVE_CR3BP_MINFUEL %s: loaded %s, tfTarget=%.4f ND, gain=%.4f (full Moon)\n', ...
    tag, artifactFile, tfTarget, fpArtifact.gain);

% --- Stage 2: eps-sharpen 1 -> 0 via the shared homotopy engine -----------
homDir = fullfile(resDir, 'homotopy');
if ~exist(homDir, 'dir'), mkdir(homDir); end

ho = struct('par', par, 'x0', X(:,1), 'tfTarget', tfTarget, 'maxIter', maxIter, ...
    'resDir', homDir, 'tag', tag, 'fp', fpExpected);

% Escalate casadi_lt_mee:boundSaturation to a hard error for the DURATION of
% this call only, so a saturation fired at ANY of the ~14 sequential
% eps-step solves is caught (a post-hoc lastwarn() check would only see the
% LAST solve's warning state and silently miss an earlier one).
warnState = warning('error', 'casadi_lt_mee:boundSaturation');
cleanupWarn = onCleanup(@() warning(warnState));
try
    [best, tbl] = homotopy_mee(sigma, X, U, dL, ho);
catch ME
    if strcmp(ME.identifier, 'casadi_lt_mee:boundSaturation')
        error('solve_cr3bp_minfuel:boundSaturationFired', ['casadi_lt_mee:' ...
            'boundSaturation fired during the eps-sharpen sweep (tag=%s): ' ...
            '%s -- this is a certification-block failure per review ' ...
            'amendment E; STOP, do not treat the returned solution as ' ...
            'certified'], tag, ME.message);
    end
    rethrow(ME);
end
clear cleanupWarn;   % restore warning state immediately (belt-and-suspenders vs onCleanup)

% --- Stage 3: certification block (review amendment E, all metrics) ------
certOk = isfield(best, 'certified') && best.certified && best.epsReached == 0 && ...
    strcmp(best.ipoptStatus, 'Solve_Succeeded') && best.maxDefect < 1e-6 && ...
    best.maxUnit < 1e-8 && best.termErr < 1e-8;
if ~certOk
    error('solve_cr3bp_minfuel:certificationFailed', ['CERTIFICATION FAILED ' ...
        '(tag=%s): certified=%d epsReached=%s status=%s defect=%.3e ' ...
        'maxUnit=%.3e termErr=%.3e -- STOP, report exactly this state'], ...
        tag, best.certified, mat2str(best.epsReached), best.ipoptStatus, ...
        best.maxDefect, best.maxUnit, best.termErr);
end

fprintf(['CAVEAT: the 2-body PMP/primer verifier (mee_primer_switch) is NOT ' ...
    'valid under lunar gravity without subtracting the zero-throttle ' ...
    'ballistic dXdt (reviewer finding) -- a CR3BP-aware primer check is a ' ...
    'recorded TODO, not silently skipped.\n']);

% --- Stage 4: compare against the 2-body certified fuel mass -------------
cert2body = table3_certified(thrustN);
dmf_kg  = best.m_f_kg - cert2body.m_f_kg;
dmf_pct = 100 * dmf_kg / cert2body.m_f_kg;

fprintf('SOLVE_CR3BP_MINFUEL %s: CERTIFIED (epsReached=0, %s)\n', tag, best.ipoptStatus);
fprintf('  m_f_kg = %.4f kg  (2-body certified m_f_kg = %.4f kg, dm_f = %+.4f kg = %+.5f%%)\n', ...
    best.m_f_kg, cert2body.m_f_kg, dmf_kg, dmf_pct);
fprintf(['  switches = %d (nodal count -- mesh-band caveat, P0 protocol: not ' ...
    'independently verified against a converged switch-mesh band at this ' ...
    'CR3BP rung)\n'], best.switches);
fprintf('  edge = %.2f%%  maxDefect = %.3e  maxUnit = %.3e  termErr = %.3e\n', ...
    100 * best.edge, best.maxDefect, best.maxUnit, best.termErr);

% --- Stage 5: save the certified min-fuel artifact + provenance ----------
fp = fpExpected;
provenance = struct('sourceArtifact', artifactFile, 'tag', tag, ...
    'thrustN', thrustN, 'phi0', phi0, 'maxIter', maxIter, ...
    'cert2body', cert2body, 'dmf_kg', dmf_kg, 'dmf_pct', dmf_pct, ...
    'epsSchedTable', tbl, 'boundSaturationChecked', true, ...
    'timestamp', datetime('now'));
outFile = fullfile(resDir, sprintf('minfuel_cr3bp_%s_phi%s.mat', rungTag, phiTag));
save(outFile, 'best', 'fp', 'provenance', 'tbl');
fprintf('SOLVE_CR3BP_MINFUEL %s: saved %s\n', tag, outFile);
end

% ---------------------------------------------------------------------------
function s = num_tag(v)
% NUM_TAG  Filename-safe numeric tag: integers -> plain digits ('10'),
% non-integers -> decimal point replaced by 'p' ('0.5' -> '0p5'), negative
% sign replaced by 'm'. Mirrors bridge_mu_continuation.m's helper of the
% same name (kept local so this file has no cross-file helper dependency).
if abs(v - round(v)) < 1e-9
    s = sprintf('%d', round(v));
else
    s = strrep(sprintf('%g', v), '.', 'p');
end
s = strrep(s, '-', 'm');
end

% ---------------------------------------------------------------------------
function check_cache_fp(S, fp, file, tag)
% CHECK_CACHE_FP  Fail-loud cache-fingerprint guard, replicated locally from
% homotopy_mee.m (E2B/core/homotopy_mee.m:119 onward), as bridge_mu_
% continuation.m already does for its own tag namespace. Here it validates
% the LOADED Task-4 artifact's .fp against the config this driver is about
% to run under (thrustN/m0kg/ispS/tfTarget/muM/DM/nM/phi0/gain) -- a stale
% or mismatched artifact under the same rung/phi tag must never be silently
% sharpened. Two backward-compat cases, mirroring the source helper:
%   (1) NO .fp AT ALL -- WARN and trust as-is, no per-field comparison
%       possible.
%   (2) SCHEMA-OLDER .fp: a field present in the CURRENT fp is absent from
%       the loaded one -- schema evolution, not a configuration mismatch.
%       WARN and treat as compatible. The hard error is reserved for fields
%       present on BOTH sides with different values.
%
% INPUTS:  S [struct, loaded artifact, must expose S.fp] fp [struct,
%          current fingerprint] file [char, artifact path, for messages]
%          tag [char, cache/run tag, for messages]
% OUTPUTS: none (WARNs or errors)
%
% REFERENCES: [1] earth_elliptic_to_geo/core/homotopy_mee.m (pattern source).
%             [2] earth_elliptic_to_geo_CR3BP/direct/bridge_mu_continuation.m
%                 (sibling local replica, Task 4).
if ~isfield(S, 'fp')
    warning('solve_cr3bp_minfuel:noCachedFingerprint', ['%s has no cached ' ...
        'config fingerprint (pre-fix artifact) -- trusting it because ' ...
        'tag=''%s'' matches; regenerate the artifact to regain fingerprint ' ...
        'protection'], file, tag);
    return;
end
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f)
        warning('solve_cr3bp_minfuel:fpSchemaOlder', ['%s: field ''%s'' ' ...
            'present in current fp but absent from the artifact (schema ' ...
            'evolution) -- trusting as compatible under tag=''%s'''], ...
            file, f, tag);
        continue;
    end
    if ~isequal(S.fp.(f), fp.(f))
        error('solve_cr3bp_minfuel:fingerprintMismatch', ['config ' ...
            'fingerprint mismatch vs %s: field ''%s'' differs between the ' ...
            'artifact and the current config (tag=''%s'') -- regenerate ' ...
            'the artifact or fix opts'], file, f, tag);
    end
end
end
