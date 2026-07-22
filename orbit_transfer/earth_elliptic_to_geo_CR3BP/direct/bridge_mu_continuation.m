function out = bridge_mu_continuation(opts)
% BRIDGE_MU_CONTINUATION  D5 stage-1 mu_M-continuation bridge: certified
% 2-body ENERGY (eps=1) solution -> CR3BP ENERGY solution at a given thrust
% rung, walking the lunar mass-continuation knob par.pert.gain from 0 to 1
% in warm-started, checkpointed steps.
%
% Pipeline (spec sec 6 Phase-1 step 3; plan Task 4):
%   1. Certified lookup (table3_certified) sets tfTarget = 1.5*tfmin (same
%      physical t_f convention as the 2-body ladder, spec D4).
%   2. Seed: two-pass mee_seed protocol, mirrored VERBATIM from
%      run_transfer_mee.m lines 132-161 (cheap N=50 revs probe with its
%      nRev window assert, then a full-density N=round(nodesPerRev*nRev)
%      sample), using opts.thrustN's OWN recipe values from table3_recipes
%      (package review A10: per-rung seedThr/nodesPerRev, not the 10 N
%      rung's values hardcoded regardless of thrustN; betaMode='tangential'
%      and initElems=[] stay fixed, same as the front door).
%   3. GATE 1 (pert absent): plain 2-body eps=1 energy solve at tfTarget --
%      confirms the seed/solver stack is healthy before any lunar term is
%      switched on (spec sec 8 gate 1 half; the byte-identical nominal RHS
%      path is gate 1's OTHER half, already proven in
%      test_lt_mee_rhs_pert.m).
%   4. GATE 2 (gain=0, pert PRESENT): re-solve tight-warm from gate 1 with
%      par.pert attached but gain=0 -- must reproduce gate 1's X to solver
%      tolerance (spec sec 8 gate 2: "the hook itself introduces no
%      drift").
%   5. Gain walk: par.pert.gain stepped through opts.gainSched (default
%      [0.25 0.5 0.75 1.0]), each step a tight warm solve from the previous
%      accepted iterate, gated on the same four-metric criterion. A failed
%      step bisects toward the last accepted gain and retries (never
%      propagates a non-certified iterate forward); step floor 0.05, below
%      which the walk errors 'bridge:stuck' rather than silently stalling.
%      Every ACCEPTED step is checkpointed (resume-safe: a second call
%      picks up from the last accepted gain; a failed step is never
%      cached, so resume can never skip past a failure).
%   6. Saves the final artifact and prints the one-line BRIDGE summary.
%
% All intermediate + final caches live under a fresh 'cr3bp_bridge_*'/
% 'energy_cr3bp_*' tag namespace (never 'MEE_M2_*', the 2-body campaign's
% own tags) so a legacy 2-body cache can never be silently reused here, and
% vice versa; every cache carries a config fingerprint (.fp) checked
% field-by-field on load (check_cache_fp, local subfunction below, mirrors
% homotopy_mee.m's helper of the same name).
%
% INPUTS:
%   opts - (optional) struct, all fields optional:
%     .thrustN   - max thrust level [N], must be a table3_certified rung
%                  (default 10)                                   [scalar]
%     .phi0      - lunar phase at t=0 [rad] (spec D6; default 0)  [scalar]
%     .gainSched - mu_M continuation waypoints in (0,1], ascending
%                  (default [0.25 0.5 0.75 1.0])                  [1xK]
%     .maxIter   - IPOPT max iterations per solve (default 1500)  [scalar]
%     .resume    - true -> reuse any existing per-stage cache whose
%                  fingerprint matches current config (default true) [logical]
%
% OUTPUTS:
%   out - the final gain-walk solver output (casadi_lt_mee out-struct: .X
%         .U .dL .success .ipoptStatus .maxDefect .maxUnit .termErr .mf
%         .m_f_kg .dV_kms .tf .switches .edge ...), PLUS:
%     .gainReached - last certified gain value reached (1.0 on full success) [scalar]
%     .gate1       - struct .ok .out  (gate-1 baseline solve, pert absent)
%     .gate2       - struct .ok .out .maxDrift  (gate-2 solve, gain=0)
%     .history     - struct array, one accepted gain-walk step per row:
%                     .gain .X .U .dL .out
%     .artifactFile - full path of the saved E3B/results/energy_cr3bp_*.mat
%
% Saves E3B/results/energy_cr3bp_T<thrustTag>N_phi<phiTag>.mat holding
% sigma, X, U, dL, tfTarget, fp (fp: thrustN, m0kg, ispS, tfTarget, muM,
% DM, nM, phi0, gain=gainReached).
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md
%       sec 2 D5 (mu-continuation bridge strategy), sec 8 (gates 1-3).
%   [2] docs/superpowers/plans/2026-07-22-elliptic-geo-cr3bp-phase1.md Task 4.
%   [3] Bonnard, Caillau, Picot, "Geometric and numerical techniques in
%       three-body low thrust transfers" (mu-continuation provenance).
%   [4] earth_elliptic_to_geo/direct/drivers/run_transfer_mee.m lines
%       132-161 (two-pass seed protocol, mirrored verbatim).
%   [5] earth_elliptic_to_geo/reproduce/table3_certified.m (certified
%       per-rung tfmin/m_f_kg/switches/revs).
%   [6] earth_elliptic_to_geo/core/homotopy_mee.m (check_cache_fp pattern,
%       replicated locally below).

if nargin < 1 || isempty(opts), opts = struct(); end
setup_paths();   % adds lib/optdef.m (E2B) to the path -- must run before optdef is used

d = @(f,v) optdef(opts, f, v);
thrustN   = d('thrustN', 10);
phi0      = d('phi0', 0);
gainSched = d('gainSched', [0.25 0.5 0.75 1.0]);
maxIter   = d('maxIter', 1500);
resumeOn  = d('resume', true);

% A10: gainSched validation -- must be finite, strictly ascending, in (0,1],
% and end at 1 (the walk's final accepted step must reach the full physical
% Moon; a schedule that overshoots/omits 1 or is non-monotonic is a config bug,
% not something the bisection-on-failure logic below is meant to paper over).
gainSched = gainSched(:).';
assert(~isempty(gainSched) && all(isfinite(gainSched)) && all(gainSched > 0) && ...
    all(gainSched <= 1) && all(diff(gainSched) > 0) && gainSched(end) == 1, ...
    'bridge:badGainSched', ['opts.gainSched must be finite, strictly ascending, ' ...
    'every element in (0,1], and end at 1 (got %s)'], mat2str(gainSched));

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

par      = kepler_lt_params(thrustN, 1500, 2000);
cert     = table3_certified(thrustN);
tfTarget = 1.5 * cert.tfmin;
xf       = [1;0;0;0;0];

thrTag = num_tag(thrustN);
phiTag = num_tag(phi0);
tag    = sprintf('cr3bp_bridge_T%sN_phi%s', thrTag, phiTag);

% --- Stage 1: seed (mirrors run_transfer_mee.m lines 132-161 VERBATIM) -----
% A10: seed knobs sourced from the certified per-rung recipe registry
% (table3_recipes), same as the front door (run_cr3bp_geo.m) -- this used to
% hardcode the 10 N rung's own values regardless of opts.thrustN, silently
% mis-seeding every OTHER rung. Off-table thrusts fall back to the campaign
% defaults (0.4 / 25), identical fallback to the front door's.
seedThr     = 0.4;  nodesPerRev = 25;
try
    rec = table3_recipes(thrustN);
    seedThr = rec.fuel.seedThr;  nodesPerRev = rec.fuel.npr;
catch ME
    if ~strcmp(ME.identifier, 'table3_recipes:unknownThrust'), rethrow(ME); end
end
betaMode    = 'tangential';
initElems   = [];

fpSeed = struct('thrustN', thrustN, 'm0kg', par.m0kg, 'ispS', par.ispS, ...
    'seedThr', seedThr, 'betaMode', betaMode, 'nodesPerRev', nodesPerRev, ...
    'xf', xf, 'initElems_isset', ~isempty(initElems));

probeFile = fullfile(resDir, [tag '_seed_probe.mat']);
if resumeOn && exist(probeFile, 'file')
    S = load(probeFile);  infoP = S.infoP;
    check_cache_fp(S, fpSeed, probeFile, tag);
else
    optsP = struct('thr', seedThr, 'betaMode', betaMode, 'N', 50, 'stopP', xf(1), ...
        'initElems', initElems);
    [~, ~, ~, ~, infoP] = mee_seed(par, optsP);
    fp = fpSeed;
    save(probeFile, 'infoP', 'fp');
end
assert(infoP.nRev >= 6.5 && infoP.nRev <= 9, 'bridge_mu_continuation:revsOutOfRange', ...
    'seedThr=%.3f gives nRev=%.3f, outside the required [6.5,9] window -- adjust seedThr', ...
    seedThr, infoP.nRev);
N = round(nodesPerRev * infoP.nRev);
fpSeedN = fpSeed;  fpSeedN.N = N;

seedFile = fullfile(resDir, [tag '_seed.mat']);
if resumeOn && exist(seedFile, 'file')
    S = load(seedFile);
    sigma = S.sigma;  X0 = S.X0;  U0 = S.U0;  dL0 = S.dL0;  seedInfo = S.seedInfo;
    check_cache_fp(S, fpSeedN, seedFile, tag);
else
    optsS = struct('thr', seedThr, 'betaMode', betaMode, 'N', N, 'stopP', xf(1), ...
        'initElems', initElems);
    [sigma, X0, U0, dL0, seedInfo] = mee_seed(par, optsS);
    fp = fpSeedN;
    save(seedFile, 'sigma', 'X0', 'U0', 'dL0', 'seedInfo', 'fp');
end
x0 = X0(:,1);
fprintf('BRIDGE %s: T=%g N, tfTarget=%.4f ND, seed N=%d nodes, seed revs=%.4f\n', ...
    tag, thrustN, tfTarget, N, seedInfo.nRev);

% --- Stage 2: GATE 1 -- pert absent, plain 2-body eps=1 energy solve ------
fpGate1 = fpSeedN;  fpGate1.tfTarget = tfTarget;  fpGate1.maxIter = maxIter;

gate1File = fullfile(resDir, [tag '_gate1.mat']);
if resumeOn && exist(gate1File, 'file')
    S = load(gate1File);  out1 = S.out1;
    check_cache_fp(S, fpGate1, gate1File, tag);
else
    out1 = casadi_lt_mee(sigma, X0, U0, dL0, struct('par', par, 'mode', 'fixedtf', ...
        'eps', 1, 'tfTarget', tfTarget, 'x0', x0, 'maxIter', maxIter, 'warmTight', false));
    fp = fpGate1;
    save(gate1File, 'out1', 'fp');
end
gate1ok = strcmp(out1.ipoptStatus, 'Solve_Succeeded') && out1.maxDefect < 1e-6 && ...
    out1.maxUnit < 1e-8 && out1.termErr < 1e-8;
fprintf('GATE 1 (pert absent):     status=%-22s defect=%.3e maxUnit=%.3e termErr=%.3e -> %s\n', ...
    out1.ipoptStatus, out1.maxDefect, out1.maxUnit, out1.termErr, pass_fail(gate1ok));
if ~gate1ok
    error('bridge_mu_continuation:gate1Failed', ['GATE 1 (pert absent) FAILED: status=%s ' ...
        'defect=%.3e maxUnit=%.3e termErr=%.3e'], out1.ipoptStatus, out1.maxDefect, ...
        out1.maxUnit, out1.termErr);
end

% --- Stage 3: GATE 2 -- gain=0, pert PRESENT, tight warm from gate 1 ------
par2 = par;
par2.pert = lunar_params(par, phi0, 0);

fpGate2 = fpGate1;
fpGate2.phi0 = phi0;  fpGate2.gain = 0;
fpGate2.muM  = par2.pert.muM;  fpGate2.DM = par2.pert.DM;  fpGate2.nM = par2.pert.nM;

gate2File = fullfile(resDir, [tag '_gate2.mat']);
if resumeOn && exist(gate2File, 'file')
    S = load(gate2File);  out2 = S.out2;
    check_cache_fp(S, fpGate2, gate2File, tag);
else
    out2 = casadi_lt_mee(sigma, out1.X, out1.U, out1.dL, struct('par', par2, ...
        'mode', 'fixedtf', 'eps', 1, 'tfTarget', tfTarget, 'x0', x0, 'maxIter', maxIter, ...
        'warmTight', true));
    fp = fpGate2;
    save(gate2File, 'out2', 'fp');
end
maxDriftG2 = max(abs(out2.X(:) - out1.X(:)));
gate2ok = strcmp(out2.ipoptStatus, 'Solve_Succeeded') && out2.maxDefect < 1e-6 && ...
    out2.maxUnit < 1e-8 && out2.termErr < 1e-8 && maxDriftG2 < 1e-8;
fprintf('GATE 2 (gain=0, pert on): status=%-22s defect=%.3e maxUnit=%.3e termErr=%.3e drift=%.3e -> %s\n', ...
    out2.ipoptStatus, out2.maxDefect, out2.maxUnit, out2.termErr, maxDriftG2, pass_fail(gate2ok));
if ~gate2ok
    error('bridge_mu_continuation:gate2Failed', ['GATE 2 (gain=0, pert present) FAILED: ' ...
        'status=%s defect=%.3e maxUnit=%.3e termErr=%.3e drift-vs-gate1=%.3e'], ...
        out2.ipoptStatus, out2.maxDefect, out2.maxUnit, out2.termErr, maxDriftG2);
end

% --- Stage 4: gain walk 0 -> 1, adaptive bisection-on-failure, resume-safe --
% Static identity fields (everything except .gain, which is the resumed
% walk's own progress marker, not a configuration field): a genuine
% thrustN/tfTarget/phi0/... drift under the same tag is still caught.
fpGWStatic = struct('thrustN', thrustN, 'm0kg', par.m0kg, 'ispS', par.ispS, ...
    'tfTarget', tfTarget, 'muM', par2.pert.muM, 'DM', par2.pert.DM, ...
    'nM', par2.pert.nM, 'phi0', phi0, 'gainSched', gainSched, 'maxIter', maxIter);

gainwalkFile = fullfile(resDir, [tag '_gainwalk.mat']);
history  = struct('gain', {}, 'X', {}, 'U', {}, 'dL', {}, 'out', {});
lastGood = 0;  Xw = out2.X;  Uw = out2.U;  dLw = out2.dL;
if resumeOn && exist(gainwalkFile, 'file')
    S = load(gainwalkFile);
    Scheck = S;
    if isfield(Scheck.fp, 'gain'), Scheck.fp = rmfield(Scheck.fp, 'gain'); end
    check_cache_fp(Scheck, fpGWStatic, gainwalkFile, tag);
    history = S.history;
    if ~isempty(history)
        lastGood = history(end).gain;
        Xw = history(end).X;  Uw = history(end).U;  dLw = history(end).dL;
        fprintf('BRIDGE: resuming gain walk from cached gain=%.4f\n', lastGood);
    end
end

queue     = gainSched(:).';
queue     = queue(queue > lastGood + 1e-12);
qi        = 1;
stepFloor = 0.05;
while qi <= numel(queue)
    g = queue(qi);
    par2.pert.gain = g;
    o = casadi_lt_mee(sigma, Xw, Uw, dLw, struct('par', par2, 'mode', 'fixedtf', ...
        'eps', 1, 'tfTarget', tfTarget, 'x0', x0, 'maxIter', maxIter, 'warmTight', true));
    ok = strcmp(o.ipoptStatus, 'Solve_Succeeded') && o.maxDefect < 1e-6 && ...
        o.maxUnit < 1e-8 && o.termErr < 1e-8;
    fprintf('  [gain walk] gain=%.4f (from %.4f) ok=%d status=%-22s defect=%.3e mf=%.2f kg\n', ...
        g, lastGood, ok, o.ipoptStatus, o.maxDefect, o.m_f_kg);
    if ok
        Xw = o.X;  Uw = o.U;  dLw = o.dL;  lastGood = g;
        history(end+1) = struct('gain', g, 'X', Xw, 'U', Uw, 'dL', dLw, 'out', o); %#ok<AGROW>
        fp = fpGWStatic;  fp.gain = g;
        save(gainwalkFile, 'history', 'fp');
        qi = qi + 1;
    else
        step = g - lastGood;
        mid  = lastGood + step/2;
        if step/2 < stepFloor
            error('bridge:stuck', ['Gain walk stuck: cannot advance past gain=%.6f from ' ...
                'last-good gain=%.6f (bisected step %.4f below the %.2f floor)'], ...
                g, lastGood, step/2, stepFloor);
        end
        queue = [queue(1:qi-1), mid, queue(qi:end)];
        fprintf('  [gain walk] gain=%.4f FAILED -- inserting midpoint gain=%.4f\n', g, mid);
    end
end

% --- Stage 5: save the final artifact + report ----------------------------
finalOut = history(end).out;
X  = finalOut.X;
U  = finalOut.U;
dL = finalOut.dL;
fp = struct('thrustN', thrustN, 'm0kg', par.m0kg, 'ispS', par.ispS, ...
    'tfTarget', tfTarget, 'muM', par2.pert.muM, 'DM', par2.pert.DM, ...
    'nM', par2.pert.nM, 'phi0', phi0, 'gain', lastGood);
artifactFile = fullfile(resDir, sprintf('energy_cr3bp_T%sN_phi%s.mat', thrTag, phiTag));
save(artifactFile, 'sigma', 'X', 'U', 'dL', 'tfTarget', 'fp');

fprintf('BRIDGE: gain=%g reached, defect=%.3e, m_f(energy)=%.4f kg\n', ...
    lastGood, finalOut.maxDefect, finalOut.m_f_kg);

out              = finalOut;
out.gainReached  = lastGood;
out.gate1        = struct('ok', gate1ok, 'out', out1);
out.gate2        = struct('ok', gate2ok, 'out', out2, 'maxDrift', maxDriftG2);
out.history      = history;
out.artifactFile = artifactFile;
end

% ---------------------------------------------------------------------------
function s = num_tag(v)
% NUM_TAG  Filename-safe numeric tag: integers -> plain digits ('10'),
% non-integers -> decimal point replaced by 'p' ('0.5' -> '0p5'), negative
% sign replaced by 'm'. Mirrors run_transfer_mee.m's mee_fuel_tag.m
% convention (integer/non-integer split), extended with a sign rule since
% phi0 may be negative.
if abs(v - round(v)) < 1e-9
    s = sprintf('%d', round(v));
else
    s = strrep(sprintf('%g', v), '.', 'p');
end
s = strrep(s, '-', 'm');
end

% ---------------------------------------------------------------------------
function s = pass_fail(ok)
% PASS_FAIL  'PASS'/'FAIL' label for a boolean gate verdict (print helper).
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

% ---------------------------------------------------------------------------
function check_cache_fp(S, fp, file, tag)
% CHECK_CACHE_FP  Fail-loud cache-fingerprint guard, replicated locally from
% homotopy_mee.m (E2B/core/homotopy_mee.m:119 onward) for this campaign's
% fresh 'cr3bp_*'/'energy_cr3bp_*' tag namespace. If loaded cache struct S
% carries a .fp field, compare it field-by-field against the current config
% fingerprint fp and error, naming the first mismatched field and the
% offending file, on any disagreement -- a stale cache built under a
% different thrustN/phi0/tfTarget/... must never be silently reused just
% because it happens to share tag (the only cache key). Two backward-compat
% cases, mirroring the E2B helper:
%   (1) NO .fp AT ALL -- WARN and trust as-is, no per-field comparison
%       possible.
%   (2) SCHEMA-OLDER .fp: the cache HAS a .fp, but a field present in the
%       CURRENT fp is simply ABSENT from the cached one -- schema
%       evolution, not a configuration mismatch. WARN and treat as
%       compatible. The hard error is reserved for fields present on BOTH
%       sides with different values.
%
% INPUTS:  S [struct, loaded cache] fp [struct, current fingerprint]
%          file [char, cache path, for messages] tag [char, cache tag]
% OUTPUTS: none (WARNs or errors)
%
% REFERENCES: [1] earth_elliptic_to_geo/core/homotopy_mee.m (pattern source).
if ~isfield(S, 'fp')
    warning('bridge_mu_continuation:noCachedFingerprint', ['%s has no cached ' ...
        'config fingerprint (pre-fix cache) -- trusting it because ' ...
        'tag=''%s'' matches; use a new tag to regain fingerprint protection ' ...
        'for this run'], file, tag);
    return;
end
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f)
        warning('bridge_mu_continuation:fpSchemaOlder', ['%s: field ''%s'' ' ...
            'present in current fp but absent from cache (schema evolution) ' ...
            '-- trusting as compatible under tag=''%s'''], file, f, tag);
        continue;
    end
    if ~isequal(S.fp.(f), fp.(f))
        error('bridge_mu_continuation:fingerprintMismatch', ['cached config ' ...
            'fingerprint mismatch in %s: field ''%s'' differs between the ' ...
            'cache and the current config -- stale cache from a different ' ...
            'configuration under the same tag=''%s''; delete the file or ' ...
            'use a new tag'], file, f, tag);
    end
end
end
