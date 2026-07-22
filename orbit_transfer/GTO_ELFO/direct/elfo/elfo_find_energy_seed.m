function [seedFile, tfSeed, relErr] = elfo_find_energy_seed(resDir, tfTarget, relTol, fpNow)
% ELFO_FIND_ENERGY_SEED  Locate the ELFO energy seed nearest a target t_f.
%
% Scans resDir for energy_elfo_*.mat seeds and selects by the seed's PHYSICAL
% final time X(8,end) — not by filename key. This makes seed lookup immune to
% factor-convention changes: legacy seeds keyed by tulip-anchored factors
% (energy_elfo_f####.mat, #### = 1000*t_f/tfMin_tulip) and new ELFO-anchored
% keys coexist safely, because the stored t_f is the authoritative coordinate
% (2026-07-21 review triage C1).
%
% Optional thrust-rung filter (2026-07-21 ladder-prep, spec sec 4): when fpNow
% is given, a candidate seed's stored fp.Tmax is checked against fpNow.Tmax --
% a mismatch skips the seed. A LEGACY seed (no stored .fp at all) is a POLICY
% call, not a schema gap: it is treated as NOMINAL (25 mN), so it is eligible
% only when fpNow.thrustN==0.025 (one aggregate warning per scan) and skipped
% under any off-nominal fpNow. Omitting fpNow entirely keeps today's behavior
% exactly (no filter at all -- every candidate is eligible by t_f alone).
%
% INPUTS:
%   resDir   - directory holding energy_elfo_*.mat seed files [char]
%   tfTarget - requested transfer time [ND, scalar]
%   relTol   - max |tfSeed-tfTarget|/tfTarget to accept [scalar, default 0.02]
%   fpNow    - (optional) current fingerprint (cr3bp_fingerprint) [struct];
%              omit for the legacy no-filter behavior
%
% OUTPUTS:
%   seedFile - path of the accepted seed ('' if none within relTol, or none
%              eligible under fpNow) [char]
%   tfSeed   - that seed's stored t_f [ND] (NaN if none)
%   relErr   - |tfSeed-tfTarget|/tfTarget for the returned seed (Inf if none)
%
% REFERENCES:
%   [1] GTO_ELFO/doc/reviews/2026-07-21_triage.md (C1 anchor rebase).
%   [2] docs/superpowers/specs/2026-07-21-ladder-prep-design.md sec 4 (seed
%       fp filter policy: legacy seeds are nominal, not schema-agnostic).
if nargin < 3 || isempty(relTol), relTol = 0.02; end
if nargin < 4, fpNow = []; end
files = dir(fullfile(resDir, 'energy_elfo_*.mat'));
seedFile = '';  tfSeed = NaN;  relErr = Inf;
ws = warning('off', 'MATLAB:load:variableNotFound');   % non-seed matches skip quietly
cleanupW = onCleanup(@() warning(ws));
legacyWarned = false;
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    try
        T = load(fpath, 'X', 'fp');
    catch
        continue;                        % unreadable/foreign file: skip
    end
    if ~isfield(T, 'X') || size(T.X,1) < 8, continue; end
    if ~isempty(fpNow)
        if isfield(T, 'fp') && isfield(T.fp, 'Tmax')
            if abs(T.fp.Tmax - fpNow.Tmax) > 1e-12
                continue;                 % off-rung fingerprint: skip
            end
        else
            % legacy (fingerprint-less) seed: POLICY = treat as NOMINAL (25 mN)
            if abs(fpNow.thrustN - 0.025) > 1e-12
                continue;                 % off-nominal request: skip
            end
            if ~legacyWarned
                warning('elfo_find_energy_seed:legacySeed', ...
                    '%s: legacy (fingerprint-less) seed(s) in %s treated as nominal (25 mN)', ...
                    mfilename, resDir);
                legacyWarned = true;
            end
        end
    end
    tfk = T.X(8, end);
    rk  = abs(tfk - tfTarget) / tfTarget;
    if rk < relErr
        relErr = rk;  tfSeed = tfk;  seedFile = fpath;
    end
end
if relErr > relTol
    seedFile = '';                       % nearest seed too far: honest miss
end
end
