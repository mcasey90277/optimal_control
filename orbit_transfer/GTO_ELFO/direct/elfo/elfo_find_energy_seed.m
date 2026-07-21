function [seedFile, tfSeed, relErr] = elfo_find_energy_seed(resDir, tfTarget, relTol)
% ELFO_FIND_ENERGY_SEED  Locate the ELFO energy seed nearest a target t_f.
%
% Scans resDir for energy_elfo_*.mat seeds and selects by the seed's PHYSICAL
% final time X(8,end) — not by filename key. This makes seed lookup immune to
% factor-convention changes: legacy seeds keyed by tulip-anchored factors
% (energy_elfo_f####.mat, #### = 1000*t_f/tfMin_tulip) and new ELFO-anchored
% keys coexist safely, because the stored t_f is the authoritative coordinate
% (2026-07-21 review triage C1).
%
% INPUTS:
%   resDir   - directory holding energy_elfo_*.mat seed files [char]
%   tfTarget - requested transfer time [ND, scalar]
%   relTol   - max |tfSeed-tfTarget|/tfTarget to accept [scalar, default 0.02]
%
% OUTPUTS:
%   seedFile - path of the accepted seed ('' if none within relTol) [char]
%   tfSeed   - that seed's stored t_f [ND] (NaN if none)
%   relErr   - |tfSeed-tfTarget|/tfTarget for the returned seed (Inf if none)
%
% REFERENCES:
%   [1] GTO_ELFO/doc/reviews/2026-07-21_triage.md (C1 anchor rebase).
if nargin < 3 || isempty(relTol), relTol = 0.02; end
files = dir(fullfile(resDir, 'energy_elfo_*.mat'));
seedFile = '';  tfSeed = NaN;  relErr = Inf;
ws = warning('off', 'MATLAB:load:variableNotFound');   % non-seed matches skip quietly
cleanupW = onCleanup(@() warning(ws));
for k = 1:numel(files)
    fp = fullfile(files(k).folder, files(k).name);
    try
        T = load(fp, 'X');
    catch
        continue;                        % unreadable/foreign file: skip
    end
    if ~isfield(T, 'X') || size(T.X,1) < 8, continue; end
    tfk = T.X(8, end);
    rk  = abs(tfk - tfTarget) / tfTarget;
    if rk < relErr
        relErr = rk;  tfSeed = tfk;  seedFile = fp;
    end
end
if relErr > relTol
    seedFile = '';                       % nearest seed too far: honest miss
end
end
