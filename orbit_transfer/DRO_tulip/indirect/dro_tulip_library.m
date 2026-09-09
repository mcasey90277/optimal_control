function lib = dro_tulip_library(here)
%% Purpose:
%
%   THE certified DRO -> tulip minimum-time solutions on disk at 70 mN,
%   gathered from every result file that holds one, as a single struct
%   array. One home, so the single-transfer front door and the sheet
%   builder seed themselves from the same list instead of each keeping its
%   own idea of what is already solved.
%
%   These are solutions, not verdicts: a consumer that intends to rely on
%   one re-certifies it (certify_root). The list carries where each came
%   from so a stale entry can be traced.
%
%% Inputs:
%
%  here                     char (optional)         the indirect/ folder
%                                                   [this file's folder]
%
%% Outputs:
%
%  lib                      struct array            .sD .sA (phase
%                                                   fractions) .tfDays .src
%                                                   .file .z [8x1] .Y
%                                                   [14 x K] .K
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1 || isempty(here), here = fileparts(mfilename('fullpath')); end
tStar = 382981.289129055;
lib = struct('sD', {}, 'sA', {}, 'tfDays', {}, 'src', {}, 'file', {}, 'z', {}, 'Y', {}, 'K', {});
f = fullfile(here, 'results', 'mintime_70mN_anchor.mat');
if isfile(f)
    L = load(f);
    lib(end+1) = struct('sD', 0, 'sA', 0.0754, 'tfDays', L.z(8)*tStar/86400, ...
                        'src', 'anchor', 'file', f, 'z', L.z(:), 'Y', L.it.Y, 'K', size(L.it.Y,2));
end
f = fullfile(here, 'results', 'mintime_70mN_certified.mat');
if isfile(f)
    L = load(f);
    lib(end+1) = struct('sD', 0, 'sA', mod(0.0754 + 10/12, 1), 'tfDays', L.best.z(8)*tStar/86400, ...
                        'src', 'cell(1,11)', 'file', f, 'z', L.best.z(:), 'Y', L.best.it.Y, ...
                        'K', size(L.best.it.Y,2));
end
for nm = {'sweep_phase_mintime.mat', 'sweep_phase_70mN.mat'}
    f = fullfile(here, 'results', nm{1});
    if ~isfile(f), continue, end
    L = load(f);  S = L.S;
    [ii, jj] = find(isfinite(S.TF));
    for k = 1:numel(ii)
        sDk = S.sD(ii(k));  sAk = S.sA(jj(k));
        if any(abs([lib.sD] - sDk) < 1e-9 & abs([lib.sA] - sAk) < 1e-9), continue, end
        Yk = S.Yj{ii(k), jj(k)};
        lib(end+1) = struct('sD', sDk, 'sA', sAk, 'tfDays', S.TF(ii(k),jj(k))*tStar/86400, ...
            'src', nm{1}, 'file', f, 'z', S.Z8(:, ii(k), jj(k)), 'Y', Yk, 'K', size(Yk, 2)); %#ok<AGROW>
    end
end
end
