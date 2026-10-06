function lib = dro_tulip_library(here, opts)
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
%  opts                     struct (optional)
%   .includeCatalog [false] also list the LIBRARY OF RECORD's primary entries
%   (library_70mN_24x48_v2) that no result file already holds. They carry
%   z8 AND the stored junction states (.Y [14 x K], .K). Off by default, because build_arrival_sheet seeds from
%   this list and seeding a sheet with its own previous output would be
%   circular.
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
%  09/11/2026  opts.includeCatalog: the front door knew 10 of the catalog's
%              115 certified entries and walked to the other 105
%  10/06/2026  anchor and cell(1,11) sourced from the library of record (the
%              2026-09-08 files missed their keys by 49 / 22 km under the
%              2026-09-12 endpoint rule); direct_certified skips duplicate
%              phase pairs like the other sources
%  10/06/2026  includeCatalog reads the library of record (with its junctions);
%              the old catalog missed one key by 2.1 km
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1 || isempty(here), here = fileparts(mfilename('fullpath')); end
tStar = 382981.289129055;
lib = struct('sD', {}, 'sA', {}, 'tfDays', {}, 'src', {}, 'file', {}, 'z', {}, 'Y', {}, 'K', {});
% The anchor (0, 0.0754) and cell(1,11) (0, 0.9087) come from the library of
% record, NOT from results/mintime_70mN_anchor.mat / mintime_70mN_certified.mat:
% those were solved on 2026-09-08 under the pre-2026-09-12 endpoint rule (before
% the periodic-spline phase_state) and sit 49 km / 22 km off their keys under
% the current rule. The files stay on disk untouched; they are never loaded here.
% src 'anchor' is kept because the certification tests select the entry by it.
rec = fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat');
if isfile(rec)
    L = load(rec);  fn = fieldnames(L);  s = L.(fn{1}).sheets(1);
    cells = {'anchor', 0, 0.0754; 'cell(1,11)', 0, mod(0.0754 + 10/12, 1)};
    for k = 1:size(cells, 1)
        ii = find(abs(s.sD_frac - cells{k,2}) < 1e-9);
        jj = find(abs(s.sA_frac - cells{k,3}) < 1e-9);
        if isempty(ii) || isempty(jj) || ~s.has_solution(ii, jj, 1)
            warning('dro_tulip_library:noRecordCell', ...
                    'library of record has no solution at (%.4f, %.4f); %s skipped', cells{k,2}, cells{k,3}, cells{k,1});
            continue
        end
        e  = s.entry_index(ii, jj, 1);
        zk = s.z8(:, e);  Yk = s.junctions{e};
        lib(end+1) = struct('sD', cells{k,2}, 'sA', cells{k,3}, 'tfDays', zk(8)*tStar/86400, ...
                            'src', cells{k,1}, 'file', rec, 'z', zk(:), 'Y', Yk, 'K', size(Yk, 2));
    end
else
    warning('dro_tulip_library:noRecord', ...
            'library of record missing (%s); anchor and cell(1,11) skipped', rec);
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
% solutions found by DIRECT solves and certified (FINDINGS 61): the family
% the continuation never reached. z8 only; the consumer rebuilds the seed
% from it (seed_from_z8), which puts the polish at a root.
f = fullfile(here, 'results', 'mintime_70mN_direct_certified.mat');
if isfile(f)
    L = load(f);
    for k = 1:numel(L.direct)
        e = L.direct(k);
        if any(abs([lib.sD] - e.sD) < 1e-9 & abs([lib.sA] - e.sA) < 1e-9), continue, end
        lib(end+1) = struct('sD', e.sD, 'sA', e.sA, 'tfDays', e.tfDays, 'src', 'direct_certified', ...
                            'file', f, 'z', e.z(:), 'Y', [], 'K', []); %#ok<AGROW>
    end
end
% the library-of-record CATALOG, on request (see opts above). Primaries only;
% the stored junction states ride along. The older results/costate_catalog_
% dro_tulip_70mN.mat predates the 2026-09-12 endpoint rule and is never read.
if nargin >= 2 && isstruct(opts) && isfield(opts, 'includeCatalog') && opts.includeCatalog
    if isfile(rec)
        L = load(rec);  fn = fieldnames(L);  s = L.(fn{1}).sheets(1);
        [ii, jj] = find(s.has_solution(:,:,1));
        for k = 1:numel(ii)
            sDk = s.sD_frac(ii(k));  sAk = s.sA_frac(jj(k));
            if any(abs([lib.sD] - sDk) < 1e-9 & abs([lib.sA] - sAk) < 1e-9), continue, end
            e  = s.entry_index(ii(k), jj(k), 1);
            zk = s.z8(:, e);  Yk = s.junctions{e};
            lib(end+1) = struct('sD', sDk, 'sA', sAk, 'tfDays', zk(8)*tStar/86400, ...
                                'src', 'catalog', 'file', rec, 'z', zk, 'Y', Yk, 'K', size(Yk, 2));
        end
    else
        warning('dro_tulip_library:noRecord', ...
                'library of record missing (%s); catalog entries skipped', rec);
    end
end
end
