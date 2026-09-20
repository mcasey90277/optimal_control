function c = arrival_sheet_as_catalog(S, ref, sD, outMat)
%% Purpose:
%
%   Wrap an ARRIVAL SHEET -- the certified winners at one departure phase and
%   many arrival phases (build_arrival_sheet) -- as a one-row phase catalog,
%   so that everything that reads a catalog (interp_study,
%   phase_catalog_interp) can read a sheet. A sheet is cheap to build at a
%   resolution a full library is not (96 arrival phases in six hours), which
%   makes it the place to study interpolation along arrival phase.
%
%  ASSUMPTIONS / NOTES:
%
% • The problem identity (orbits, engine, constants) is COPIED from a
%   reference catalog: a sheet does not carry it. The caller vouches that
%   the sheet was built for that problem.
% • Family codes: a column's family is the ARC its winner came from, the two
%   walking directions of one anchor counted as one family; winners that
%   came from a seed root share one further code. 0 where nothing certified.
%   These are this file's own codes, not family_map's.
% • It is a VIEW for study, not a deliverable: no verdict fields are carried.
%
%% Inputs:
%
%  S                        struct                  arrival sheet: .sA [1 x nA]
%                                                   .TF [1 x nA] days (NaN =
%                                                   none) .Z8 [8 x nA] .cand
%                                                   {1 x nA} .arcs {names}
%
%  ref                      struct                  a phase catalog of the same
%                                                   problem (.constants
%                                                   .thruster .rungs_N
%                                                   .sheets(1) orbit keys)
%
%  sD                       double                  the sheet's departure
%                                                   phase [0]
%
%  outMat                   char (optional)         save the catalog here (as
%                                                   variable `catalog`)
%
%% Outputs:
%
%  c                        struct                  .sheets(1) with .sD_frac
%                                                   .sA_frac .has_solution
%                                                   .tf_nd .entry_index .z8
%                                                   .family_index + the orbit
%                                                   keys; .constants .thruster
%                                                   .rungs_N .n_entries
%                                                   .families.labels .note
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(sD), sD = 0; end
nA = numel(S.sA);
has = isfinite(S.TF(:).') & all(isfinite(S.Z8), 1);

% the family of each column: the arc its winner came from
labels = {};  fam = zeros(1, nA);
for jc = find(has)
    cands = S.cand{jc};  name = 'seed';
    k = find([cands.ok] & abs([cands.tfDays] - S.TF(jc)) < 1e-9, 1);
    if ~isempty(k) && isfield(cands, 'arc') && ~isempty(cands(k).arc) && cands(k).arc >= 1
        tok = regexp(S.arcs{cands(k).arc}, 'arrival_arc_[^_]+_(.+?)_(dn|up)', 'tokens', 'once');
        if ~isempty(tok), name = tok{1}; end
    end
    [known, where] = ismember(name, labels);
    if ~known, labels{end+1} = name;  where = numel(labels); end %#ok<AGROW>
    fam(jc) = where;
end

r = ref.sheets(1);
sheet = struct('sD_frac', sD, 'sA_frac', S.sA(:).', 'has_solution', has, 'tf_nd', NaN(1, nA), ...
               'entry_index', zeros(1, nA), 'z8', S.Z8(:, has), 'family_index', fam, ...
               'tauDRO', r.tauDRO, 'Np', r.Np, 'pm', r.pm, 'period_tulip_nd', r.period_tulip_nd);
sheet.entry_index(has) = 1:nnz(has);
sheet.tf_nd(has) = S.Z8(8, has);
c = struct('sheets', sheet, 'constants', ref.constants, 'thruster', ref.thruster, 'rungs_N', ref.rungs_N, ...
           'n_entries', nnz(has), 'families', struct('labels', {labels}), ...
           'note', sprintf('an arrival sheet at s_D = %g wrapped as a one-row catalog (arrival_sheet_as_catalog); a view for study', sD));
if nargin >= 4 && ~isempty(outMat), catalog = c;  save(outMat, 'catalog'); end
end
