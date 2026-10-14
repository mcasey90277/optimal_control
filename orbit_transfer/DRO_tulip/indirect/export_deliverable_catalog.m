function [c, info] = export_deliverable_catalog(recordMat, outMat, opts)
%% Purpose:
%
%   Export the library of record as the cleaned DELIVERABLE catalog: the
%   internal bookkeeping (top-level merge and status_layer) is removed, a
%   build struct and a provenance paragraph are added, and every
%   alternative's source string loses its file paths. Both content keys
%   (catalog_content_key, alternatives_content_key) are asserted unchanged,
%   so no datum is touched.
%
%  ASSUMPTIONS / NOTES:
%
% • alternatives_content_key leaves .source out of its hash, which is why
%   the rewrite keeps the key; the assert proves it on every export.
% • The catalog is saved under its variable name (opts.varName), so a bare
%   load gives the struct back.
%
%% Inputs:
%
%  recordMat                char                    path of the record .mat
%
%  outMat                   char                    path of the export .mat
%
%  opts                     struct   (optional)     .varName
%                                                   ['costate_catalog_dro_tulip_70mN']
%
%% Outputs:
%
%  c                        struct                  the exported catalog
%
%  info                     struct                  .primKey .altKey (char)
%                                                   .removed (cellstr of
%                                                   top-level fields removed)
%                                                   .nSourceRewritten
%                                                   (scalar)
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if nargin < 3, opts = struct(); end
varName = fieldd(opts, 'varName', 'costate_catalog_dro_tulip_70mN');
R = load(recordMat);  c = R.(char(fieldnames(R)));
primKey = catalog_content_key(c);  altKey = alternatives_content_key(c);
removed = {};
sl = struct();  if isfield(c, 'status_layer'), sl = c.status_layer; end
for f = {'merge', 'status_layer'}
    if isfield(c, f{1}), c = rmfield(c, f{1});  removed{end+1} = f{1}; end
end
c.build = struct( ...
    'dates', '24 x 48 build 2026-09-20/21; refill 2026-10-02/03; merge with the 24 x 24 record 2026-10-03; status-layer backfill, re-certification and audit 2026-10-04/05', ...
    'nCandidates', fieldd(sl, 'nCandidates', NaN), 'nRecert', fieldd(sl, 'nRecert', NaN), 'nMoved', fieldd(sl, 'nMoved', NaN), ...
    'nBorderline', nnz(arrayfun(@(a) isfield(a, 'borderline') && isequal(a.borderline, true), c.alternatives)), ...
    'borderlineNote', fieldd(sl, 'borderlineNote', ''), ...
    'auditVerdict', 'CLEAN', ...
    'note', 'Every alternative was re-certified under the shipped certifier and audited fail-closed (1,808 ok / 0 bad); see README.');
c.provenance = sprintf(['%s Status layer (FINDINGS 98): candidates harvested from the builds, every inferred one ' ...
    're-certified, primaries matched to their junctions by root, assembled beside the unchanged primaries, ' ...
    'audited fail-closed (CLEAN), borderline statuses relabelled to the meet of two runs. Exported for deliverable 8 on %s.'], ...
    c.provenance, char(datetime('now', 'Format', 'yyyy-MM-dd')));
n = 0;
for k = 1:numel(c.alternatives)
    s0 = c.alternatives(k).source;  s1 = cleanSource(s0);
    if ~strcmp(s0, s1), n = n + 1; end
    c.alternatives(k).source = s1;
end
assert(strcmp(catalog_content_key(c), primKey) && strcmp(alternatives_content_key(c), altKey), 'export_deliverable_catalog:keys', 'the export changed a content key');
p = catalog_schema('validate', c);
assert(isempty(p), 'export_deliverable_catalog:schema', 'the export does not validate: %s', strjoin(p, '; '));
S.(varName) = c;  save(outMat, '-struct', 'S');
info = struct('primKey', primKey, 'altKey', altKey, 'removed', {removed}, 'nSourceRewritten', n);
end

function s = cleanSource(s)
% CLEANSOURCE  Drop file paths and file names from a source string, keep the
% kind words.  INPUTS: s [char].  OUTPUTS: s [char], never empty.
s = regexprep(s, '[\w./-]*\.mat', '');
s = regexprep(s, '(reproduce_|sheet96_resolution_test|departure_rib96_test|round_\d+)\S*', '');
s = regexprep(s, '^\s*cand(?=\W|$)', 'sheet candidate');
s = strtrim(regexprep(s, '\s+', ' '));
s = regexprep(s, '^\s*,\s*', '');
if isempty(s), s = 'build candidate'; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field of s, or the default when absent or empty.
% INPUTS:  s [struct]; f [char]; d_ default.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
