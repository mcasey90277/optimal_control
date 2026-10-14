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
% • build.auditVerdict / build.note are what the CALLER passes
%   (opts.auditVerdict / opts.auditNote): this function sees no audit, so
%   it stamps 'not audited at export' unless told otherwise. The builder
%   gates first and passes 'CLEAN' with a note built from the audit counts.
% • Source rewriting keeps the sheet kind: a source naming
%   sheet96_resolution_test reads '96-phase sheet candidate, column N
%   (...)', the 48-phase sheet's 'sheet candidate, column N (...)'.
% • export_deliverable_catalog('cleanSource', s) returns the rewritten
%   source string s (char or cellstr) without exporting (the test seam).
%
%% Inputs:
%
%  recordMat                char                    path of the record .mat
%
%  outMat                   char                    path of the export .mat
%
%  opts                     struct   (optional)     .varName
%                                                   ['costate_catalog_dro_tulip_70mN']
%                                                   .auditVerdict (char)
%                                                   ['not audited at export']
%                                                   .auditNote (char) ['']
%                                                   .deliverable (char)
%                                                   ['deliverable 8'] the
%                                                   label in the provenance
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
%  M. Casey  absolute paths in top-level text (families.arcDir, any
%            char/cellstr outside sheets/alternatives) cut to their last
%            component; conj_test counts recomputed from the per-cell
%            verdict grid                                   (c) 10/13/2026
%  M. Casey  opts.auditVerdict / opts.auditNote (no CLEAN stamp without
%            an audit); 96-phase sheet sources keep their kind; the
%            'cleanSource' mode                             (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if ischar(recordMat) && strcmp(recordMat, 'cleanSource')   % test seam
    if iscell(outMat), c = cellfun(@cleanSource, outMat, 'UniformOutput', false);
    else, c = cleanSource(outMat); end
    info = struct();
    return
end
if nargin < 3, opts = struct(); end
varName = fieldd(opts, 'varName', 'costate_catalog_dro_tulip_70mN');
verdict = fieldd(opts, 'auditVerdict', 'not audited at export');
aNote   = fieldd(opts, 'auditNote', '');
delivLabel = fieldd(opts, 'deliverable', 'deliverable 8');
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
    'auditVerdict', verdict, ...
    'note', aNote);
c.provenance = sprintf(['%s Status layer (FINDINGS 98): candidates harvested from the builds, every inferred one ' ...
    're-certified, primaries matched to their junctions by root, assembled beside the unchanged primaries, ' ...
    'audited fail-closed (verdict at export: %s), borderline statuses relabelled to the meet of two runs. Exported for %s on %s.'], ...
    c.provenance, verdict, delivLabel, char(datetime('now', 'Format', 'yyyy-MM-dd')));
n = 0;
for k = 1:numel(c.alternatives)
    s0 = c.alternatives(k).source;  s1 = cleanSource(s0);
    if ~strcmp(s0, s1), n = n + 1; end
    c.alternatives(k).source = s1;
end
top = setdiff(fieldnames(c), {'sheets', 'alternatives'}, 'stable');
for k = 1:numel(top), c.(top{k}) = stripPaths(c.(top{k})); end
if isfield(c, 'conj_test') && isfield(c.sheets, 'conj_pass')
    vc = cell(numel(c.sheets), 1);
    for k = 1:numel(c.sheets)          % the verdict of every solved cell
        h = logical(c.sheets(k).has_solution);  g = c.sheets(k).conj_pass;
        vc{k} = g(h);
    end
    v = vertcat(vc{:});
    cnt = struct('nPass', nnz(v == 1), 'nFail', nnz(v == 0), 'nNotrun', nnz(v == -1), ...
                 'nUndecided', nnz(~ismember(v, [1 0 -1])));
    fc = fieldnames(cnt);
    for k = 1:numel(fc)
        if k <= 2 || isfield(c.conj_test, fc{k}), c.conj_test.(fc{k}) = cnt.(fc{k}); end
    end
    c.conj_test.note = 'counts recomputed from the per-cell verdicts at export';
end
assert(strcmp(catalog_content_key(c), primKey) && strcmp(alternatives_content_key(c), altKey), 'export_deliverable_catalog:keys', 'the export changed a content key');
p = catalog_schema('validate', c);
assert(isempty(p), 'export_deliverable_catalog:schema', 'the export does not validate: %s', strjoin(p, '; '));
S.(varName) = c;  save(outMat, '-struct', 'S');
info = struct('primKey', primKey, 'altKey', altKey, 'removed', {removed}, 'nSourceRewritten', n);
end

function x = stripPaths(x)
% STRIPPATHS  Recursively replace every absolute path inside char / string
% values of a struct or cell with its bare last path component.
% INPUTS:  x (any).  OUTPUTS: x (same class and size).
if ischar(x)
    x = regexprep(x, '(?<![\w.])/(?:[\w.\-]+/)+([\w.\-]+)/?', '$1');
elseif isstring(x)
    for k = 1:numel(x), x(k) = string(stripPaths(char(x(k)))); end
elseif iscell(x)
    for k = 1:numel(x), x{k} = stripPaths(x{k}); end
elseif isstruct(x)
    fn = fieldnames(x);
    for e = 1:numel(x)
        for q = 1:numel(fn), x(e).(fn{q}) = stripPaths(x(e).(fn{q})); end
    end
end
end

function s = cleanSource(s)
% CLEANSOURCE  Drop file paths and file names from a source string, keep the
% kind words and the sheet kind (96-phase sheet vs the 48-phase grid).
% INPUTS: s [char].  OUTPUTS: s [char], never empty.
is96 = contains(s, 'sheet96_resolution_test');
s = regexprep(s, '[\w./-]*\.mat', '');
s = regexprep(s, '(reproduce_|sheet96_resolution_test|departure_rib96_test|round_\d+)\S*', '');
s = regexprep(s, '^\s*cand(?=\W|$)', 'sheet candidate');
s = strtrim(regexprep(s, '\s+', ' '));
s = regexprep(s, '^\s*,\s*', '');
if is96, s = regexprep(s, '^sheet candidate', '96-phase sheet candidate'); end
if isempty(s), s = 'build candidate'; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field of s, or the default when absent or empty.
% INPUTS:  s [struct]; f [char]; d_ default.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
