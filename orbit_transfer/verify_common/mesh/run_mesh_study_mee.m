function out = run_mesh_study_mee(rows, factors, opts)
% RUN_MESH_STUDY_MEE  Tier A driver: the earth order study.
%
% The study's core measurement. For each banked earth MEE row, re-solves it on
% a four-level mesh ladder, collects the per-quantity series, and writes a
% report. This is what discriminates between the two competing predictions the
% study's external reviewers made (H1 p~1 vs H2 p~2 for final mass and switch
% times).
%
% CACHE ISOLATION IS A HARD CONSTRAINT. Every artifact this writes goes to
% results/mesh_study/ under a MESH_ prefix. A study solve must not be able to
% load or clobber a production .mat -- the campaign's certified rows are the
% repo's headline numbers and this is a report-only study. (Precedent: the
% Table-3 reproducer's REPRO_ isolation.)
%
% THE LIVE MODEL IS STRIPPED BEFORE SAVING. mesh_ladder_mee attaches
% out.model, which holds live CasADi Opti handles; save() cannot serialize them
% ("Warning: Serializing SWIG objects not supported") and would write a .mat
% whose model field is silently useless. It is removed rather than saved
% broken.
%
% SAVING IS PER LEVEL, NOT PER ROW. The 1 N row at x8 is ~13900 nodes and
% projects to hours; a crash or a kill part-way through must not discard the
% levels already solved.
%
% INPUTS:
%   rows    - row tags to study [1xR cellstr]. Default
%             {'MEE_M2_10N','MEE_M2_2p5N','MEE_M2_1N'}. 2.5 N is deliberately
%             in the default set: it carried the transversality question that
%             Task 0 retracted, and it is the mid-thrust check on whether the
%             order established at 10 N survives more switches.
%   factors - node-count multipliers [1xL]. Default [1 2 4 8] -- four levels,
%             so the sliding three-level slopes can be checked against each
%             other; three would give a single fragile slope.
%   opts    - struct (optional):
%             .modRoot     earth campaign root [default: resolved from here]
%             .maxIterBase IPOPT cap at factor 1 [default: each row's own]
%             .resDir      output directory [default <modRoot>/results/mesh_study]
%             .verbose     [default true]
%
% OUTPUTS:
%   out - struct array [1xR]: .row .S (mesh_collect) .T (report table)
%         .wall .ok .why
%
% SIDE EFFECTS: writes, per row,
%   <resDir>/MESH_<row>_x<factor>.mat   one per solved level (model stripped)
%   <resDir>/mesh_MESH_<row>.mat        the collected study + table
%   <resDir>/mesh_study_log.txt         appended progress log
%
% REFERENCES:
%   [1] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Task 5.
%   [2] mesh_ladder_mee.m, mesh_collect.m, mesh_report.m.

if nargin < 1 || isempty(rows),    rows = {'MEE_M2_10N','MEE_M2_2p5N','MEE_M2_1N'}; end
if nargin < 2 || isempty(factors), factors = [1 2 4 8]; end
if nargin < 3, opts = struct(); end
if ischar(rows), rows = {rows}; end

here    = fileparts(mfilename('fullpath'));
addpath(here);
modRoot = local_default(opts, 'modRoot', ...
    fullfile(fileparts(fileparts(here)), 'earth_elliptic_to_geo', 'direct'));
resDir  = local_default(opts, 'resDir', fullfile(modRoot, 'results', 'mesh_study'));
verbose = local_default(opts, 'verbose', true);
if ~isfolder(resDir), mkdir(resDir); end
logF = fullfile(resDir, 'mesh_study_log.txt');

out = struct('row',{},'S',{},'T',{},'wall',{},'ok',{},'why',{});
local_log(logF, sprintf('\n===== Tier A start: rows {%s}, factors [%s] =====', ...
    strjoin(rows,', '), num2str(factors)));

for r = 1:numel(rows)
    tag  = rows{r};
    mat  = fullfile(modRoot, 'results', [tag '.mat']);
    tRow = tic;
    if ~isfile(mat)
        local_log(logF, sprintf('%-14s SKIP -- no such row at %s', tag, mat));
        out(end+1) = struct('row',tag,'S',[],'T',[],'wall',0,'ok',false, ...
                            'why','row .mat not found'); %#ok<AGROW>
        continue
    end

    local_log(logF, sprintf('%-14s starting ladder [%s]', tag, num2str(factors)));
    lopts = struct('verbose', verbose);
    if isfield(opts,'maxIterBase') && ~isempty(opts.maxIterBase)
        lopts.maxIterBase = opts.maxIterBase;
    end

    ok = true;  why = '';  S = [];  T = [];
    try
        L = mesh_ladder_mee(mat, factors, lopts);
        for k = 1:numel(L)
            local_log(logF, sprintf('  %-12s x%-3g N=%-6d %-7s wall %7.1f s  m_f %.6f  sw %d', ...
                tag, L(k).factor, L(k).N, local_tf(L(k).ok,'OK','FAILED'), L(k).wall, ...
                L(k).out.m_f_kg, L(k).out.switches));
            local_save_level(resDir, tag, L(k));
        end
        nOk = sum([L.ok]);
        if nOk < 2
            ok = false;  why = sprintf('only %d level(s) converged', nOk);
        else
            S = mesh_collect(L);
            T = mesh_report(S, ['MESH_' tag], resDir);
            if nOk < numel(factors)
                why = sprintf('ladder stopped after %d of %d levels', nOk, numel(factors));
                local_log(logF, sprintf('  %-12s PARTIAL: %s', tag, why));
            end
        end
    catch err
        ok = false;  why = sprintf('%s: %s', err.identifier, err.message);
        local_log(logF, sprintf('  %-12s ERROR %s', tag, why));
    end

    w = toc(tRow);
    out(end+1) = struct('row',tag,'S',S,'T',T,'wall',w,'ok',ok,'why',why); %#ok<AGROW>
    local_log(logF, sprintf('%-14s done in %.1f s (%s)', tag, w, local_tf(ok,'ok',why)));
end

local_log(logF, sprintf('===== Tier A end: %d/%d rows ok =====', sum([out.ok]), numel(out)));
end

% ---------------------------------------------------------------------------
function local_save_level(resDir, tag, Lk)
% LOCAL_SAVE_LEVEL  Persist one ladder level under the MESH_ isolation prefix.
% The live CasADi model is removed first: save() cannot serialize SWIG handles
% and would otherwise write a silently broken field.
% INPUTS: resDir, tag [char]; Lk - one ladder level    OUTPUTS: none
lvl = Lk;
if isfield(lvl.out, 'model'), lvl.out = rmfield(lvl.out, 'model'); end
f = fullfile(resDir, sprintf('MESH_%s_x%g.mat', tag, Lk.factor));
save(f, 'lvl', '-v7.3');
end

% ---------------------------------------------------------------------------
function local_log(logF, msg)
% LOCAL_LOG  Echo to stdout and append to the run log (stdout is buffered
% under matlab -batch, so the file is the reliable progress record).
% INPUTS: logF [char], msg [char]    OUTPUTS: none
fprintf('%s\n', msg);
fid = fopen(logF, 'a');
if fid > 0
    fprintf(fid, '%s  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); %#ok<TNOW1,DATST>
    fclose(fid);
end
end

% ---------------------------------------------------------------------------
function s = local_tf(c, a, b)
% LOCAL_TF  Pick a label by condition. INPUTS: c,a,b   OUTPUTS: s [char]
if c, s = a; else, s = b; end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
