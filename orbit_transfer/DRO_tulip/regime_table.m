function [T, F] = regime_table(cfg)
%% Purpose:
%
%   Aggregate every finished REGIME MAP arm into one feature table: load the
%   per-job race outputs in direct/results/regime, extract features with
%   regime_features, print the "who arrived" matrix per (cell, gamma), and
%   list the ONE-FAMILY-ONLY cases -- the goal's headline output.
%
%   A one-family-only case is only a claim about the FAMILY when the losers
%   were given the same budget; this table therefore reports the budget's
%   identity (all jobs come from run_regime_map's single budget) and flags
%   any case whose losers have not yet been retried at double budget.
%
%% Inputs:
%
%  cfg                      struct (optional)       .dir [direct/results/regime]
%                                                   .print [true]
%                                                   .outMat ['' = no save]
%
%% Outputs:
%
%  T                        table                   one row per arm
%  F                        struct array            the raw feature structs
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, cfg = struct(); end
here   = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
d      = getf(cfg, 'dir', fullfile(resDir, 'regime'));
doPrint = getf(cfg, 'print', true);
partial = getf(cfg, 'includePartial', false);

Lp = load(fullfile(resDir, 'minenergy_pilot.mat'));
Lh = load(fullfile(resDir, 'minenergy_highgamma.mat'));
files = dir(fullfile(d, 'regime_*.mat'));
F = struct([]);
for k = 1:numel(files)
    L = load(fullfile(d, files(k).name));
    if ~isfield(L, 'out') || ~isfield(L.out, 'arms'), continue, end
    fam = fieldnames(L.out.arms);
    tok = regexp(files(k).name, 'regime_(\w+?)_c(\d)(\d)_g(\d+)_', 'tokens', 'once');
    src = tok{1};
    if strcmp(src, 'grid'), R = Lp.R; else, R = Lh.R; end
    cellIdx = [str2double(tok{2}) str2double(tok{3})];
    gTarget = str2double(tok{4})/100;
    ki = find(arrayfun(@(r) isequal([r.iD r.iA], cellIdx) && ...
                            abs(r.gam - gTarget) < 5e-3, R), 1);
    if isempty(ki), warning('no seed record for %s', files(k).name); continue, end
    A = L.out.arms.(fam{1});
    % A race saves after EVERY rung, so a file exists while the arm is still
    % walking; only a FINISHED arm carries wallTotal.  A mid-walk arm would
    % read as "walled" and corrupt the who-arrived matrix.
    if ~isfield(A, 'wallTotal') && ~partial, continue, end
    f = regime_features(A, R(ki), fam{1});
    f.src = src;
    if isempty(F), F = f; else, F(end+1) = f; end %#ok<AGROW>
end
if isempty(F), T = table(); return, end

T = struct2table(rmfield(F, {}), 'AsArray', true);

if doPrint
    fprintf('\n=== REGIME MAP: %d arms, one identical budget ===\n', numel(F));
    fprintf('%-5s %-7s %-7s %6s  %-7s %8s %9s %8s %5s %5s %7s %6s %6s %7s\n', ...
        'cell', 'gamma', 'family', 'tf[d]', 'outcome', 'p_floor', 'ramp', 'm_f', ...
        'fail', 'bis', 'wall[m]', 'nX', 'revs', 'min|Qd|');
    [~, ord] = sortrows([arrayfun(@(f) f.cellIdx(1)*10+f.cellIdx(2), F)', [F.gamma]']);
    for k = ord(:)'
        f = F(k);
        fprintf('(%d,%d) %-7.4g %-7s %6.1f  %-7s %8.4g %9.4g %8.6f %5d %5d %7.1f %6d %6.2f %7.3f\n', ...
            f.cellIdx(1), f.cellIdx(2), f.gamma, f.family, f.tfDays, f.outcome, ...
            f.pFloor, f.rampWidth, f.mf, f.nFail, f.nBisect, f.wallMin, ...
            f.nCross, f.revs, f.minAbsDQdt);
    end
    printVerdicts(F);
end
if ~isempty(getf(cfg, 'outMat', '')), save(cfg.outMat, 'T', 'F'); end
end

% ------------------------------------------------------------------------
function printVerdicts(F)
% PRINTVERDICTS  The who-solved-it matrix on the PHYSICAL criterion (mass
% agreement), with the procedural p-floor shown beside it.
% INPUTS: F (feature struct array).  OUTPUTS: none (prints).
V = regime_verdicts(F);
fprintf('\n--- who SOLVED it (dm_f <= 1e-3 vs the best arm of the case) ---\n');
fprintf('%-14s %-22s %-26s %10s %s\n', 'case', 'verdict', 'failed (dm_f)', 'nX solved', 'nX failed');
for k = 1:numel(V)
    v = V(k);
    fl = '';
    for m = 1:numel(v.failed)
        d = v.dMf(strcmp(v.families, v.failed{m}));
        fl = [fl sprintf('%s %.1e  ', v.failed{m}, d)]; %#ok<AGROW>
    end
    fprintf('(%d,%d)@%-8.4g %-22s %-26s %10d %s\n', v.cellIdx(1), v.cellIdx(2), v.gamma, ...
        v.verdict, fl, v.nSwitchSolved, mat2str(v.nSwitchFailed));
end
n1 = sum(arrayfun(@(v) numel(v.solved) == 1 && v.nArms == 3 && v.winnerAtFloor, V));
nf = sum(arrayfun(@(v) numel(v.solved) == 1 && v.nArms == 3 && ~v.winnerAtFloor, V));
nsc = sum(arrayfun(@(v) v.structureChange, V));
fprintf(['\n%d of %d complete cases are ONE-FAMILY-ONLY (winner AT the bang-bang limit); ' ...
         '%d more have a single best arm that never reached the limit (one family merely got ' ...
         'furthest). %d cases show a switch-structure change at the failure.\n'], ...
        n1, sum([V.nArms] == 3), nf, nsc);
end

function printOnlyOne(F)
% PRINTONLYONE  List (cell, gamma) cases where exactly one family arrived.
% INPUTS: F (feature struct array).  OUTPUTS: none (prints).
key = arrayfun(@(f) sprintf('%d%d_%.4f', f.cellIdx(1), f.cellIdx(2), f.gamma), ...
               F, 'UniformOutput', false);
[u, ~, ic] = unique(key);
fprintf('\n--- who arrived (needs all three arms present) ---\n');
for k = 1:numel(u)
    g = F(ic == k);
    if numel(g) < 3, fprintf('%-14s INCOMPLETE (%d/3 arms)\n', u{k}, numel(g)); continue, end
    arrived = {g(strcmp({g.outcome}, 'floor')).family};
    walled  = {g(~strcmp({g.outcome}, 'floor')).family};
    tag = 'all three';
    if numel(arrived) == 1, tag = sprintf('*** %s ONLY ***', upper(arrived{1}));
    elseif isempty(arrived), tag = 'NONE arrived';
    elseif numel(arrived) == 2, tag = sprintf('%s + %s', arrived{1}, arrived{2});
    end
    fprintf('%-14s %-28s  walled: %s\n', u{k}, tag, strjoin(walled, ' '));
end
end

function v = getf(s, f, d)
% GETF  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
