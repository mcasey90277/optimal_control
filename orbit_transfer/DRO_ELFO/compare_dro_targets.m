% compare_dro_targets
%   The matched comparison of two arrival targets for the same departure orbit: DRO -> tulip (the control)
%   and DRO -> ELFO at one or more plane orientations, each a torus sheet of the same grid and the same thrust
%   ladder written by run_dro_elfo_torus. The script loads the sheets from out_dir, puts a difficulty scorecard
%   (transfer_difficulty, revolutions about the Moon) on each, prints one combined table per rung, the ladder
%   summary of every target and the per-rung ELFO-minus-tulip differences, and saves the scorecards and one
%   figure. With solve = false a missing sheet is an error naming the command that makes it; with solve = true
%   the missing sheets are run first (short runs only: a 6 x 6 torus is hours). A file <tag>_conj.mat beside a
%   sheet, holding conjPass [nD x nA x rungs], adds the refuted share; deep-probe sheets in out_dir/deep_probe/
%   (<tag>.mat) add a second table below 0.5 N. A struct compare_opts in the workspace overrides knobs field by
%   field for ONE run (test seam): the script clears it when it has read it.

%% knobs
dro_tau        = 1;                   % DRO period [ND], 1 = 4.43 d
n_phase        = 6;                   % torus grid, n_phase x n_phase cells per sheet
rungs          = [15 12 10 7 5 3 2 1.5 1 0.75 0.5];   % thrust ladder [N], high to low
elfo_raan_deg  = [0 90];              % ELFO plane orientations, one sheet each [deg, integer]
deep_rungs     = [0.3 0.2 0.15 0.1 0.07];             % thrust of the deep probe [N] (the probe itself is run outside)
n_deep_cells   = 4;                   % full-ladder cells per sheet in the deep probe
out_dir        = '';                  % '' = results/compare beside this script
solve          = false;               % true = run run_dro_elfo_torus for the missing sheets first
make_plot      = true;                % draw the figure and save the PNG next to the .mat
fig_visible    = 'on';                % figure visibility, 'off' for unattended runs
torus_fcn      = @run_dro_elfo_torus; % test seam, leave as is

knob_names = {'dro_tau', 'n_phase', 'rungs', 'elfo_raan_deg', 'deep_rungs', 'n_deep_cells', 'out_dir', ...
              'solve', 'make_plot', 'fig_visible', 'torus_fcn'};
if exist('compare_opts', 'var')
    ov = compare_opts;
    clear compare_opts                % cleared BEFORE validating: a refused field must not poison the next run
    for ff = fieldnames(ov).'
        if ~ismember(ff{1}, knob_names), error('compare_dro_targets:knob', 'unknown compare_opts field ''%s''', ff{1}); end
        eval([ff{1} ' = ov.(ff{1});']);
    end
    clear ov ff
end
here = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = fullfile(here, 'results', 'compare'); end
setup_dro_elfo_paths();
t_star_s = 382981.289129055;          % time unit [s]

%% the sheets
% target list: the tulip first (the control), then one ELFO sheet per orientation
tags   = {sprintf('dro_tulip_tau%g', dro_tau)};
labels = {'tulip'};
for raan = elfo_raan_deg(:).'
    tags{end + 1}   = sprintf('dro_elfo_tau%g_raan%03d', dro_tau, mod(raan, 360));
    labels{end + 1} = sprintf('elfo raan %g', raan);
end
nT = numel(tags);
sheet_files = fullfile(out_dir, strcat(tags, '.mat'));
missing = find(~cellfun(@isfile, sheet_files));
if ~isempty(missing)
    sheets = struct('target', {}, 'droTau', {}, 'raanDeg', {});
    for km = missing
        if km == 1
            sheets(end + 1) = struct('target', 'tulip', 'droTau', dro_tau, 'raanDeg', NaN);
        else
            sheets(end + 1) = struct('target', 'elfo', 'droTau', dro_tau, 'raanDeg', elfo_raan_deg(km - 1));
        end
    end
    torus_opts = struct('sheets', sheets(:), 'nD', n_phase, 'nA', n_phase, 'rungs', rungs, 'outDir', out_dir);
    if ~solve
        error('compare_dro_targets:missingSheet', ...
              ['sheet file %s not found; create it with solve = true or with run_dro_elfo_torus([], ' ...
               'struct(''sheets'', <sheets>, ''nD'', %d, ''nA'', %d, ''rungs'', %s, ''outDir'', ''%s'')) ' ...
               'with <sheets> = struct(''target'', %s, ''droTau'', %g, ''raanDeg'', %s)'], ...
              sheet_files{missing(1)}, n_phase, n_phase, mat2str(rungs), out_dir, ...
              ['''' sheets(1).target ''''], dro_tau, mat2str(sheets(1).raanDeg));
    end
    if ~isfolder(out_dir), mkdir(out_dir); end
    fprintf('solving %d missing sheet(s): %s\n', numel(missing), strjoin(tags(missing), ' '));
    torus_fcn([], torus_opts);
    still = find(~cellfun(@isfile, sheet_files));
    if ~isempty(still)
        error('compare_dro_targets:missingSheet', 'sheet file %s still missing after the solve', sheet_files{still(1)});
    end
end
attempted = zeros(1, nT);  total = zeros(1, nT);
for kt = 1:nT                         % grid, ladder and identity must be the knobs'
    S = load(sheet_files{kt}, 'rungs', 'sD', 'sA', 'ATT', 'sheet_key');
    if numel(S.sD) ~= n_phase || numel(S.sA) ~= n_phase || numel(S.rungs) ~= numel(rungs) ...
            || any(abs(S.rungs(:) - rungs(:)) > 1e-12)
        error('compare_dro_targets:sheetMismatch', ...
              'sheet %s is %d x %d with rungs %s; the knobs ask %d x %d with rungs %s', sheet_files{kt}, ...
              numel(S.sD), numel(S.sA), mat2str(S.rungs(:).'), n_phase, n_phase, mat2str(rungs(:).'));
    end
    if ~isfield(S, 'sheet_key')
        error('compare_dro_targets:sheetMismatch', 'sheet %s has no sheet_key variable: its identity cannot be verified', ...
              sheet_files{kt});
    end
    if kt == 1, want = struct('target', 'tulip', 'depSize', dro_tau, 'raanDeg', NaN);
    else,       want = struct('target', 'elfo',  'depSize', dro_tau, 'raanDeg', mod(elfo_raan_deg(kt - 1), 360)); end
    got = S.sheet_key;
    if ~(strcmp(char(got.target), want.target) && got.depSize == want.depSize ...
            && (isequaln(got.raanDeg, want.raanDeg)))
        error('compare_dro_targets:sheetMismatch', ...
              'sheet %s is loaded as %s but its sheet_key says target %s, depSize %g, raanDeg %g; expected target %s, depSize %g, raanDeg %g', ...
              sheet_files{kt}, tags{kt}, char(got.target), got.depSize, got.raanDeg, want.target, want.depSize, want.raanDeg);
    end
    attempted(kt) = nnz(S.ATT > 0);  total(kt) = numel(S.ATT);
end
complete = attempted == total;
fprintf('sheets (%d x %d cells, rungs %s N)\n', n_phase, n_phase, mat2str(rungs));
for kt = 1:nT
    fprintf('%s\n', attempt_line(labels{kt}, attempted(kt), total(kt)));
end
notice = '';
if ~all(complete)
    notice = ['INCOMPLETE: ' strjoin(arrayfun(@(k) sprintf('%s (%d of %d cells attempted)', labels{k}, attempted(k), total(k)), ...
              find(~complete), 'UniformOutput', false), '; ')];
    warning('compare_dro_targets:incomplete', '%s', notice);
end

%% scorecards
% D{k}.perRung  table per rung (thrustN, shareAccepted, tfMedianND, revMedian, dVmedianKms, roughness,
%               shareRefuted when a conjPass file is present), D{k}.ladder, D{k}.meta, D{k}.lines
D = cell(1, nT);
for kt = 1:nT
    o = struct('primary', 'moon', 'print', false);
    conj_file = fullfile(out_dir, [tags{kt} '_conj.mat']);
    if isfile(conj_file)
        Q = load(conj_file, 'conjPass');
        if ~isequal(size(Q.conjPass), [n_phase n_phase numel(rungs)])
            error('compare_dro_targets:sheetMismatch', 'conjPass in %s has size %s, the sheet is %d x %d x %d', ...
                  conj_file, mat2str(size(Q.conjPass)), n_phase, n_phase, numel(rungs));
        end
        o.conjPass = Q.conjPass;
    end
    D{kt} = transfer_difficulty(sheet_files{kt}, o);
end

%% combined table
% combined_table: one row per rung; per target the columns <label>_share, _tf, _rev, _dV, _rough, _refuted
nR = numel(rungs);
combined_table = table(rungs(:), 'VariableNames', {'thrustN'});
for kt = 1:nT
    combined_table = [combined_table, target_columns(labels{kt}, D{kt}.perRung)];
end
% every sheet's rungs equal the knob rungs (enforced by sheetMismatch), so row k is the same thrust in all scorecards
if ~isempty(notice), fprintf('\n%s\n', notice); end
fprintf('%s\n', strjoin(table_lines(labels, D, true(nR, 1)), newline));

%% ladders
fprintf('\nladder summary\n');
fprintf('%-14s %10s %10s %6s %6s %16s %26s\n', 'target', 'lowest90', 'lowestAny', 'full', 'empty', 'cells attempted', ...
        'median attempts per attempted cell');
for kt = 1:nT
    L = D{kt}.ladder;
    fprintf('%-14s %10s %10s %6d %6d %16s %26s\n', labels{kt}, num_or_dash(L.lowestRungAllCells, '%.4g'), ...
            num_or_dash(L.lowestRungAnyCell, '%.4g'), L.nFullLadders, L.nEmptyCells, ...
            sprintf('%d of %d', L.nAttempted, total(kt)), num_or_dash(L.attMedian, '%.4g'));
end

%% differences
% diffs: one row per ELFO sheet and rung: dShare = share(elfo) - share(tulip), tfRatio = median t_f(elfo) / (tulip)
diffs = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, ...
              'VariableNames', {'target', 'thrustN', 'dShare', 'tfRatio'});
fprintf('\ndifferences (elfo minus tulip)\n');
if ~isempty(notice), fprintf('%s\n', notice); end
for kt = 2:nT
    P0 = D{1}.perRung;  P1 = D{kt}.perRung;
    dsh = P1.shareAccepted - P0.shareAccepted;
    rat = P1.tfMedianND ./ P0.tfMedianND;
    diffs = [diffs; table(repmat(string(labels{kt}), nR, 1), rungs(:), dsh, rat, 'VariableNames', diffs.Properties.VariableNames)];
    fprintf('%s\n%8s %10s %10s\n', labels{kt}, 'thrust N', 'd share', 'tf ratio');
    for kr = 1:nR
        fprintf('%8g %10s %10s\n', rungs(kr), num_or_dash(dsh(kr), '%+.3f'), num_or_dash(rat(kr), '%.3f'));
    end
end

%% deep probe
% deep(k): tag, label, thrustN (rungs below 0.5 N), perRung (the scorecard rows of those rungs)
deep = struct('tag', {}, 'label', {}, 'thrustN', {}, 'perRung', {});
deep_D = {};
for kt = 1:nT
    deep_file = fullfile(out_dir, 'deep_probe', [tags{kt} '.mat']);
    if ~isfile(deep_file), continue, end
    Dd = transfer_difficulty(deep_file, struct('primary', 'moon', 'print', false));
    keep = Dd.perRung.thrustN < 0.5;
    deep(end + 1) = struct('tag', tags{kt}, 'label', labels{kt}, 'thrustN', Dd.perRung.thrustN(keep), ...
                           'perRung', Dd.perRung(keep, :));
    deep_D{end + 1} = Dd;
end
if isempty(deep)
    fprintf('\ndeep probe: not run\n');
else
    fprintf('\nbelow 0.5 N\n');
    for kd = 1:numel(deep)
        mask = deep_D{kd}.perRung.thrustN < 0.5;
        fprintf('%s\n%s\n', deep(kd).label, strjoin(table_lines(deep(kd).label, deep_D(kd), mask), newline));
    end
end

%% save
if ~isfolder(out_dir), mkdir(out_dir); end
cards = struct('tag', tags, 'label', labels, 'D', D, 'attempted', num2cell(attempted), 'total', num2cell(total), ...
               'complete', num2cell(complete));
knobs = struct('dro_tau', dro_tau, 'n_phase', n_phase, 'rungs', rungs, 'elfo_raan_deg', elfo_raan_deg, ...
               'deep_rungs', deep_rungs, 'n_deep_cells', n_deep_cells);
mat_file = fullfile(out_dir, 'compare_dro_targets.mat');
save(mat_file, 'cards', 'combined_table', 'diffs', 'deep', 'knobs', 'attempted', 'total', 'complete');
fprintf('\nsaved: %s\n', mat_file);

%% plot
if make_plot
    png_file = fullfile(out_dir, 'compare_dro_targets.png');
    fig = figure('Visible', fig_visible, 'Color', 'w', 'Position', [100 100 1000 720]);
    series = {'shareAccepted', 'share accepted', 1; 'tfMedianND', 'median t_f [d]', t_star_s/86400; ...
              'revMedian', 'median revolutions', 1; 'roughness', 'roughness', 1};
    marks  = {'o-', 's-', '^-', 'd-', 'v-'};
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
    for ks = 1:4
        ax = nexttile(tl);  hold(ax, 'on');
        for kt = 1:nT
            plot(ax, rungs, D{kt}.perRung.(series{ks, 1})*series{ks, 3}, marks{mod(kt - 1, 5) + 1}, ...
                 'LineWidth', 1.5, 'DisplayName', labels{kt});
        end
        set(ax, 'XScale', 'log', 'XDir', 'reverse');  grid(ax, 'on');
        xlabel(ax, 'thrust [N]');  ylabel(ax, series{ks, 2});
        if ks == 1, legend(ax, 'Location', 'best'); end
    end
    title(tl, sprintf('DRO (period %g ND) to tulip and ELFO, %d x %d torus', dro_tau, n_phase, n_phase));
    exportgraphics(fig, png_file, 'Resolution', 150);
    if strcmp(fig_visible, 'off'), close(fig); end
    fprintf('figure: %s\n', png_file);
end

%% local functions
function Tc = target_columns(label, P)
% TARGET_COLUMNS  Table of one target's columns, prefixed by its label.
% INPUTS: label char; P table (transfer_difficulty perRung). OUTPUTS: Tc table [nRung x 6].
nr = height(P);
if ismember('shareRefuted', P.Properties.VariableNames), ref = P.shareRefuted; else, ref = nan(nr, 1); end
pre = regexprep(label, '\W+', '_');
Tc = table(P.shareAccepted, P.tfMedianND, P.revMedian, P.dVmedianKms, P.roughness, ref, 'VariableNames', ...
           strcat([pre '_'], {'share', 'tf', 'rev', 'dV', 'rough', 'refuted'}));
end

function lines = table_lines(labels, Dc, mask)
% TABLE_LINES  Side-by-side per-rung table of several scorecards.
% INPUTS: labels char | cell of char (one per scorecard); Dc cell of transfer_difficulty structs;
%   mask logical [nRung x 1] rows to print (rungs of the first scorecard). OUTPUTS: lines cell of char.
labels = cellstr(labels);
colw = [5 6 4 5 5 5];
names = {'acc', 'tf', 'rev', 'dV', 'rough', 'refut'};
l1 = sprintf('%5s', '');  l2 = sprintf('%5s', 'thr N');
for kt = 1:numel(Dc)
    l1 = [l1 ' ' sprintf('%-*s', sum(colw) + numel(colw) - 1, labels{kt}) ' |'];
    for kc = 1:numel(colw), l2 = [l2 ' ' sprintf('%*s', colw(kc), names{kc})]; end
    l2 = [l2 ' |'];
end
lines = {l1, l2};
rows = find(mask(:)).';
for kr = rows
    ln = sprintf('%5g', Dc{1}.perRung.thrustN(kr));
    for kt = 1:numel(Dc)
        P = Dc{kt}.perRung;
        if ismember('shareRefuted', P.Properties.VariableNames), ref = P.shareRefuted(kr); else, ref = NaN; end
        vals = [P.shareAccepted(kr) P.tfMedianND(kr) P.revMedian(kr) P.dVmedianKms(kr) P.roughness(kr) ref];
        fm = {'%.3f', '%.4f', '%.2f', '%.2f', '%.2f', '%.3f'};
        for kc = 1:6
            if kc == 6 && isnan(vals(kc)), s = 'n/a'; else, s = num_or_dash(vals(kc), fm{kc}); end
            ln = [ln ' ' sprintf('%*s', colw(kc), s)];
        end
        ln = [ln ' |'];
    end
    lines{end + 1} = ln;
end
lines{end + 1} = 'acc share of cells accepted, tf median flight time [ND], rev median revolutions about the Moon,';
lines{end + 1} = 'dV median delta-V [km/s], rough t_f roughness over the torus, refut refuted share (n/a without a conj file)';
end

function s = attempt_line(label, n, nTot)
% ATTEMPT_LINE  One line per sheet: cells attempted, INCOMPLETE when not all.
% INPUTS: label char; n, nTot double. OUTPUTS: s char.
s = sprintf('%s: %d of %d cells attempted', label, n, nTot);
if n < nTot, s = [s ' INCOMPLETE']; end
end

function s = num_or_dash(x, fm)
% NUM_OR_DASH  Number as text, '-' for NaN.  INPUTS: x double scalar; fm char format. OUTPUTS: s char.
if isnan(x), s = '-'; else, s = sprintf(fm, x); end
end
