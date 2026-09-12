%% BUILD_70MN_LIBRARY  The whole 70 mN DRO-to-Tulip costate library, end to end
%
%   The companion of transfer_study: that script studies ONE transfer with the
%   scaffolding exposed; this one runs the CHAIN that produces the library,
%   with the same discipline -- parameter blocks first, every stage named,
%   every stage's output a file that the next stage reads. Nothing here is
%   new machinery: each stage is one existing front door.
%
%     0  parameters and stage switches   (engine, orbits, grid, budgets)
%     1  prerequisites                   (the two certified ANCHORS must exist)
%     2  arrival-phase ARCS              (pseudo-arclength continuation, 4 arcs,
%                                         HOURS each -- normally run in batch)
%     3  the arrival SHEET              (certify every grid crossing, assemble)
%     4  departure RIBS                  (walk sD off every certified column)
%     5  PACKAGE                         (sheet + ribs -> catalog .mat)
%     6  AUDIT                           (re-derive every entry from its own keys)
%     7  second-order SWEEP              (spectrum / H6 / lift margins, sidecar;
%                                         v3 since 2026-09-11 = resolved scan)
%     8  pictures                        (branch map, phase sheet, findings torus)
%     9  deliverable                     (refuses without a clean audit; ship
%                                         decision is Mike's, so OFF by default)
%
%   Stage switches exist because the arcs and ribs take hours: a stage that
%   is off is SKIPPED and its saved output reused, so the chain can be
%   resumed from any point. Every stage after 2 checks that its input file
%   exists and names the stage that makes it.
%
%   OUTPUTS AND INPUTS ARE SEPARATE. Anchors, arcs and ribs are always read
%   from results/; the sheet, catalog, audit, sidecar and pictures go to
%   outDir, which defaults to results/ (rebuild in place). Point outDir
%   elsewhere to rebuild BESIDE the shipped files and compare. A batch
%   driver can set the struct `chainOverrides` (.outDir, .run) before
%   calling this script instead of editing it.
%
%   The package stage will not overwrite a catalog that carries the
%   second-order sweep's measurements unless the sweep stage is on
%   (guard_catalog_overwrite), and backs up whatever it overwrites.
%
%   Every candidate on either axis goes through ONE gate stack (certify_root):
%   polish, admissible flight, pointwise Pontryagin checks, foreign witness
%   and its own flight, conjugate test, hypothesis gates, H6. The audit then
%   repeats the recipient's view of the shipped file. See FINDINGS 37-41 and
%   ../../doc/CERTIFICATION_DISCIPLINE.md.
%
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.

clearvars -except chainOverrides; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
cd(here)                                    % every path below is results/...
resDir = 'results';
if ~isfolder(resDir), mkdir(resDir); end

%% ========================================================================
%  0. PARAMETERS AND STAGE SWITCHES
%% ========================================================================
engine = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150);
orbits = struct('tauDRO', 1.0, 'NpTulip', 7, 'pmTulip', -1);
grid   = struct('nA', 12, 'nD', 12, 'sD0', 0, 'sA0', 0.0754);   % sA0 = the anchor's phase
tag    = '70mN';

% the two certified anchors the arcs start from (stage 1 checks them)
anchors = { ...                       % name        seed .mat                              sA0
    'anchor', 'results/mintime_70mN_anchor.mat',     grid.sA0;          ...
    'cell11', 'results/mintime_70mN_certified.mat',  grid.sA0 + 10/12};

% arc budgets: the 2026-09-09 arcs used nStep 4000 / 5 h each
arc = struct('nStep', 4000, 'deadlineSec', 5*3600, 'span', 1.02, ...
             'levels', grid.sA0 + (-grid.nA:2*grid.nA)/grid.nA);     % unwrapped ladder

% rib budgets: one rib per certified arrival column, 11 departure points
rib = struct('direction', -1, 'nPts', grid.nD - 1, 'wallSec', 900);

run = struct('arcs', false, ...       % hours; the four arcs are normally batch jobs
             'sheet', true, 'ribs', false, 'package', true, 'audit', true, ...
             'sweep', false, 'pictures', true, 'deliverable', false);

% where the OUTPUTS go (inputs always come from resDir); see the header
outDir = resDir;
if exist('chainOverrides', 'var')
    if isfield(chainOverrides, 'outDir'), outDir = chainOverrides.outDir; end
    if isfield(chainOverrides, 'run')
        for sw = fieldnames(chainOverrides.run)'
            assert(isfield(run, sw{1}), 'unknown stage switch "%s"', sw{1});
            run.(sw{1}) = chainOverrides.run.(sw{1});
        end
    end
end
if ~isfolder(outDir), mkdir(outDir); end

% the sidecar is the THIRD sweep's. v1: lift margins differed from the
% catalog's by up to 1.35e4. v2 (2026-09-10/11): the plateau classifier --
% its 42 "near-miss" cells were plateaus under nested refinement, not
% located positive minima, and endpoint clusters were never refined
% (FINDINGS 48). v3: conj_spectrum with located minima, shifted grids, the
% UNRESOLVED class and endpoint refinement. A sidecar resumes by skipping
% its done records, so pointing at v2 would re-write the OLD verdicts and
% call it a sweep; the name change is what forces the re-measurement.
files = struct( ...
    'sheet',   fullfile(outDir, sprintf('arrival_sheet_%s.mat', tag)), ...
    'catalog', fullfile(outDir, sprintf('costate_catalog_dro_tulip_%s.mat', tag)), ...
    'audit',   fullfile(outDir, sprintf('audit_%s.mat', tag)), ...
    'sidecar', fullfile(outDir, 'second_order_progress_v3.mat'), ...
    'branch',  fullfile(outDir, 'arrival_branch_map.png'), ...
    'torus',   fullfile(outDir, sprintf('phase_torus_%s.png', tag)), ...
    'findings', fullfile(outDir, 'phase_torus_findings.png'));
setupOpts = struct('thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
                   'tauDRO', orbits.tauDRO, 'NpTulip', orbits.NpTulip, 'sD', grid.sD0);
fprintf('0. LIBRARY %s: %.0f mN / Isp %g s / %g kg, DRO tau %g -> %d-petal tulip, %d x %d grid\n', ...
        tag, engine.thrustN*1000, engine.ispS, engine.m0kg, orbits.tauDRO, orbits.NpTulip, grid.nD, grid.nA);
fprintf('   inputs from %s, outputs to %s\n', resDir, outDir);

%% ========================================================================
%  1. PREREQUISITES -- the anchors. An anchor is one CERTIFIED root (z + it.Y)
%     at a known phase pair; the first one came from the direct solve of the
%     abstract case polished by ms_tfmin (probe_abstract_case), the second
%     from the first sheet's cell (1,11). A new operating point needs a new
%     anchor: run_dro_tulip walks to one.
%% ========================================================================
for k = 1:size(anchors, 1)
    assert(isfile(anchors{k,2}), ...
        'anchor "%s" missing: %s. Make one with run_dro_tulip and save z + it.Y.', ...
        anchors{k,1}, anchors{k,2});
    fprintf('1. anchor %-8s %s  (sA0 = %.4f)\n', anchors{k,1}, anchors{k,2}, anchors{k,3});
end
pool = capped_pool();                       % the fence every external call runs behind

%% ========================================================================
%  2. ARRIVAL ARCS -- pseudo-arclength continuation in the homogeneous chart,
%     one arc per anchor per direction, recording every crossing of the grid
%     levels (folds included, so a level can hold several candidates).
%% ========================================================================
arcFiles = {};
for k = 1:size(anchors, 1)
    for direction = [-1, +1]
        nm = sprintf('%s_%s', anchors{k,1}, pick(direction < 0, 'dn', 'up'));
        f  = fullfile(resDir, sprintf('arrival_arc_%s_long.mat', nm));
        arcFiles{end+1} = f; %#ok<SAGROW>
        if ~run.arcs
            fprintf('2. arc %-10s %s\n', nm, pick(isfile(f), 'reused', 'MISSING (run.arcs is off)'));
            continue
        end
        so = setupOpts;  so.sA0 = anchors{k,3};  so.anchorMat = anchors{k,2};
        [~, anc] = arclength_arrival('setup', so);
        ao = so;  ao.direction = direction;  ao.sAStop = anchors{k,3} + direction*arc.span;
        ao.levels = arc.levels;  ao.nStep = arc.nStep;  ao.deadlineSec = arc.deadlineSec;
        ao.logFile = fullfile(resDir, sprintf('arrival_arc_%s_long.log', nm));
        A = arclength_arrival(anc, ao);
        save(f, 'A', '-v7.3');
        fprintf('2. arc %-10s %d roots, sA %.4f -> %.4f, %d folds, %d crossings, stop = %s\n', ...
                nm, numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), A.stop);
    end
end
assert(any(cellfun(@isfile, arcFiles)), 'no arrival arcs in %s: turn run.arcs on', resDir);

%% ========================================================================
%  3. THE ARRIVAL SHEET -- every grid crossing of every arc through the gate
%     stack, then assembled: per arrival phase, the fastest CERTIFIED root.
%% ========================================================================
if run.sheet
    S = build_arrival_sheet(setfield(setupOpts, 'out', files.sheet)); %#ok<SFLD>
    plot_arrival_arcs(arcFiles(cellfun(@isfile, arcFiles)), S, files.branch);
else
    assert(isfile(files.sheet), 'sheet missing: %s (stage 3 makes it)', files.sheet);
    L = load(files.sheet);  S = L.S;
    fprintf('3. sheet reused: %s\n', files.sheet);
end
certifiedCols = find(~isnan(S.TF(:).'));          % a column is certified iff it has a t_f
fprintf('3. sheet: %d of %d arrival phases certified\n', numel(certifiedCols), grid.nA);

%% ========================================================================
%  4. DEPARTURE RIBS -- off every certified column, one file per column so a
%     rib that stalls costs only its own row.
%% ========================================================================
ribFiles = {};
for j = certifiedCols
    f = fullfile(resDir, sprintf('arrival_rib_col%02d.mat', j));
    if run.ribs
        build_ribs(files.sheet, struct('only', j, 'direction', rib.direction, ...
                   'nPts', rib.nPts, 'wallSec', rib.wallSec, 'out', f));
    end
    if isfile(f), ribFiles{end+1} = f; end %#ok<SAGROW>
end
% ribs from the 2026-09-09/10 campaign were grouped by letter (arrival_ribA..E);
% carry every rib file present, once
allRibs = dir(fullfile(resDir, 'arrival_rib*.mat'));
ribFiles = unique([ribFiles, fullfile(resDir, {allRibs.name})]);
fprintf('4. ribs: %d rib files\n', numel(ribFiles));

%% ========================================================================
%  5. PACKAGE -- only CERTIFIED entries are carried; the problem identity is
%     stamped from the sheet, never from defaults.
%% ========================================================================
if run.package
    bakCat = guard_catalog_overwrite(files.catalog, run.sweep);
    if ~isempty(bakCat), fprintf('5. previous catalog backed up to %s\n', bakCat); end
    cat_ = package_phase_catalog(files.sheet, ribFiles, struct('tag', tag, ...
        'thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
        'nD', grid.nD, 'sD0', grid.sD0, 'outDir', outDir));
else
    assert(isfile(files.catalog), 'catalog missing: %s (stage 5 makes it)', files.catalog);
    L = load(files.catalog);  fn = fieldnames(L);  cat_ = L.(fn{1});   % one variable, the catalog
end
sh = cat_.sheets(1);
fprintf('5. catalog %s: %d entries, %d of %d cells, conj PASS on %d\n', files.catalog, ...
        cat_.n_entries, nnz(sh.has_solution), numel(sh.has_solution), nnz(sh.conj_pass == 1));

%% ========================================================================
%  6. AUDIT -- the recipient's view: rebuild every entry from the file's own
%     keys, fly it, re-run the witness and the second-order verdicts.
%% ========================================================================
% A BAD ENTRY IS NOT A REASON TO STOP MEASURING. Until 2026-09-12 a single
% bad audit row aborted the chain here, so the second-order sweep -- the
% stage that produces the evidence a bad row needs interpreting WITH --
% never ran. The audit is diagnostic; the SHIP gate (stage 9) is where
% badness must block. Blockers are collected and reported at the end.
blockers = {};
if run.audit
    Aud = audit_phase_catalog(files.catalog, struct('out', files.audit, 'pool', pool));
    fprintf('6. audit: %d ok / %d bad\n', Aud.nOk, Aud.nBad);
    if Aud.nBad > 0
        fprintf('   %s\n', Aud.problems{:});
        blockers{end+1} = sprintf('audit: %d of %d entries bad', Aud.nBad, Aud.nOk + Aud.nBad);
        fprintf(['   the chain CONTINUES (the sweep still measures every entry); the\n' ...
                 '   deliverable stage will refuse. See the blocker summary at the end.\n']);
    end
else
    if isfile(files.audit)
        fprintf('6. audit reused: %s\n', files.audit);
    else
        fprintf('6. audit MISSING (%s) and stage 6 is off: nothing may ship\n', files.audit);
        blockers{end+1} = 'audit never run';
    end
end

%% ========================================================================
%  7. SECOND-ORDER SWEEP -- dense singular spectrum, H6 margin, lift margin
%     for every entry; sidecar after each so it resumes for free; written
%     into the catalog only after a COMPLETE census.
%% ========================================================================
% Each stage is fenced: a stage that throws records a blocker and the chain
% goes on to the next one, so one failure never costs the measurements of
% every stage behind it.
if run.sweep
    try
        Sw = second_order_pass(files.catalog, struct('sideMat', files.sidecar));
        fprintf('7. sweep census: %d done, %d to do\n', Sw.nDone, Sw.nTodo);
        if Sw.done
            Sw = second_order_pass(files.catalog, struct('sideMat', files.sidecar, 'writeback', true));
            fprintf('7. written back: worst interior crossings %g, worst H6 %.2fx, worst lift %.0fx\n', ...
                    Sw.worstSpectrum, Sw.worstH6, Sw.worstLift);
        else
            blockers{end+1} = sprintf('sweep incomplete: %d entries still to measure', Sw.nTodo);
        end
    catch ME
        fprintf('7. sweep THREW: %s\n', ME.message);
        blockers{end+1} = ['sweep threw: ' ME.message];
    end
else
    fprintf('7. sweep skipped (run.sweep is off)\n');
end

%% ========================================================================
%  8. PICTURES
%% ========================================================================
if run.pictures
    try
        set(0, 'DefaultFigureVisible', 'off');
        Q = load(fullfile(outDir, sprintf('phase_sheet_%s', tag), ...
                          sprintf('dro_tulip_%s_tau%g_Np%d.mat', tag, orbits.tauDRO, orbits.NpTulip)));
        plot_phase_sheet(Q, files.torus);
        I = plot_torus_findings(files.catalog, files.findings);
        fprintf('8. pictures: %s, %s (%d of %d certified)\n', files.torus, files.findings, I.nCert, I.nCells);
        close all
    catch ME
        fprintf('8. pictures THREW: %s\n', ME.message);
        blockers{end+1} = ['pictures threw: ' ME.message];   % never blocks the ship gate on its own
    end
end

%% ========================================================================
%  9. DELIVERABLE -- refuses without the audit; the ship decision is Mike's.
%% ========================================================================
shipBlockers = blockers(~contains(blockers, 'pictures threw'));
if run.deliverable
    if ~isempty(shipBlockers)
        fprintf('9. deliverable REFUSED -- %d blocker(s) stand:\n', numel(shipBlockers));
        fprintf('     - %s\n', shipBlockers{:});
    else
        out = build_dro_deliverable(struct('catMat', files.catalog, 'auditMat', files.audit));
        fprintf('9. deliverable: %s\n', out.zip);
    end
else
    fprintf('9. deliverable not built (run.deliverable is off; it needs a clean audit)\n');
end

%% ========================================================================
%  BLOCKER SUMMARY -- what stands between this run and a shippable library.
%  The chain having finished is NOT the same as the library being clean, so
%  the two are printed separately and the caller can read either.
%% ========================================================================
fprintf('\n');
if isempty(blockers)
    fprintf('CHAIN CLEAN: every stage that ran, passed.\n');
else
    fprintf('CHAIN COMPLETE WITH %d BLOCKER(S):\n', numel(blockers));
    fprintf('  - %s\n', blockers{:});
    fprintf('Nothing ships until these clear; the measurements above are still valid.\n');
end
chainBlockers = blockers;      % the batch driver reads this

%% ------------------------------------------------------------------------
function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end
