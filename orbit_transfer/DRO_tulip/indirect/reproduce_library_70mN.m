function out = reproduce_library_70mN(opts)
%% Purpose:
%
%   REBUILD THE 70 mN DRO -> TULIP COSTATE LIBRARY FROM ITS ANCHORS, IN ONE
%   CALL, AND CHECK THE RESULT AGAINST THE LIBRARY OF RECORD.
%
%   The library of record (results/library_70mN_24x24_final, 576 certified
%   minimum-time transfers over a 24 x 24 grid of departure and arrival
%   phases) was built by hand over eight rounds, 2026-09-13 to 09-15
%   (FINDINGS 59-72, process/PHASE_TORUS_RUNBOOK.md section 10). This script
%   is that build as one reproducible chain:
%
%       reproduce_library_70mN                       % PLAN: print, touch nothing
%       reproduce_library_70mN(struct('go', true))   % BUILD, compare, report
%
%   What it does, in the order of the numbered sections below:
%
%     0. the switches                 what to build, where, and how
%     1. the problem and the grid     one engine, one orbit pair, 24 x 24 phases
%     2. the five families            the certified roots every arc starts from
%     3. the arcs and seed roots      adopt the ten arcs already walked (or
%                                     re-walk them) and the direct-found roots
%     4. the build                    run_phase_torus: sheet, ribs, package,
%                                     audit, sweep, holes -- to a fixed point
%     5. the comparison               new catalog against the record, cell by cell
%     6. the verdict                  one line, and everything behind it in `out`
%
%  ASSUMPTIONS / NOTES:
%
% • WHY FIVE ANCHORS AND NO DISCOVERY. The hand build FOUND its families one
%   round at a time (the gap protocol). They are known now, so all five are
%   given as anchors and discovery is off: the round loop then has nothing to
%   find, and the run is the deterministic part of the build -- sheet, ribs,
%   filler, packaging. Set .discover = true to let it search as well.
% • WHY THE ARCS ARE ADOPTED. An arc is 2 to 6 hours of continuation and the
%   ten of this library are on disk, walked from exactly these anchors.
%   Adopting them makes the rebuild about 6 to 10 hours; .adoptArcs = false
%   re-walks them, about 30 hours, and reproduces the arcs too.
% • WHAT "THE SAME LIBRARY" MEANS. Every entry is re-polished, so the rebuilt
%   catalog agrees with the record to solver tolerance, not bit for bit; and
%   a column whose two fastest roots tie within certification noise may take
%   the other one. compare_phase_catalogs reports every such cell by name.
% • RESUMABLE. Everything is on disk under .outDir; calling again with the
%   same options continues where it stopped (run_phase_torus's state file).
%
%% Inputs:
%
%  opts                     struct (optional)
%   .go                     logical                 false = PLAN ONLY: print what
%                                                   would run, create nothing [false]
%   .outDir                 char                    the campaign folder
%                                                   [results/reproduce_70mN_24x24]
%   .nD, .nA                int                     the RESOLUTION: departure and
%                                                   arrival phases [24, 24]. The
%                                                   arcs are re-scanned at any nA;
%                                                   a multiple of 24 keeps every
%                                                   record cell on the grid, and
%                                                   the comparison is then made
%                                                   on those shared cells. Cost
%                                                   grows with nD*nA: 24 x 24 is a
%                                                   day, 48 x 48 about four
%   .adoptArcs              logical                 adopt the walked arcs [true]
%   .discover               logical                 run the discovery step [false]
%   .nWorkers               int                     rib workers [4]
%   .maxRounds              int                     round budget [3]
%   .compareOnly            logical                 skip the build; compare the
%                                                   catalog already in .outDir [false]
%   .print                  logical                 print sections 0-3 [true]
%
%% Outputs:
%
%  out                      struct                  .state ('planned' | 'done' |
%                                                   'budget' | 'compared') .spec
%                                                   (the run_phase_torus spec)
%                                                   .nArcsFound .nRegistryRoots
%                                                   .torus (the
%                                                   driver's output) .catalog
%                                                   .comparison .audit .ok
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

% TEST SEAM: handles to this file's local functions (tests/test_reproduce_library)
if nargin == 1 && ischar(opts) && strcmp(opts, 'localfunctions'), out = localHandles(localfunctions);  return, end
if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
res  = fullfile(here, 'results');
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));

%% ========================================================================
%  0. THE SWITCHES
%% ========================================================================
go          = pick(opts, 'go', false);
nD          = pick(opts, 'nD', 24);              % departure phases  } the RESOLUTION. 24 x 24 is the record's;
nA          = pick(opts, 'nA', 24);              % arrival phases    } any multiple of 24 keeps every record cell on the grid
outDir      = pick(opts, 'outDir', fullfile(res, sprintf('reproduce_70mN_%dx%d', nD, nA)));
adoptArcs   = pick(opts, 'adoptArcs', true);
discover    = pick(opts, 'discover', false);
nWorkers    = pick(opts, 'nWorkers', 4);
maxRounds   = pick(opts, 'maxRounds', 3);
compareOnly = pick(opts, 'compareOnly', false);
say         = pick(opts, 'print', true);

recordDir = fullfile(res, 'library_70mN_24x24_final');               % what the rebuild is measured against
recordCat = fullfile(recordDir, 'costate_catalog_dro_tulip_70mN.mat');
tag = '70mN';

%% ========================================================================
%  1. THE PROBLEM AND THE GRID -- the same as the record's
%% ========================================================================
sA0    = 0.0754;                                   % the first anchor's arrival phase: the grid's origin
engine = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150);
spec = struct( ...
    'sD',        (0:nD-1)/nD, ...                                  % nD departure phases; sD(1) = 0 is the spine
    'sA',        sort(mod(sA0 + (0:nA-1)/nA, 1)), ...              % nA arrival phases, a lattice through the first anchor's
    'departure', struct('family', 'dro', 'tau', 1.0), ...          % the DRO of period 1.0 ND (4.43 d)
    'arrival',   struct('family', 'tulip', 'Np', 7, 'pm', -1), ... % the 7-petal southern tulip
    'engine',    engine, ...
    'outDir',    outDir, 'tag', tag, ...
    'nWorkers',  nWorkers, 'maxRounds', maxRounds, ...
    'discover',  discover, ...
    'librarySeeds', true);                         % the sheet is also seeded from the shipped certified roots, as the record's was

%% ========================================================================
%  2. THE FIVE FAMILIES -- one certified root each, all at departure phase 0.
%     Every arrival arc starts from one of these; the table is the chain
%     script's (build_70mN_library section 0), with where each came from.
%% ========================================================================
spec.anchors = { ...
  % name        certified root (.mat)                                       arrival phase    family label
    'anchor',   fullfile(res, 'mintime_70mN_anchor.mat'),                   sA0,             'fast';     ... % the 2026-09-09 family
    'cell11',   fullfile(res, 'mintime_70mN_certified.mat'),                sA0 + 10/12,     'A2';       ... % the 26 d family
    'fast2',    fullfile(res, 'mintime_70mN_anchor_fast2.mat'),             0.8671,          'fast2';    ... % 2026-09-13, faster near 0.87
    'direct18', fullfile(res, 'mintime_70mN_anchor_direct18.mat'),          sA0 + 17/24,     'direct18'; ... % found by a direct solve at 0.7837
    'direct11', fullfile(res, 'mintime_70mN_anchor_direct11.mat'),          sA0 + 10/24,     'direct11'};    % found by a direct solve at 0.4921
names = spec.anchors(:, 1);
for k = 1:numel(names)
    assert(compareOnly || isfile(spec.anchors{k, 2}), 'reproduce_library_70mN:anchor', ...
           'the %s anchor is missing: %s', names{k}, spec.anchors{k, 2});      % (a comparison needs no anchor)
end

%% ========================================================================
%  3. THE ARCS -- two per family (arrival phase walked down and up)
%% ========================================================================
walked = cellfun(@(n, dn) fullfile(res, sprintf('arrival_arc_%s_%s_long.mat', n, dn)), ...
                 repelem(names, 2), repmat({'dn'; 'up'}, numel(names), 1), 'UniformOutput', false);
nArcsFound = nnz(cellfun(@isfile, walked));

% THE SEED ROOTS. Besides its arcs, the record's sheet was seeded with
% certified roots found by direct solves at single grid phases (FINDINGS 63,
% 68) -- some at arrival phases no family's arc reaches. They live in one
% registry file; the campaign adopts it as its own registry, which the driver
% then reads for every sheet and appends to.
seedRegistry = fullfile(res, 'mintime_70mN_direct_certified.mat');
nRegistryRoots = 0;
if isfile(seedRegistry), Lr = load(seedRegistry);  nRegistryRoots = numel(Lr.direct); end

out = struct('state', 'planned', 'spec', spec, 'nArcsFound', nArcsFound, 'nRegistryRoots', nRegistryRoots, ...
             'torus', [], 'catalog', '', 'comparison', [], 'audit', [], 'ok', false);
if say, printPlan(spec, names, walked, adoptArcs, go, compareOnly, recordCat, nRegistryRoots); end
if ~go && ~compareOnly
    return                                         % PLAN ONLY: nothing was created
end
assert(isfile(recordCat), 'reproduce_library_70mN:record', 'the library of record is missing: %s', recordCat);

if ~compareOnly
    if adoptArcs
        nCopied = adopt_walked_arcs(res, names, fullfile(outDir, 'arcs'), tag);
        fprintf('3. arcs: %d adopted now, %d already in the campaign folder\n', nCopied, 2*numel(names) - nCopied);
    else
        % "re-walk" is only true of a folder that holds no arcs: the driver walks
        % an arc only when its file is ABSENT, so arcs left by an earlier call
        % would be used as they are
        have = dir(fullfile(outDir, 'arcs', sprintf('arrival_arc_%s_*_long.mat', tag)));
        assert(isempty(have), 'reproduce_library_70mN:arcs', ...
               ['.adoptArcs = false asks for the arcs to be walked again, but %s already holds %d arc(s), ' ...
                'which the driver would reuse. Use a fresh .outDir.'], fullfile(outDir, 'arcs'), numel(have));
        fprintf('3. arcs: NOT adopted -- run_phase_torus walks all %d (about 2-6 h each)\n', 2*numel(names));
    end
    campaignRegistry = fullfile(outDir, 'direct_certified.mat');
    if isfile(seedRegistry) && ~isfile(campaignRegistry)     % never overwritten: the driver appends to it
        if ~isfolder(outDir), mkdir(outDir); end
        copyfile(seedRegistry, campaignRegistry);
        fprintf('3. seed roots: %d adopted into the campaign''s registry\n', nRegistryRoots);
    end

    %% ====================================================================
    %  4. THE BUILD -- one call. Each round: sheet at the nA arrival phases
    %     from every arc, ribs down the nD - 1 other departure phases for every
    %     column, package + audit + second-order sweep, then direct solves
    %     for the holes. It stops at a fixed point (a round that registers no
    %     new root and adds no anchor).
    %% ====================================================================
    out.torus = run_phase_torus(spec);
    out.state = out.torus.state;
end

%% ========================================================================
%  5. THE COMPARISON -- the rebuilt catalog against the record
%% ========================================================================
out.catalog = fullfile(outDir, 'final', 'costate_catalog_dro_tulip_70mN.mat');
assert(isfile(out.catalog), 'reproduce_library_70mN:noCatalog', ...
       'no rebuilt catalog at %s (did the build finish?)', out.catalog);
fprintf('\n5. ');
out.comparison = compare_phase_catalogs(out.catalog, recordCat);
if compareOnly, out.state = 'compared'; end

%% ========================================================================
%  6. THE VERDICT -- the comparison, and the build's own fail-closed audit
%     (every round's catalog is audited before run_phase_torus accepts it;
%     the last round's count is read back here so it is on the page)
%% ========================================================================
Lc = load(out.catalog);  fnc = fieldnames(Lc);  nEntries = Lc.(fnc{1}).n_entries;
[out.audit, whyAudit] = finalAudit(outDir, nEntries);
fprintf('\n6. VERDICT\n');
fprintf('   campaign status       : %s\n', out.state);
if ~isempty(out.audit)
    fprintf('   audit (fail-closed)   : %d ok / %d bad, covering all %d entries of the final catalog\n', out.audit.nOk, out.audit.nBad, nEntries);
else
    fprintf('   audit (fail-closed)   : NOT ESTABLISHED -- %s\n', whyAudit);
end
fprintf('   same library as record: %s  (%d of %d shared cells agree in t_f, costates and family%s)\n', ...
        passFail(out.comparison.ok), out.comparison.nAgree, out.comparison.nBoth, ...
        ternary(out.comparison.sameGrid, '', sprintf('; %d new cells lie off the record''s grid', out.comparison.nNewOffGrid)));
out.ok = out.comparison.ok && ~isempty(out.audit) && out.audit.nBad == 0;
fprintf('   REPRODUCED            : %s\n', passFail(out.ok));
end

% ==========================================================================
function printPlan(spec, names, walked, adoptArcs, go, compareOnly, recordCat, nRegistryRoots)
% PRINTPLAN  Sections 0-3 as text: what would be built, from what, and how
% long it takes.  INPUTS: as named.  OUTPUTS: none.
fprintf('REPRODUCE THE 70 mN DRO -> TULIP LIBRARY%s\n', ternary(go || compareOnly, '', '   [PLAN ONLY -- pass struct(''go'', true) to build]'));
fprintf('  problem   : DRO tau %.1f ND -> %d-petal tulip (branch %+d); %.0f mN, Isp %g s, %g kg\n', ...
        spec.departure.tau, spec.arrival.Np, spec.arrival.pm, spec.engine.thrustN*1000, spec.engine.ispS, spec.engine.m0kg);
fprintf('  grid      : %d departure x %d arrival phases = %d cells\n', numel(spec.sD), numel(spec.sA), numel(spec.sD)*numel(spec.sA));
fprintf('  families  :\n');
for k = 1:numel(names)
    found = nnz(cellfun(@isfile, walked(2*k-1:2*k)));
    fprintf('     %-9s sA %.4f  (%-8s)  walked arcs on disk: %d of 2\n', names{k}, mod(spec.anchors{k, 3}, 1), spec.anchors{k, 4}, found);
end
nFound = nnz(cellfun(@isfile, walked));
if adoptArcs
    fprintf('  arcs      : ADOPT the %d walked arcs (%d missing would stop the run)\n', nFound, numel(walked) - nFound);
else
    fprintf('  arcs      : RE-WALK all %d (2-6 h each)\n', numel(walked));
end
hrs = estimateHours(numel(spec.sD), numel(spec.sA), spec.nWorkers, adoptArcs);
fprintf('  time      : about %.0f h (%.1f days) with %d rib workers, from the rates measured on the 24 x 24 rebuild\n', hrs, hrs/24, spec.nWorkers);
fprintf('  seed roots: %d certified direct-found roots adopted into the campaign''s registry\n', nRegistryRoots);
fprintf('  discovery : %s\n', ternary(spec.discover, 'ON (probes at empty and slow columns)', 'off (all five families are given)'));
fprintf('  output    : %s\n', spec.outDir);
fprintf('  record    : %s\n', recordCat);
end

% ==========================================================================
function hrs = estimateHours(nD, nA, nWorkers, adoptArcs)
% ESTIMATEHOURS  How long a build takes, from the rates MEASURED on the
% 24 x 24 rebuild of 2026-09-18/19 (4 rib workers, 24 h 12 min in all):
%   sheet       1 h 34 for 24 arrival phases           -> scales with nA
%   ribs        235 worker-seconds per rib point       -> nA (nD - 1) points / workers
%   finalizer   28 s per entry (package, audit, sweep) -> nD nA entries
%   filler      10 % of the cells are holes, 7.2 min each
%   re-package  15.7 s per entry (the audit runs over the whole catalog again)
%   arcs        ten arcs, about 2 h each on two at a time, when they are walked
% The first estimate quoted for this script was "6-10 h". It was a guess.
% INPUTS: nD; nA; nWorkers; adoptArcs.  OUTPUTS: hrs.
sheet    = 1.57*nA/24;
ribs     = nA*max(nD - 1, 0)*235/max(nWorkers, 1)/3600;
finalize = nD*nA*28/3600;
filler   = 0.10*nD*nA*7.2/60;
repack   = nD*nA*15.7/3600;
arcs     = 20*(~adoptArcs);
hrs = sheet + ribs + finalize + filler + repack + arcs;
end

% ==========================================================================
function [A, why] = finalAudit(outDir, nEntries)
% FINALAUDIT  The audit OF THE FINAL CATALOG, or [] with the reason. `final`
% is copied from the LAST round in the driver's state, so that round's
% audit file is the only one that counts: an earlier round's, or one that
% covers fewer entries than the catalog holds, says nothing about the library
% being judged (it used to be "the newest audit file found anywhere"). The file
% is found by PATTERN: the packaging chain names it audit_<its own tag>.mat,
% which is '70mN' whatever tag the driver was given (seen in the 3 x 3
% rehearsal, where the driver's tag was 'test3d').
% INPUTS: outDir; nEntries (entries in the final catalog).
% OUTPUTS: A (audit_phase_catalog output) | []; why (char).
A = [];
stateF = fullfile(outDir, 'torus_state.mat');
if ~isfile(stateF), why = 'the campaign has no state file';  return, end
L = load(stateF);
if isempty(L.st.rounds), why = 'the campaign finished no round';  return, end
lastRound = L.st.rounds(end).dir;
found = dir(fullfile(lastRound, 'audit_*.mat'));
if numel(found) ~= 1
    why = sprintf('the last round (%s) holds %d audit file(s), not exactly one', lastRound, numel(found));
    return
end
La = load(fullfile(found.folder, found.name));
if La.A.nOk + La.A.nBad ~= nEntries
    why = sprintf('the last round''s audit covers %d entries, the final catalog holds %d', La.A.nOk + La.A.nBad, nEntries);
    return
end
A = La.A;  why = 'the last round''s audit, covering every entry';
end

function H = localHandles(fh)
% LOCALHANDLES  This file's local functions as a struct of handles keyed by
% name -- the TEST SEAM.  INPUTS: fh (cell of handles).  OUTPUTS: H struct.
H = struct();
for k = 1:numel(fh), H.(func2str(fh{k})) = fh{k}; end
end

% ==========================================================================
function v = pick(s, f, d_)
% PICK  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function t = passFail(c)
% PASSFAIL  'PASS' or 'FAIL'.  INPUTS: c (logical).  OUTPUTS: t (char).
if c, t = 'PASS'; else, t = 'FAIL'; end
end

function v = ternary(c, a, b)
% TERNARY  a when c, else b.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end
