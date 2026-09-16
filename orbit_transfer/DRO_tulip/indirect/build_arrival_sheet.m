function S = build_arrival_sheet(opts)
%% Purpose:
%
%   Front door for the ARRIVAL-PHASE sheet at one operating point: gather
%   every arclength arc saved in results/ (arrival_arc_*.mat), certify all
%   grid crossings with certify_crossing, assemble the sheet with
%   sheet_from_arcs, save it, and print the table (grid phase, minimum
%   certified t_f, dV, fuel, number of candidates / certified, first
%   failure reason when nothing certified).
%
%% Inputs:
%
%  opts                     struct (optional)
%   .pattern ['arrival_arc_*.mat'] .arcDir [results/] .out ['results/arrival_sheet_70mN.mat']
%   .sA [] an explicit list of arrival phases (any spacing) -- the sheet's
%   columns; else the lattice .sA0 [0.0754] + (0:nA-1)/nA, .nA [12]
%   .copts (certify_crossing options)
%   .seedFiles {} .mat files holding a `direct` struct array (sD sA tfDays
%   z src) of certified roots to seed with; .librarySeeds [true] also seed
%   from dro_tulip_library (false for a campaign at another operating point)
%   .rescan [true] take the arcs' crossings from a RE-SCAN of their stored
%   path at this sheet's levels, instead of the crossings each walk happened
%   to record. An arc records crossings of the levels it was GIVEN, so
%   refining the grid would otherwise mean walking the arcs again -- hours
%   each -- when the path they already stored contains every crossing there
%   is (crossings_from_arc). The re-scan reproduces the recorded crossings
%   exactly and adds the arcs' own start points, which duplicate the library
%   seeds below and are merged by sheet_from_arcs.
%   plus arclength_arrival setup options (.thrustN .ispS .sD ...)
%
%% Outputs:
%
%  S                        struct                  sheet_from_arcs output +
%                                                   .arcs (file list) .B
%                                                   .anc .opts .built
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
pat = d('pattern', 'arrival_arc_*.mat');
out = d('out', fullfile(here, 'results', 'arrival_sheet_70mN.mat'));
lStar = 389703.264829278;  tStar = 382981.289129055;

% the fence's pool, created AFTER startup has set the path (workers inherit
% the client path at pool creation)
pool = capped_pool();
[B, anc] = arclength_arrival('setup', opts);
arcDir = d('arcDir', fullfile(here, 'results'));   % a campaign may keep its own arcs
files = dir(fullfile(arcDir, pat));
% a walk in progress leaves <arc>.partial.mat beside the finished arcs
% (arclength_ms .partialFile); it is not an arc of record and would enter
% the sheet twice once the walk finishes
files = files(~contains({files.name}, '.partial.'));
assert(~isempty(files), 'no arcs match %s', pat);
arcs = cell(1, numel(files));
nAwant = d('nA', 12);  sA0want = d('sA0', 0.0754);
if isfield(opts, 'sA') && ~isempty(opts.sA)
    sAlist = mod(opts.sA(:).', 1);  nAwant = numel(sAlist);  sA0want = sAlist(1);
else
    sAlist = mod(sA0want + (0:nAwant-1)/nAwant, 1);
end
rescan = d('rescan', true);
levels = sAlist;                            % crossings_from_arc scans each level mod 1
for k = 1:numel(files)
    L = load(fullfile(files(k).folder, files(k).name));
    A = L.A;
    nRec = numel(A.crossings);
    if rescan
        A.crossings = crossings_from_arc(A, levels, B);
        A.rescannedAt = levels;
    end
    arcs{k} = A;
    fprintf('arc %d: %-32s %4d roots, sA %.4f -> %.4f, %d folds, %d crossings%s, stop = %s\n', ...
        k, files(k).name, numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), ...
        pick(rescan, sprintf(' (re-scanned at %d levels; the walk recorded %d)', nAwant, nRec), ''), ...
        A.stop);
end

% SEED the sheet with the certified library solutions at this departure
% phase: they are the arcs' start points, so they are not crossings of any
% arc and a crossings-only sheet reports NaN at exactly the phases whose
% solutions launched it.
sD0 = anc.sD;
policy = d('copts', struct());
policy.pool = pool;
if ~isfield(policy, 'wallSec'), policy.wallSec = 600; end
lib = struct('sD', {}, 'sA', {}, 'tfDays', {}, 'src', {}, 'z', {}, 'Y', {});
if d('librarySeeds', true)
    lib0 = dro_tulip_library(here);
    for k = 1:numel(lib0), lib(end+1) = seedRow(lib0(k)); end
end
for sf = d('seedFiles', {})
    Lf = load(sf{1});
    assert(isfield(Lf, 'direct'), 'seed file %s holds no `direct` struct array', sf{1});
    for k = 1:numel(Lf.direct), lib(end+1) = seedRow(Lf.direct(k)); end
end
if ~isempty(lib), lib = lib(abs(mod([lib.sD] - sD0 + 0.5, 1) - 0.5) < 1e-8); end
seeds = struct([]);
for k = 1:numel(lib)
    if min(abs(mod(sAlist - lib(k).sA + 0.5, 1) - 0.5)) > 1e-6, continue, end   % off-grid
    rv0 = B.stateD(sD0);
    % THE shared builder (costate_common/seed_from_entry)
    seed = seed_from_entry(lib(k), rv0(1:6), ...
                           struct('Tnd', B.Tnd, 'cnd', B.cnd, 'muStar', B.mu, 'K', 24));
    % ONE certification policy for both routes. Seeds used to be certified
    % with a fresh hardcoded struct while opts.copts applied only to
    % crossings, so a stricter requested gate silently did not reach the
    % seeds that populate the same sheet. (Astra chain review 2026-09-10.)
    C = certify_root(seed, rv0, B.stateA(lib(k).sA), B, ...
                     setfield(setfield(policy, 'sA', lib(k).sA), 'sD', sD0)); %#ok<SFLD>
    C.note = join_note(sprintf('seed: %s (%.3f d)', lib(k).src, lib(k).tfDays), C);
    fprintf('library seed (%.4f, %.4f) [%s]: %s\n', sD0, lib(k).sA, lib(k).src, C.reason);
    if isempty(seeds), seeds = C; else, seeds(end+1) = C; end %#ok<AGROW>
end

S = sheet_from_arcs(arcs, struct('sA', sAlist, 'seeds', seeds, ...
                                 'B', B, 'anc', anc, 'copts', policy));
S.arcs = {files.name};  S.B = B;  S.anc = anc;  S.opts = opts;  S.built = datestr(now);
S.problem = B.problem;          % the identity packaging must use
S.policy = rmfield(policy, 'pool');
S.B = rmfield(S.B, {'res', 'dRdq', 'stateA', 'stateD'});
% published in one rename: an interrupted save must not leave a sheet
% file that the chain would then reuse as if it were complete
tmpS = sprintf('%s.%s.part', out, char(java.util.UUID.randomUUID()));
save(tmpS, 'S');
publish_atomic(tmpS, out);

fprintf('\n  j    sA      t_f [d]   dV [km/s]  fuel [kg]  cand  cert   note\n');
for j = 1:numel(S.sA)
    c = S.cand{j};  nc = numel(c);  ncert = 0;  note = '';
    if nc > 0, ncert = nnz([c.ok]); end
    if isfinite(S.TF(j))
        k = find([c.ok] & abs([c.tfDays] - S.TF(j)) < 1e-9, 1);
        fprintf(' %2d  %.4f  %9.4f  %9.4f  %9.2f   %2d    %2d\n', j, S.sA(j), S.TF(j), c(k).dvKms, c(k).propellantKg, nc, ncert);
    else
        if nc > 0, note = c(1).reason; else, note = 'no crossing'; end
        fprintf(' %2d  %.4f        ---        ---        ---   %2d    %2d   %s\n', j, S.sA(j), nc, ncert, note);
    end
end
fprintf('saved %s\n', out);
end

function s = join_note(prefix, C)
% JOIN_NOTE  The entry's provenance note: this producer's clause in front
% of whatever the certifier (or an earlier producer) already wrote.
% INPUTS: prefix (char); C (struct, .note optional).  OUTPUTS: s.
parts = {prefix};
if isfield(C, 'note') && ~isempty(C.note), parts{end+1} = C.note; end
s = strjoin(parts, ' | ');
end

function r = seedRow(e)
% SEEDROW  One seed in the library's row layout (Y optional).  INPUTS: e.
% OUTPUTS: r.
r = struct('sD', e.sD, 'sA', e.sA, 'tfDays', e.tfDays, 'src', e.src, 'z', e.z(:), 'Y', []);
if isfield(e, 'Y'), r.Y = e.Y; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

% ------------------------------------------------------------------------
function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end
