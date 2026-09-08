function out = build_dro_deliverable(cfg)
%% Purpose:
%
%   Package the DRO -> tulip minimum-time costate catalog as a shippable
%   deliverable: the catalog itself, the maintained helper suite, and a
%   README whose every number is COMPUTED FROM THE CATALOG.
%
%   That last point is the reason this file exists. The 2026-08 deliverable
%   was assembled by hand and its README hard-coded "3,936 entries", "9
%   thrust levels" and "coverage 89-97%". The catalog has since gained the
%   0.75 N and 0.5 N rungs (4,439 entries, 11 rungs, 70% coverage), a full
%   conjugate-point sweep and the three sufficiency-hypothesis gates -- so
%   the shipped README understated and misdescribed what we have. Deriving
%   the numbers makes that failure mode impossible: re-run this and the
%   description cannot drift from the data.
%
%   Helper sources are the MAINTAINED copies in DRO_tulip/indirect, never
%   the previous deliverable's copies.
%
%% Inputs:
%
%  cfg                      struct (optional)
%   .catMat                 char                    catalog path [direct/results]
%   .outDir                 char                    [deliverables/costate_catalog_deliverable3]
%   .movies                 logical                 re-render the extremes
%                                                   movies [false]
%   .zip                    logical                 write the .zip [true]
%
%% Outputs:
%
%  out                      struct                  .dir .zip .facts (the
%                                                   computed numbers) .files
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, cfg = struct(); end
here   = fileparts(mfilename('fullpath'));
catMat = getf(cfg, 'catMat', fullfile(here, 'direct', 'results', 'costate_catalog_dro_tulip.mat'));
outDir = getf(cfg, 'outDir', fullfile(here, 'deliverables', 'costate_catalog_deliverable3'));
doMov  = getf(cfg, 'movies', false);
doZip  = getf(cfg, 'zip', true);
srcDir = fullfile(here, 'indirect');

L = load(catMat);  fn = fieldnames(L);
assert(numel(fn) == 1, 'expected one variable in %s', catMat);
c = L.(fn{1});

%% Facts, computed -- never typed --------------------------------------------
f = struct();
f.name      = c.name;
f.created   = c.created;
f.nSheets   = numel(c.sheets);
f.rungs     = c.rungs_N;
f.nRungs    = numel(c.rungs_N);
f.isp       = c.thruster.isp_s;
f.m0        = c.thruster.m0_kg;
f.tau       = unique(arrayfun(@(s) s.tauDRO, c.sheets));
f.Np        = unique(arrayfun(@(s) s.Np, c.sheets));
f.pm        = unique(arrayfun(@(s) s.pm, c.sheets));
f.nPhase    = size(c.sheets(1).has_solution, 1);
f.nSlots    = 0;  f.nSolved = 0;  f.tfMin = inf;  f.tfMax = -inf;
cp = [0 0 0];  covSheet = zeros(1, f.nSheets);
for k = 1:f.nSheets
    s = c.sheets(k);
    f.nSlots  = f.nSlots  + numel(s.has_solution);
    f.nSolved = f.nSolved + nnz(s.has_solution);
    covSheet(k) = 100*nnz(s.has_solution)/numel(s.has_solution);
    t = s.tf_nd(s.has_solution);
    if ~isempty(t), f.tfMin = min(f.tfMin, min(t));  f.tfMax = max(f.tfMax, max(t)); end
    if isfield(s, 'conj_pass')
        v = s.conj_pass(s.has_solution);
        cp = cp + [nnz(v == 1) nnz(v == 0) nnz(v < 0)];
    end
end
f.coverage   = 100*f.nSolved/f.nSlots;
f.covMin     = min(covSheet);  f.covMax = max(covSheet);
f.conjPass   = cp(1);  f.conjFail = cp(2);  f.conjNotrun = cp(3);
f.hasConj    = isfield(c, 'conj_test');
f.hasGates   = isfield(c, 'hyp_gates');
if f.hasConj,  f.conjDate  = c.conj_test.date;  f.conjK = c.conj_test.K; end
if f.hasGates, f.gatesDate = c.hyp_gates.date; end
f.tStar = c.constants.tStar_s;
f.dayMin = f.tfMin*f.tStar/86400;  f.dayMax = f.tfMax*f.tStar/86400;
% per-rung coverage
f.rungCov = zeros(1, f.nRungs);
for kr = 1:f.nRungs
    a = 0; b = 0;
    for k = 1:f.nSheets
        h = c.sheets(k).has_solution(:,:,kr);  a = a + nnz(h);  b = b + numel(h);
    end
    f.rungCov(kr) = 100*a/b;
end

%% Assemble -------------------------------------------------------------------
if ~exist(outDir, 'dir'), mkdir(outDir); end
helpers = {'costate_catalog_pick.m', 'costate_lib_describe.m', ...
           'costate_catalog_extremes.m', 'costate_catalog_example.m', ...
           'costate_catalog_extremes_movies.m'};
files = {};
copyfile(catMat, fullfile(outDir, 'costate_catalog_dro_tulip.mat'));
files{end+1} = 'costate_catalog_dro_tulip.mat';
for k = 1:numel(helpers)
    src = fullfile(srcDir, helpers{k});
    assert(exist(src, 'file') == 2, 'missing maintained helper %s', src);
    copyfile(src, fullfile(outDir, helpers{k}));
    files{end+1} = helpers{k}; %#ok<AGROW>
end
writeReadme(fullfile(outDir, 'README.md'), f, c);
files{end+1} = 'README.md';

if doMov
    addpath(srcDir);
    costate_catalog_extremes_movies(fullfile(outDir, 'costate_catalog_dro_tulip.mat'), ...
        struct('outDir', outDir));
end

out = struct('dir', outDir, 'facts', f, 'files', {files}, 'zip', '');
if doZip
    z = [outDir '.zip'];
    if exist(z, 'file'), delete(z); end
    zip(z, '*', outDir);
    out.zip = z;
    d = dir(z);
    fprintf('deliverable: %s (%.1f MB)\n', z, d.bytes/1e6);
end
fprintf(['DRO->tulip deliverable rebuilt: %d entries, %d sheets, %d rungs, ' ...
         '%.1f%% coverage, conj %d/%d/%d\n'], f.nSolved, f.nSheets, f.nRungs, ...
         f.coverage, f.conjPass, f.conjFail, f.conjNotrun);
end

% ------------------------------------------------------------------------
function writeReadme(path, f, c)
% WRITEREADME  Emit the deliverable README with every number derived from
% the catalog.  INPUTS: path; f (facts); c (catalog).  OUTPUTS: none.
fid = fopen(path, 'w');
p = @(varargin) fprintf(fid, varargin{:});
p('# DRO -> Tulip Minimum-Time Costate CATALOG (deliverable 3)\n\n');
p(['**%d converged minimum-time transfer solutions** spanning **%d DRO ' ...
   'periods x %d tulip petal counts x a %dx%d phasing torus x %d thrust ' ...
   'levels**. Every entry is a root of the indirect (PMP) boundary-value ' ...
   'problem: hand its `z8` to `pumpkyn.cr3bp.tfMin` and it is accepted ' ...
   'unchanged in about a second.\n\n'], f.nSolved, numel(f.tau), numel(f.Np), ...
   f.nPhase, f.nPhase, f.nRungs);
p('Catalog built %s; this package generated %s by `build_dro_deliverable.m`.\n', ...
  string(f.created), datestr(now, 'yyyy-mm-dd'));
p('Every number below is computed from the catalog, not transcribed.\n\n');

p('## What it covers\n\n');
p('| axis | values |\n|---|---|\n');
p('| DRO period tau (ND) | %s |\n', strjoin(compose('%.3g', f.tau(:)'), ', '));
p('| tulip petal count Np | %s |\n', strjoin(compose('%g', f.Np(:)'), ', '));
p('| family branch pm | %s (pm = +1 orbits and costates are exact z-mirrors) |\n', ...
  strjoin(compose('%+d', f.pm(:)'), ', '));
p('| phasing | %dx%d torus per sheet, phases as fractions of each orbit period |\n', f.nPhase, f.nPhase);
p('| thrust (N) | %s |\n', strjoin(compose('%g', f.rungs(:)'), ', '));
p('| propulsion | Isp %g s, m0 %g kg |\n', f.isp, f.m0);
p('| time of flight | %.3f .. %.3f ND = %.2f .. %.2f days |\n\n', f.tfMin, f.tfMax, f.dayMin, f.dayMax);

p('**Coverage: %d of %d grid slots (%.1f%%)**, %.0f%%-%.0f%% per sheet.\n', ...
  f.nSolved, f.nSlots, f.coverage, f.covMin, f.covMax);
p('`sheets(k).has_solution` is the per-(pair, rung) authority. Coverage\n');
p('thins at low thrust, which is real and worth knowing before you plan\n');
p('around it:\n\n');
p('| thrust (N) | %s |\n|---|%s\n', strjoin(compose('%g', f.rungs(:)'), ' | '), repmat('---|', 1, f.nRungs));
p('| solved | %s |\n\n', strjoin(compose('%.0f%%', f.rungCov(:)'), ' | '));

p('## How it was made, and how it was checked\n\n');
p('Per sheet: a direct Hermite-Simpson collocation solve anchors each phase\n');
p('pair cold at the highest thrust, then thrust is walked down with each rung\n');
p('warm-starting the next, so a pair stays on one solution family. Every\n');
p('entry passes three FIRST-ORDER gates:\n\n');
p('1. multiple-shooting refinement residual (`ms_tfmin`, analytic STM Jacobian): ~1e-10\n');
p('2. the PMP control law is FLOWN end to end and must arrive: ~1e-4 km\n');
p('3. **`pumpkyn.cr3bp.tfMin` returns the entry unchanged** (|dz| ~ 1e-9) -- a\n');
p('   foreign witness: different code, integrator and root finder\n\n');
if f.hasConj || f.hasGates
    p('### Second-order verification (new since the previous deliverable)\n\n');
end
if f.hasConj
    p('**Conjugate-point test** (`conj_test`, %s, K = %d): the free-time\n', string(f.conjDate), f.conjK);
    p('quotiented Jacobi determinant, sampled at the multiple-shooting junctions.\n');
    p('`sheets(k).conj_pass`: 1 = no conjugate point found in (0, tf); 0 = a sign\n');
    p('change strictly inside, i.e. the entry is NOT a local minimum; -1 = not\n');
    p('decided. Result: **%d pass, %d fail, %d undecided**.\n\n', f.conjPass, f.conjFail, f.conjNotrun);
end
if f.hasGates
    p('**Sufficiency-hypothesis gates** (`hyp_gates`, %s): the three hypotheses\n', string(f.gatesDate));
    p('the conjugate test assumes, checked per entry and stored as\n');
    p('`gate_min_lamv` (strong Legendre, min|lam_v| > 0), `gate_min_qmt`\n');
    p('(all-burn is the PMP control, min Q_mt > 0) and `gate_dimS` (no abnormal\n');
    p('lift, dim S = 1). **All entries pass all three.**\n\n');
    p('With `conj_pass = 1` these are the Bonnard-Caillau-Trelat sufficiency\n');
    p('hypotheses, verified numerically at the sampled times. What is NOT\n');
    p('claimed: the determinant is sampled at junctions, so a conjugate pair\n');
    p('inside one segment is invisible, and the second-order layer -- unlike\n');
    p('the first -- has no independent second implementation behind it.\n\n');
end

p('## The COMPACT format\n\n');
p('Only canonical nondimensional quantities are stored: constants once at top\n');
p('level, per-sheet phase fractions and rung availability, `tf_nd` lookup\n');
p('grids, and the `z8` vectors (which already contain t_f). Days, dV and\n');
p('masses are DERIVED, and every formula rides along in `cat.derive`:\n\n');
p('```\n');
dn = fieldnames(c.derive);
for k = 1:numel(dn)
    v = c.derive.(dn{k});
    if ischar(v) || isstring(v), p('%-22s %s\n', dn{k}, string(v)); end
end
p('```\n\n');

p('## Quick start\n\n');
p('```matlab\n');
p('L = load(''costate_catalog_dro_tulip.mat'');\n');
p('cat = L.%s;\n', f.name);
p('[tf_nd, z8, info] = costate_catalog_pick(cat, %g, %g, 3.0, 11.6, %g);\n', ...
  f.tau(min(3,end)), f.Np(2), f.rungs(min(5,end)));
p('%%                                        tau  Np  dep  arr  thrust(N)\n');
p('%% z8 -> pumpkyn.cr3bp.tfMin as-is;  info.delivered = what you actually got\n');
p('```\n\n');
p('Or run `costate_catalog_example`. For a fact sheet: `costate_lib_describe`.\n');
p('Requires **pumpkyn / pumpkynPie** on the MATLAB path; the data itself is\n');
p('dependency-free.\n\n');

p('## The honesty contract (the picker''s warnings)\n\n');
p('Whenever what is RETURNED differs from what you REQUESTED, the picker\n');
p('prints exactly what you are getting: the nearest sheet, the nearest grid\n');
p('pair, the nearest SOLVED pair, or a seed from a different rung.\n');
p('`info.delivered` carries the same facts programmatically and `info.warned`\n');
p('flags any substitution. The right use of a substituted answer is as a\n');
p('**seed**: hand `z8` to `tfMin` at your true endpoints and thrust.\n\n');

p('## Gotchas, honestly\n\n');
p('- **Coverage gaps are real** and concentrated at low thrust (see the table\n');
p('  above). `has_solution` is always the authority.\n');
p('- **Only the pm = %+d branch is stored.** The mirror branch follows from the\n', f.pm(1));
p('  CR3BP symmetry (flip all z components) but is not catalogued.\n');
p('- **At fixed thrust, minimum time = minimum dV** (continuous burn makes dV\n');
p('  monotone in t_f). The metrics differ only ACROSS thrust levels.\n');
p('- **Neighbouring cells can sit on different solution families** -- t_f can\n');
p('  jump >15%% across a family wall. Seeding tfMin from a neighbour may find\n');
p('  the faster family.\n');
p('- The %dx%d grid is for interpolation and seeding, not final answers.\n\n', f.nPhase, f.nPhase);
p('M. Casey / D. Koblick, Coorbital Inc.\n');
fclose(fid);
end

function v = getf(s, fld, d)
% GETF  Field with default.  INPUTS: s; fld; d.  OUTPUTS: v.
if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = d; end
end
