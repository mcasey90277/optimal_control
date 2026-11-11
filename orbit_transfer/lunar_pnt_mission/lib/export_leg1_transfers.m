function E = export_leg1_transfers(leg1Mat, opts)
%% Purpose:
%
%   CAPTURE THE LEG-1 TRANSFERS (GTO -> DRO slot, 5 N, minimum time) in a
%   SIDECAR file beside the leg-1 study results. Each accepted transfer of
%   the study's P1 becomes one entry with the fields a GTO -> DRO sheet
%   entry has (costates z8 = [lambda0(7); t_f], t_f, the thrust regime)
%   plus what makes it OFF-GRID: its true GTO orientation, departure phase,
%   arrival DRO phase, epochs and the endpoint STATES themselves (the
%   orientation and the arrival phase are not on the package grid, so the
%   phases alone do not rebuild the endpoints without the leg-1 handles).
%   Nothing is written into the GTO -> DRO package, and no certifier is
%   run here (standing rule 2026-11-10: no new certifier or schema).
%
%   STATUS HONESTY. Every entry records exactly what the study checked:
%   the flown arrival miss (fly_transfer from rv0 to rvf), the LOCATED
%   minimum Earth and Moon altitudes and their floors, the one-GTO and
%   hand-off gates, the second-order status ('pass', 'refuted',
%   'notReproduced', 'error') with the conjugate verdict and its covered
%   flag and the hypothesis gates (minLamV, minQmt, dimS) -- and what it
%   did NOT: no catalog acceptance by pumpkyn tfMin beyond the polish's own
%   ms_tfmin solve, no Earth-floor constraint inside the optimisation (the
%   clearance is checked after the fact), no pointwise PMP checks on a
%   tight flight, no admission to any library or package. SCOPE: the
%   second-order verdict and "minimum-time" refer to the fixed-endpoint
%   problem to the slot's arrival point (a local test at the root); the
%   moving-target rendezvous condition is not checked (said in the note
%   and first in every notChecked list).
%
%   Transfers with P1(k).accepted false are not exported; results with
%   none accepted throw export_leg1_transfers:noTransfers. An ACCEPTED
%   transfer whose stored gates did not pass (G.pass false) is refused:
%   export_leg1_transfers:notPassed, nothing written. The note is built
%   from the stored verdicts (counts per gate, statuses), not fixed text.
%
%   OUTPUT PATH: resolved to a canonical absolute path and REFUSED
%   (export_leg1_transfers:outFile, nothing created) unless it lies inside
%   this campaign's results/ folder. Written atomically (temporary file in
%   the same folder, then movefile); ONE variable, leg1_transfers (= E).
%
%   Self-demo (nargin == 0): the default results file, table printed.
%
%% Inputs:
%
%  leg1Mat                  char (optional)         leg-1 study results
%                                                   [results/lunar_pnt_leg1_
%                                                   study.mat]
%
%  opts                     struct (optional)       .outFile [results/
%                                                   leg1_transfers_gto_dro_
%                                                   5N.mat]
%                                                   .verbose [true]
%
%% Outputs:
%
%  E                        struct                  leg1_transfers:
%                                                   .entries [1 x n] (sat
%                                                   slot z8 [8 x 1] tf_nd
%                                                   rv0 rvf [6 x 1] thrustN
%                                                   ispS m0kg orientDeg sD
%                                                   arrPhase tDep tArr
%                                                   seededFrom [iD iA iS]
%                                                   nWait route flyKm flyVms
%                                                   verdicts checked
%                                                   notChecked)
%                                                   .regime .package
%                                                   .provenance .note .usage
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  fix round 1: refuse G.pass false, note from verdicts (c) 11/11/2026
%  M. Casey  final fix wave: scope sentence, canon_path         (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));
resDir = fullfile(here, 'results');
if nargin == 0
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    leg1Mat = fullfile(resDir, 'lunar_pnt_leg1_study.mat');
    opts = struct('verbose', true);
end
if nargin < 2, opts = struct(); end
for f = fieldnames(opts).'
    if ~any(strcmp(f{1}, {'outFile', 'verbose'}))
        error('export_leg1_transfers:opts', 'unknown option "%s"', f{1});
    end
end
outFile = fdef(opts, 'outFile', fullfile(resDir, 'leg1_transfers_gto_dro_5N.mat'));
verbose = fdef(opts, 'verbose', true);
tExport = tic;

%% Output path: inside results/ only
outFile = canon_path(outFile);
if ~startsWith(outFile, [canon_path(resDir) filesep])
    error('export_leg1_transfers:outFile', 'refusing to write outside %s: %s', canon_path(resDir), outFile);
end
[outDir, outName] = fileparts(outFile);

%% Read the results:
if ~isfile(leg1Mat)
    error('export_leg1_transfers:results', 'results file not found: %s', leg1Mat);
end
need = {'P1', 'T0', 'run_meta', 'S'};
w = whos('-file', leg1Mat);
miss = need(~ismember(need, {w.name}));
if ~isempty(miss)
    error('export_leg1_transfers:results', '%s lacks %s', leg1Mat, strjoin(miss, ', '));
end
R = load(leg1Mat, need{:});
acc = find([R.P1.accepted]);
if isempty(acc)
    error('export_leg1_transfers:noTransfers', '%s holds no accepted leg-1 transfer', leg1Mat);
end
bad = acc(arrayfun(@(q) ~(isstruct(q.G) && isequal(q.G.pass, true)), R.P1(acc)));
if ~isempty(bad)
    error('export_leg1_transfers:notPassed', ['satellite(s) %s are marked accepted but their stored gates ' ...
          'did not pass (G.pass false); refusing to export'], strtrim(sprintf('%d ', bad)));
end
rg = R.run_meta.regime;
dayNd = rg.tStar/86400;

%% One entry per accepted transfer:
entries = cell(1, numel(acc));
for ka = 1:numel(acc)
    k = acc(ka);
    q = R.P1(k);
    G = q.G;
    sec = G.second;
    v = struct('arrivalMissKm', G.arrival.missKm, 'arrivalMissMps', G.arrival.missMps, ...
               'arrivalPass', G.arrival.pass, 'flightAdmissible', G.flight.ok, ...
               'oneGtoDStateNd', G.oneGto.dStateNd, 'oneGtoPass', G.oneGto.pass, ...
               'handoffErrPhase', G.handoff.errPhase, 'handoffPass', G.handoff.pass, ...
               'earthMinKm', G.earth.hMinKm, 'earthTAtMinNd', G.earth.tAtMin, ...
               'earthFloorKm', G.earth.floorKm, 'earthPass', G.earth.pass, ...
               'moonMinKm', G.moon.hMinKm, 'moonTAtMinNd', G.moon.tAtMin, ...
               'moonFloorKm', G.moon.floorKm, 'moonPass', G.moon.pass, ...
               'secondStatus', sec.status, 'conjPass', sec.conjPass, 'conjCovered', sec.conjCovered, ...
               'conjVerdict', char(string(sec.conjVerdict)), 'minLamV', sec.minLamV, ...
               'minQmt', sec.minQmt, 'dimS', sec.dimS, 'reproducedDz', sec.dz, ...
               'secondReason', char(string(sec.reason)), 'gatesPass', G.pass, ...
               'refuted', G.refuted, 'secondUnknown', G.secondUnknown);
    checked = {sprintf('arrival: re-flown miss %.3g km / %.3g m/s (gate 100 km, 10 m/s)', ...
                       G.arrival.missKm, G.arrival.missMps), ...
               sprintf('flight admissible (fly_transfer validate_flight): %d', G.flight.ok), ...
               sprintf('one GTO: departure state recomputed from (orient0, gtoPhase0, tDep), diff %.2g ND', ...
                       G.oneGto.dStateNd), ...
               sprintf('hand-off: DRO phase of the flown arrival at T0 vs (k-1)/6, error %.2g of a period', ...
                       G.handoff.errPhase), ...
               sprintf('clearance on the LOCATED minimum: Earth %.1f km (floor %g), Moon %.1f km (floor %g)', ...
                       G.earth.hMinKm, G.earth.floorKm, G.moon.hMinKm, G.moon.floorKm), ...
               sprintf(['second order: ms_tfmin re-solve with the conjugate test (%s, conjPass %g, covered %g) ' ...
                        'and mintime_hypothesis_gates (minLamV %.3g, minQmt %.3g, dimS %g): status %s'], ...
                       char(string(sec.conjVerdict)), sec.conjPass, sec.conjCovered, sec.minLamV, ...
                       sec.minQmt, sec.dimS, sec.status)};
    notChecked = {['the moving-target rendezvous condition: the second-order verdict and "minimum-time" refer to ' ...
                   'the fixed-endpoint problem to the slot''s arrival point (a local test at the root)'], ...
                  'catalog acceptance by pumpkyn tfMin (only the polish''s own ms_tfmin solve and the second-order re-solve)', ...
                  'an Earth floor (or any path constraint) inside the optimisation: clearance is checked after the fact', ...
                  'pmp_pointwise_checks on a tight-tolerance flight, lift_margin, the H6 margin gate, a dense conjugate scan', ...
                  'admission to any library or package: this sidecar is not part of the GTO -> DRO package'};
    entries{ka} = struct('sat', k, 'slot', q.slot, 'z8', q.z8(:), 'tf_nd', q.tf, ...
                         'rv0', q.rv0(:), 'rvf', q.rvf(:), 'thrustN', rg.thrustN, 'ispS', rg.ispS, ...
                         'm0kg', rg.m0kg, 'orientDeg', q.orientDeg, 'sD', q.sD, 'arrPhase', q.arrPhase, ...
                         'tDep', q.tDep, 'tArr', q.tArr, 'tDepDays', q.tDep*dayNd, ...
                         'tArrDays', q.tArr*dayNd, 'seededFrom', [q.iD q.iA q.iS], 'nWait', q.n, ...
                         'route', q.route, 'flyKm', G.arrival.missKm, 'flyVms', G.arrival.missMps, ...
                         'verdicts', v, 'checked', {checked}, 'notChecked', {notChecked});
end
entries = [entries{:}];

%% Package, provenance, note:
idn = R.run_meta.identity;
package = struct('file', idn.file, 'sha256', idn.sha256, 'created', idn.created, ...
                 'n_entries', idn.n_entries, 'orientations', idn.orientations, 'rung', idn.rung);
provenance = struct('run_meta', R.run_meta, 'resultsFile', leg1Mat, ...
                    'exportDate', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
                    'exportFunction', 'export_leg1_transfers', 'T0', R.T0, ...
                    'orient0', R.S.orient0, 'gtoPhase0', R.S.gtoPhase0, 'delta', R.S.delta, ...
                    'exportWallS', NaN);
V = [entries.verdicts];
nE = numel(entries);
scope = ['The second-order verdict and "minimum-time" refer to the fixed-endpoint problem to the slot''s ' ...
         'arrival point (a local test at the root); the moving-target rendezvous condition is not checked.'];
cnt = @(f) nnz(arrayfun(@(v) isequal(v.(f), true), V));
st = {V.secondStatus};
note = sprintf(['%d OFF-GRID minimum-time GTO -> DRO transfers (5 N per %g kg, Isp %g s) captured from the ' ...
                'lunar_pnt_leg1_study run of %s (optimal_control %s%s), study verdict "%s". Each was polished ' ...
                'from a clean cell of the GTO -> DRO package (SHA-256 %s..., created %s, %d entries) to its ' ...
                'true GTO orientation and onto its exact DRO slot. Stored verdicts: flight admissible %d of %d, ' ...
                'arrival %d of %d (max miss %.3g km), one GTO %d of %d, hand-off at T0 %d of %d, Earth floor %d ' ...
                'of %d (min located %.1f km), Moon floor %d of %d (min located %.1f km), all gates %d of %d. ' ...
                'Second-order statuses: %s. %s NOT checked: catalog acceptance by pumpkyn tfMin, an Earth floor ' ...
                'inside the optimisation, tight-flight pointwise checks, lift_margin, H6, a dense conjugate ' ...
                'scan. Not admitted to any library; the package is not modified.'], ...
               nE, rg.m0kg, rg.ispS, R.run_meta.date, R.run_meta.optimal_control, ...
               tern(isequal(R.run_meta.optimal_controlDirty, true), ' (dirty)', ''), fdef(R.run_meta, 'verdict', '?'), ...
               idn.sha256(1:12), idn.created, idn.n_entries, cnt('flightAdmissible'), nE, cnt('arrivalPass'), nE, ...
               max([V.arrivalMissKm]), cnt('oneGtoPass'), nE, cnt('handoffPass'), nE, cnt('earthPass'), nE, ...
               min([V.earthMinKm]), cnt('moonPass'), nE, min([V.moonMinKm]), cnt('gatesPass'), nE, ...
               strjoin(st, ', '), scope);
usage = {'load(''leg1_transfers_gto_dro_5N.mat'');   % one variable: leg1_transfers'; ...
         'e  = leg1_transfers.entries(1);  r = leg1_transfers.regime;'; ...
         'pr = nd_propulsion(r.thrustN, r.ispS, r.m0kg, r.lStar, r.tStar);'; ...
         'phys = struct(''Tnd'', pr.Tnd, ''cnd'', pr.cnd, ''muStar'', r.muStar, ''lStar'', r.lStar, ''tStar'', r.tStar, ''m0kg'', r.m0kg);'; ...
         'F  = fly_transfer(e.z8, e.rv0, e.rvf, phys);   % lands within e.flyKm of e.rvf'};
E = struct('entries', entries, 'regime', rg, 'package', package, 'provenance', provenance, ...
           'note', note, 'usage', {usage});
E.provenance.exportWallS = toc(tExport);

%% Atomic write:
if ~isfolder(outDir), mkdir(outDir); end
tmpFile = fullfile(outDir, sprintf('.%s.partial_%s.mat', outName, char(java.util.UUID.randomUUID())));
cleanTmp = onCleanup(@() delete_if(tmpFile));
leg1_transfers = E;
save(tmpFile, 'leg1_transfers');
movefile(tmpFile, outFile, 'f');
clear cleanTmp

if verbose
    fprintf('  wrote %s (%d entries)\n', outFile, numel(entries));
    fprintf('  %-4s %-4s %-10s %8s %9s %9s %9s %10s %10s %s\n', 'sat', 'slot', 'seed', 'orient', ...
            't_dep d', 't_f d', 'miss km', 'Earth km', 'Moon km', 'second order');
    for ke = 1:numel(entries)
        e = entries(ke);
        fprintf('  %-4d %-4d %-10s %8.3f %9.4f %9.4f %9.3g %10.1f %10.1f %s\n', e.sat, e.slot, ...
                sprintf('(%d,%d,%d)', e.seededFrom), e.orientDeg, e.tDepDays, e.tf_nd*dayNd, e.flyKm, ...
                e.verdicts.earthMinKm, e.verdicts.moonMinKm, e.verdicts.secondStatus);
    end
end
end

% ------------------------------------------------------------------------
function delete_if(f)
% DELETE_IF  Delete a file if it exists.  INPUTS: f char.  OUTPUTS: none.
if isfile(f), delete(f); end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c logical; a, b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
