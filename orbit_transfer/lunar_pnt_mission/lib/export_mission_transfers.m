function E = export_mission_transfers(resultsMat, opts)
%% Purpose:
%
%   CAPTURE THE MISSION'S NEW TRANSFERS AS COSTATE-LIBRARY ENTRIES, in a
%   SIDECAR file. The study solves six new minimum-time transfers (each a
%   library primary polished onto an exact tulip slot state at a pinned
%   arrival epoch). This function reads the study's saved results and
%   writes them as off-grid entries with EXACTLY the fields of the 70 mN
%   library's cat_.alternatives, so any entry can be handed to the
%   deliverable-9 kit as-is: catalog_transfer(cat_, entry). The library of
%   record is NOT modified (an outFile inside the kit folder is refused).
%
%   ENTRY (one per accepted polished transfer, spacecraft order):
%     sD = mod(droPhase0(k) + tDep/tauD, 1)   DRO phase at departure
%     sA = mod(slotPhase0(j) + tArr/tauT, 1)  tulip phase at arrival
%     iD = iA = 0 (off grid), sheet = 1, source names the study, the
%     spacecraft, the slot and the seed cell.
%     FULL-STATUS ENTRY (the library's convention, make_alternative.m:
%     z8 = C.z, junctions = C.Y, flyKm = C.flyKm): z8 = the CERTIFIER's
%     root, junctions [14 x 24] = its multiple-shooting junctions (state +
%     costate at the START of each of K = 24 uniform segments, column 1 =
%     [rv0; 1; z8(1:7)]), flyKm / flyVms = its flight. The kit flies the
%     junctions over h = z8(8)/K and the z8 from [rv0; 1; z8(1:7)]
%     (catalog_transfer.m 171-184), so z8 and junctions must come from the
%     same solve: both are the certifier's. The root the MISSION flew (the
%     polish's) is kept in .mission.polishZ8; the two differ by
%     .mission.certifyDz (< 1e-6 required; measured 2.3e-13 .. 4.8e-11).
%     FALLBACK ENTRY: z8 = the polish's root, junctions = the K = 24 cut of
%     its own single flight from rv0 (seed_from_z8, first 24 columns),
%     flyKm / flyVms = the polish's flight. The cut matches the library's
%     ms junctions in column 1 exactly and elsewhere to 3e-6 relative
%     (measured on primaries (3,5), (17,40), (16,31)).
%   The endpoints are checked: M.depState(sD) and M.arrState(sA) must
%   equal the polish's rv0 and rvf to 1e-9 (export_mission_transfers:
%   endpoint otherwise).
%
%   STATUS. opts.fullStatus = true (default; measured 19-63 s per entry,
%   mean 40 s, on the six 2026-11-10 transfers, unfenced): the
%   library's OWN certifier, certify_root (DRO_tulip/indirect), runs on
%   each (seed_from_z8 seed, rv0, rvf) and the library's own packager,
%   make_alternative, turns the result into the status, status_reason,
%   inferred, conj, conjVerdict, minLamV, minQmt, dimS, h6Margin and
%   liftMargin fields -- the same rule (optimality_status) and the same
%   reason text as every library entry. The certifier re-polishes; its root
%   must reproduce the stored z8 to 1e-6, otherwise its status is not
%   attributed to this z8 and the fallback applies (reason says so).
%   opts.fullStatus = false (or a certifier throw / non-reproduction):
%   optimality_status is applied to the evidence the polish holds, with
%   certify_root's stage numbering: 1 polish (the polish's ms_tfmin
%   reproduction |R| < 1e-8, = certify_root's tolRelax), 2 flight (100 km,
%   10 m/s). Stage 3 needs the pointwise checks on the TIGHT-tolerance
%   flight and stage 4 the tfMin witness's own flight, which the polish
%   does not run, so the highest stage established is 2 and the status is
%   1 ("neither: necessary conditions not established"), with the reason
%   listing every certify_root check that was not run. h6Margin is then the
%   polish's computed (not gated) value; liftMargin is NaN.
%   inferred = false and borderline = false on every entry.
%
%   REFUSALS: a transfer with P.ok false (export_mission_transfers:
%   notAccepted, before anything is computed); results with no P
%   (noTransfers); a results fingerprint whose catalog identity (name,
%   created, MD5 of z8 and tf_nd, regime constants, slot phasing, DRO
%   phases, orbit-handle MD5s) differs from the catalog now on the path
%   (fingerprint); a schedule row inconsistent with its polish (schedule).
%   A polish REUSED by the study from an earlier results file is still a
%   certified transfer: it is exported, and .mission records reused = true
%   and its origin (file, date, git head).
%
%   PATHS: certify_root, make_alternative and optimality_status live in
%   DRO_tulip/indirect. When certify_root does not resolve, that folder is
%   appended with addpath -end (never in front, so nothing the campaign
%   asserts can be shadowed); the campaign's which-assertions are re-run
%   and the certifier's instruments are asserted to resolve in
%   costate_common. The session path found on entry is restored exactly on
%   return (onCleanup), error or not.
%
%   OUTPUT PATH: resolved to a canonical absolute path (relative to pwd,
%   '..' and links resolved, a not-yet-existing tail included) and refused
%   (export_mission_transfers:outFile) when it lies inside the kit folder.
%
%   The file is written atomically (a temporary file in the same folder,
%   then movefile) and holds ONE variable, mission_transfers (= E).
%
%   Self-demo (nargin == 0): the default results file, full status, the
%   table printed.
%
%% Inputs:
%
%  resultsMat               char (optional)         the study's saved results
%                                                   [results/lunar_pnt_mission_
%                                                   study.mat]
%
%  opts                     struct (optional)       .outFile [mission_transfers_
%                                                   70mN.mat beside resultsMat]
%                                                   .fullStatus [true]
%                                                   .pool [[]] a connected
%                                                   parallel.Pool for
%                                                   certify_root's fence ([] =
%                                                   UNFENCED, as the polish:
%                                                   an unattended run can hang
%                                                   in the certifier;
%                                                   provenance.fenced says
%                                                   which)
%                                                   .verbose [true] print table
%                                                   .only [[]] TEST SEAM:
%                                                   spacecraft indices to
%                                                   export ([] = all)
%
%% Outputs:
%
%  E                        struct                  mission_transfers:
%                                                   .entries [1 x n] (the
%                                                   library's alternatives
%                                                   fields) .mission [1 x n]
%                                                   .catalog .provenance
%                                                   .usage {lines}
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  fix round 1: certifier's z8 + junctions + flight for
%            full-status entries, catalog_identity helper, fenced
%            flag, canonical output path, path restored, .only    (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));
if nargin == 0
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    resultsMat = fullfile(here, 'results', 'lunar_pnt_mission_study.mat');
    opts = struct('verbose', true);
end
if nargin < 2, opts = struct(); end
if isempty(resultsMat), resultsMat = fullfile(here, 'results', 'lunar_pnt_mission_study.mat'); end
known = {'outFile', 'fullStatus', 'pool', 'verbose', 'only'};
for f = fieldnames(opts).'
    if ~any(strcmp(f{1}, known))
        error('export_mission_transfers:opts', 'unknown option "%s"', f{1});
    end
end
outFile    = fdef(opts, 'outFile',    fullfile(fileparts(resultsMat), 'mission_transfers_70mN.mat'));
fullStatus = fdef(opts, 'fullStatus', true);
pool       = fdef(opts, 'pool',       []);
verbose    = fdef(opts, 'verbose',    true);
only       = fdef(opts, 'only',       []);
if ~(islogical(fullStatus) && isscalar(fullStatus))
    error('export_mission_transfers:opts', 'opts.fullStatus must be a logical scalar');
end
if ~isempty(pool) && ~(isa(pool, 'parallel.Pool') && isvalid(pool) && pool.Connected)
    error('export_mission_transfers:opts', 'opts.pool must be [] or a connected parallel.Pool');
end
fenced     = ~isempty(pool);
p0         = path;                     % restored exactly on return
restorePath = onCleanup(@() path(p0));
Kseg       = 24;                       % the library's junction count
reproTol   = 1e-6;                     % certifier root vs stored z8
endTol     = 1e-9;                     % sD / sA vs the polish's endpoints [ND]
tExport    = tic;

% ---- read the results; refuse before anything is computed ----
if ~isfile(resultsMat)
    error('export_mission_transfers:results', 'results file not found: %s', resultsMat);
end
w = whos('-file', resultsMat);
need = {'P', 'S', 'run_meta', 'fingerprint'};
miss = need(~ismember(need, {w.name}));
if ~isempty(miss)
    error('export_mission_transfers:results', '%s lacks %s', resultsMat, strjoin(miss, ', '));
end
R = load(resultsMat, need{:});
P = R.P;
if isempty(P)
    error('export_mission_transfers:noTransfers', ...
          '%s holds no accepted transfers (the study stopped before all six were accepted)', resultsMat);
end
for k = 1:numel(P)
    if ~isequal(P(k).ok, true)
        error('export_mission_transfers:notAccepted', ...
              'transfer of spacecraft %d was not accepted (%s); refusing to export', k, P(k).reason);
    end
end
kitDir = fullfile(fileparts(here), 'DRO_tulip', 'deliverables', 'costate_library_70mN_user_deliverable9');
outFile = canon(outFile);
if startsWith(outFile, [canon(kitDir) filesep])
    error('export_mission_transfers:outFile', 'refusing to write inside the library kit folder: %s', outFile);
end
[outDir, outName] = fileparts(outFile);

% ---- the catalog now on the path must be the one the seeds came from ----
M = mission_inputs();
fpNow = catalog_identity(M);
fn = fieldnames(fpNow);
bad = fn(~cellfun(@(f) isfield(R.fingerprint, f) && isequal(R.fingerprint.(f), fpNow.(f)), fn));
if ~isempty(bad)
    error('export_mission_transfers:fingerprint', ...
          'the results'' catalog fingerprint differs from the catalog on the path in: %s', strjoin(bad.', ', '));
end
S = R.S;
nE = numel(P);
if isempty(only), only = 1:nE; end
src = struct('origin', repmat({'fresh'}, 1, nE), 'file', '', 'date', '', 'head', '');
if isfield(R.run_meta, 'polishSource')
    for k = 1:nE
        q = R.run_meta.polishSource(k);
        src(k) = struct('origin', q.origin, 'file', q.file, 'date', q.date, 'head', q.head);
    end
end

% ---- the library's certifier on the path (appended, never in front) ----
indirect = fullfile(fileparts(here), 'DRO_tulip', 'indirect');
if fullStatus
    if isempty(which('certify_root'))
        addpath(indirect, '-end');
        setup_lunar_pnt_paths();               % the campaign's which-list still holds
    end
    cc = fullfile(fileparts(here), 'costate_common');
    want = {'certify_root', indirect; 'make_alternative', indirect; 'optimality_status', indirect; ...
            'ms_tfmin', cc; 'validate_flight', cc; 'pmp_pointwise_checks', cc; ...
            'mintime_hypothesis_gates', cc; 'lift_margin', cc; 'conj_spectrum', cc; ...
            'run_capped', cc; 'current_pool', cc; 'seed_from_z8', cc};
    for kw = 1:size(want, 1)
        got = fileparts(which(want{kw, 1}));
        if ~strcmp(got, want{kw, 2})
            error('export_mission_transfers:shadow', '%s resolves to "%s", expected "%s"', ...
                  want{kw, 1}, got, want{kw, 2});
        end
    end
elseif isempty(which('optimality_status'))
    addpath(indirect, '-end');
    setup_lunar_pnt_paths();
end
B = struct('Tnd', M.Tnd, 'cnd', M.cnd, 'mu', M.muStar);
dayNd = M.tStar/86400;

% ---- one entry per transfer ----
entries = cell(1, numel(only));
mission = cell(1, numel(only));
for ko = 1:numel(only)
    k = only(ko);
    Pk = P(k);
    jSlot = S.assign(k);
    iD = S.iD(k);  iA = S.iA(k);
    if ~(abs(Pk.tArr - S.tArr(k)) <= 1e-12 && abs(Pk.tfCell - M.tf(iD, iA)) <= 1e-12)
        error('export_mission_transfers:schedule', ...
              'spacecraft %d: the polish does not belong to schedule row (slot %d, cell (%d,%d))', k, jSlot, iD, iA);
    end
    sD = mod(M.droPhase0(k) + Pk.tDep/M.tauD, 1);
    sA = mod(M.slotPhase0(jSlot) + Pk.tArr/M.tauT, 1);
    eD = max(abs(M.depState(sD) - Pk.rv0));
    eA = max(abs(M.arrState(sA) - Pk.rvf));
    if ~(eD < endTol && eA < endTol)
        error('export_mission_transfers:endpoint', ...
              'spacecraft %d: depState(sD) / arrState(sA) differ from rv0 / rvf by %.2e / %.2e ND', k, eD, eA);
    end
    z8   = Pk.z8(:);
    seed = seed_from_z8(z8, Pk.rv0, Kseg, M.Tnd, M.cnd, M.muStar);
    J    = seed.Y(:, 1:Kseg);
    dateTxt = src(k).date;
    if isempty(dateTxt), dateTxt = R.run_meta.date; end
    source = sprintf('lunar_pnt_mission %s: sc %d -> slot %d, polished from primary (%d,%d)', ...
                     dateTxt(1:min(10, end)), k, jSlot, iD, iA);
    if strcmp(src(k).origin, 'reused'), source = [source ' (polish reused)']; end

    % status
    tS = tic;
    C = [];  dzC = NaN;  why = '';
    if fullStatus
        cOpts = struct('pool', pool, 'allowUnfenced', isempty(pool), 'm0kg', M.m0kg, 'sD', sD, 'sA', sA);
        ws = warning('off', 'certify_root:unfenced');
        try
            C = quiet_call(@certify_root, seed, Pk.rv0, Pk.rvf, B, cOpts);
        catch cerr
            why = ['certify_root threw: ' cerr.message];
        end
        warning(ws);
        if isempty(why)
            dzC = max(abs(C.z(:) - z8));
            A = make_alternative(C, sD, sA, 0, 0, source);
            if isempty(A)
                why = sprintf('certify_root result below the library floor (%s)', C.reason);
            elseif ~(dzC < reproTol)
                why = sprintf('certify_root converged to a root %.1e from the stored z8 (limit %g)', dzC, reproTol);
            elseif ~isequal(size(A.junctions), [14 Kseg])
                why = sprintf('certify_root junctions are [%s], not [14 x %d]', num2str(size(A.junctions)), Kseg);
            end
        end
        if isempty(why)
            statusRun = 'certify_root (full)';
        else
            statusRun = ['fallback: ' why];
        end
    else
        statusRun = 'fallback: fullStatus = false';
    end
    if startsWith(statusRun, 'fallback')
        A = evidence_status(Pk, z8, sD, sA, source, statusRun, fdefs(R.fingerprint, 'polishGates', true));
        zE = z8;  JE = J;  kmE = Pk.flight.flyKm;  vmsE = Pk.flight.flyVms;
    else
        % the library's convention: the packager's z8, junctions and flight
        zE = A.z8(:);  JE = A.junctions;  kmE = A.flyKm;  vmsE = A.flyVms;
    end
    statusWallS = toc(tS);

    e = struct('sD', sD, 'sA', sA, 'iD', 0, 'iA', 0, 'z8', zE, 'junctions', JE, ...
               'tf_nd', zE(8), 'status', double(A.status), 'status_reason', char(A.status_reason), ...
               'inferred', logical(A.inferred), 'conj', double(A.conj), 'conjVerdict', char(A.conjVerdict), ...
               'minLamV', double(A.minLamV), 'minQmt', double(A.minQmt), 'dimS', double(A.dimS), ...
               'h6Margin', double(A.h6Margin), 'liftMargin', double(A.liftMargin), ...
               'flyKm', double(kmE), 'flyVms', double(vmsE), 'source', source, ...
               'sheet', 1, 'borderline', false);
    entries{ko} = e;
    mission{ko} = struct('sc', k, 'slot', jSlot, 'seedCell', [iD iA], 'nRev', S.nRev(k), ...
                        'tDepNd', Pk.tDep, 'tArrNd', Pk.tArr, 'tDepDays', Pk.tDep*dayNd, ...
                        'tArrDays', Pk.tArr*dayNd, 'tfNd', Pk.tf, 'tfDays', Pk.tf*dayNd, ...
                        'tfCellNd', Pk.tfCell, 'familyChanged', Pk.familyChanged, ...
                        'propKg', Pk.propKg, 'dvKms', Pk.dvKms, 'flyKm', Pk.flight.flyKm, ...
                        'flyVms', Pk.flight.flyVms, 'gates', Pk.gates, 'thresholds', Pk.thresholds, ...
                        'nEval', Pk.nEval, 'reused', strcmp(src(k).origin, 'reused'), ...
                        'polishOrigin', src(k), 'statusRun', statusRun, 'statusWallS', statusWallS, ...
                        'polishZ8', z8, 'certifyDz', dzC, 'certify', certify_summary(C));
    if verbose
        fprintf('  entry %d: sc %d -> slot %d, status %d (%s, %.1f s)\n', k, k, jSlot, e.status, statusRun, statusWallS);
    end
end
entries = [entries{:}];
mission = [mission{:}];

% ---- catalog, provenance, usage ----
sh = M.sh;
catalog = struct('fingerprint', R.fingerprint, 'file', 'costate_catalog_dro_tulip_70mN.mat', ...
                 'kitFolder', 'orbit_transfer/DRO_tulip/deliverables/costate_library_70mN_user_deliverable9', ...
                 'name', char(string(M.cat.name)), 'created', char(string(M.cat.created)), ...
                 'regime', struct('thrustN', M.thrustN, 'ispS', M.ispS, 'm0kg', M.m0kg, ...
                                  'muStar', M.muStar, 'lStar', M.lStar, 'tStar', M.tStar, ...
                                  'Tnd', M.Tnd, 'cnd', M.cnd), ...
                 'dep_family', sh.dep_family, 'dep_params', sh.dep_params, ...
                 'arr_family', sh.arr_family, 'arr_params', sh.arr_params, ...
                 'status_key', M.cat.status_key);
provenance = struct('run_meta', R.run_meta, 'resultsFile', resultsMat, ...
                    'exportDate', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
                    'exportFunction', 'export_mission_transfers', ...
                    'statusMode', tern(fullStatus, 'full', 'fallback'), ...
                    'fenced', fenced, 'exported', only, ...
                    'exportWallS', NaN, ...
                    'note', provenance_note(R.run_meta, entries, mission, fullStatus, fenced));
usage = {'L  = load(''costate_catalog_dro_tulip_70mN.mat'');  cat_ = L.costate_catalog_dro_tulip_70mN;'; ...
         'load(''mission_transfers_70mN.mat'');               % one variable: mission_transfers'; ...
         'F  = catalog_transfer(cat_, mission_transfers.entries(1));                     % auto: the junctions'; ...
         'F  = catalog_transfer(cat_, mission_transfers.entries(1), struct(''method'', ''z8''));  % one flight from z8'; ...
         'fprintf(''%s: %.3f d, miss %.3g km, status %d\n'', mission_transfers.entries(1).source, F.tfDays, F.flyKm, F.status);'; ...
         '% the kit builds its endpoints with a plain spline of get_family_orbit, this campaign with the'; ...
         '% periodic phase_state; measured 2026-11-10 on the six entries, both methods land on the kit''s'; ...
         '% arrival state: z8 flights 0.07-2.2 m, junction flights (the certifier''s ms junctions) below 0.1 mm.'};
E = struct('entries', entries, 'mission', mission, 'catalog', catalog, ...
           'provenance', provenance, 'usage', {usage});
E.provenance.exportWallS = toc(tExport);

% ---- atomic write: temporary file in the target folder, then movefile ----
if ~isfolder(outDir), mkdir(outDir); end
tmpFile = fullfile(outDir, sprintf('.%s.partial_%s.mat', outName, char(java.util.UUID.randomUUID())));
cleanTmp = onCleanup(@() delete_if(tmpFile));
mission_transfers = E;
save(tmpFile, 'mission_transfers');
movefile(tmpFile, outFile, 'f');
clear cleanTmp

if verbose
    fprintf('  wrote %s (%d entries, %.0f s)\n', outFile, numel(only), E.provenance.exportWallS);
    fprintf('  %-5s %-3s %-4s %-8s %10s %10s %9s %6s %10s  %s\n', 'entry', 'sc', 'slot', 'cell', 'sD', 'sA', ...
            't_f d', 'status', 'certifyDz', 'status run: reason');
    for ke = 1:numel(entries)
        q = mission(ke);  e = entries(ke);
        fprintf('  %-5d %-3d %-4d %-8s %10.6f %10.6f %9.4f %6d %10.2e  %s: %s\n', ke, q.sc, q.slot, ...
                sprintf('(%d,%d)', q.seedCell), e.sD, e.sA, e.tf_nd*dayNd, e.status, q.certifyDz, ...
                q.statusRun(1:min(40, end)), e.status_reason(1:min(50, end)));
    end
end
end

% ------------------------------------------------------------------------
function A = evidence_status(Pk, z8, sD, sA, source, statusRun, polishGates)
% EVIDENCE_STATUS  The status the polish's own evidence supports under the
%   library's rule (optimality_status, stamped branch), with certify_root's
%   stage numbering: stage 1 the polish's ms_tfmin reproduction (|R| < 1e-8,
%   certify_root's tolRelax), stage 2 the flight gate; stages 3 (pointwise
%   checks on the TIGHT flight) and 4 (the witness's own flight) are not
%   run by the polish, so the highest stage is 2 and the status is 1.
%   INPUTS: Pk polish struct; z8 [8 x 1]; sD, sA phases; source char;
%   statusRun char (why the full status did not run); polishGates logical.
%   OUTPUTS: A struct with the alternatives' status fields.
g = Pk.gates;
stage = 0;
if isequal(g.passResid, true), stage = 1; end
if stage == 1 && isequal(g.passFly, true), stage = 2; end
passed = {};
if isequal(g.passResid, true), passed{end+1} = sprintf('ms_tfmin reproduction |R| = %.1e', g.resid); end
if isequal(g.passFly, true),   passed{end+1} = sprintf('flight %.2g km / %.2g m/s', g.flyKm, g.flyVms); end
if isequal(g.passPW, true),    passed{end+1} = 'pointwise checks on the loose flight'; end
if isequal(g.passHyp, true),   passed{end+1} = 'min|lam_v| > 0, min Q_mt > 0, dim S = 1'; end
if isequal(g.passConj, true),  passed{end+1} = 'junction conjugate test PASS (covered)'; end
if isequal(g.passAcc, true),   passed{end+1} = sprintf('tfMin acceptance |dz| = %.1e', g.accDz); end
if ~polishGates, passed{end+1} = 'polish ran with gates off (fly gate only)'; end
r = sprintf(['library status run not attributable (%s); the polish established: %s; NOT run: pointwise ' ...
             'checks on the tight-tolerance flight, the tfMin witness''s own flight, the X2 field cross-check ' ...
             'and H2/H3 lower-bound floors as gates, lift_margin, the H6 margin gate, the dense conjugate ' ...
             'scan (conj_spectrum)'], statusRun, strjoin(passed, '; '));
C = struct('z', z8, 'flyKm', g.flyKm, 'flyVms', g.flyVms, 'ok', false, 'stage', stage, 'reason', r);
[code, reason, inferred] = optimality_status(C);
conj = -1;  verdict = '';
if isequal(g.reproduced, true) && g.conjCovered == 1
    conj = double(g.conjPass);
    tok = regexp(g.conjReason, '^reproduced; (\w+):', 'tokens', 'once');
    if ~isempty(tok), verdict = tok{1}; end
end
hy = g.hyp;
A = struct('status', code, 'status_reason', reason, 'inferred', inferred, 'conj', conj, ...
           'conjVerdict', verdict, 'minLamV', fdefs(hy, 'minLamV', NaN), 'minQmt', fdefs(hy, 'minQmt', NaN), ...
           'dimS', fdefs(hy, 'dimS', NaN), 'h6Margin', fdefs(hy, 'h6Margin', NaN), 'liftMargin', NaN, ...
           'sD', sD, 'sA', sA, 'source', source);
end

% ------------------------------------------------------------------------
function s = certify_summary(C)
% CERTIFY_SUMMARY  The certifier's numbers behind a status (NaN / '' when
%   it did not run), one fixed field set so the .mission rows concatenate.
%   INPUTS: C certify_root result or [].  OUTPUTS: s struct.
num = {'stage', 'normR', 'dz', 'flyKm', 'flyVms', 'flyKmWitness', 'flyVmsWitness', 'Hmax', 'lamMf', ...
       'lamMfUnc', 'adjErr', 'fullGap', 'throttleErr', 'fieldErr', 'adjErrRef', 'nullResidRel', 'Hresid', ...
       'minLamVEstimate', 'minQmtEstimate', 'liftMargin', 'h6Margin', 'conj', 'wallSec'};
s = struct('ran', ~isempty(C), 'ok', false);
for kf = 1:numel(num)
    s.(num{kf}) = NaN;
    if ~isempty(C) && isfield(C, num{kf}) && ~isempty(C.(num{kf})), s.(num{kf}) = double(C.(num{kf})); end
end
s.ok = ~isempty(C) && isequal(C.ok, true);
s.reason = '';  s.note = '';  s.conjVerdict = '';  s.conjDense = [];
if ~isempty(C)
    s.reason = C.reason;  s.note = C.note;  s.conjVerdict = C.conjVerdict;  s.conjDense = C.conjDense;
end
end

% ------------------------------------------------------------------------
function txt = provenance_note(rm, entries, mission, fullStatus, fenced)
% PROVENANCE_NOTE  The one-paragraph disclosure stored with the sidecar.
%   INPUTS: rm run_meta; entries, mission [1 x n]; fullStatus, fenced
%   logical.
%   OUTPUTS: txt char.
n = numel(entries);
gitTxt = sprintf('optimal_control %s%s', fdefs(rm, 'optimal_control', 'unknown'), ...
                 tern(isequal(fdefs(rm, 'optimal_controlDirty', false), true), ' (dirty)', ''));
txt = sprintf(['These %d entries are mission-specific OFF-GRID minimum-time transfers (iD = iA = 0) captured from ' ...
               'the lunar_pnt_mission study run of %s (%s): each is a 70 mN library primary polished onto an exact ' ...
               'tulip slot state at a pinned arrival epoch. They are NOT part of the library of record ' ...
               '(costate_catalog_dro_tulip_70mN.mat is not modified) and the library''s tools do not schedule them. ' ...
               'Each passed the polish''s acceptance gates: the ms_tfmin re-solve reproduces the root with |R| < 1e-8, ' ...
               'the flown miss is < 100 km and < 10 m/s, pmp_pointwise_checks, mintime_hypothesis_gates (min|lam_v| > 0, ' ...
               'min Q_mt > 0, dim S = 1), the junction conjugate test PASS and covered, and pumpkyn tfMin acceptance ' ...
               '|dz| < 1e-6. '], n, rm.date, gitTxt);
full = strcmp({mission.statusRun}, 'certify_root (full)');
if fullStatus && all(full)
    txt = [txt sprintf(['Each was then run through the library''s own certifier, certify_root (DRO_tulip/indirect), ' ...
                        'and carries the status the library''s rule (optimality_status, via make_alternative) assigns: ' ...
                        'statuses [%s]; every check status 4 requires was run (polish, flight, tight-flight pointwise ' ...
                        'checks, tfMin witness and its flight, conjugate test, hypothesis gates with X2 and lift_margin, ' ...
                        'H6, dense conjugate scan). '], strtrim(sprintf('%d ', [entries.status])))];
else
    txt = [txt sprintf(['The library''s certifier was NOT run on %d of them (%s); by the library''s rule ' ...
                        '(optimality_status) an entry whose evidence stops before the tfMin witness stage is status 1, ' ...
                        'and its status_reason lists the certify_root checks not run (tight-flight pointwise checks, the ' ...
                        'witness flight, X2 and H2/H3 floors as gates, lift_margin, the H6 margin gate, the dense ' ...
                        'conjugate scan). Statuses [%s]. '], nnz(~full), ...
                        strjoin(unique({mission(~full).statusRun}), '; '), strtrim(sprintf('%d ', [entries.status])))];
end
if any(full)
    txt = [txt sprintf(['Certified entries store the certifier''s own z8, multiple-shooting junctions and flown miss ' ...
                        '(the library''s convention, make_alternative); the root the mission flew is .mission.polishZ8 ' ...
                        '(max difference %.1e). The certifier ran %s. '], max([mission(full).certifyDz]), ...
                        tern(fenced, 'FENCED (run_capped on a worker pool)', ...
                             'UNFENCED (no worker pool; an unattended run can hang in it, as an unfenced library certification can)'))];
end
if any(~full)
    txt = [txt 'Fallback entries store the polish''s z8, the K = 24 cut of its single flight as junctions and the polish''s flown miss. '];
end
reu = [mission.reused];
if any(reu)
    o = [mission(reu).polishOrigin];
    srcTxt = strjoin(unique(arrayfun(@(q) sprintf('%s (certified %s at %s)', q.file, q.date, q.head), o, ...
                                     'UniformOutput', false), 'stable'), ' and ');
    txt = [txt sprintf('The polishes for spacecraft %s were REUSED by the study from %s, not re-polished in that run.', ...
                       strtrim(sprintf('%d ', [mission(reu).sc])), srcTxt)];
else
    txt = [txt 'Every polish was solved fresh in that run.'];
end
end

% ------------------------------------------------------------------------
function c = canon(p)
% CANON  Canonical absolute path of a file or folder that need not exist:
%   relative to pwd, '.', '..' and symbolic links resolved (java.io.File
%   getCanonicalPath resolves the existing prefix and keeps the tail), no
%   trailing separator.  INPUTS: p char.  OUTPUTS: c char.
p = char(p);
if ~java.io.File(p).isAbsolute(), p = fullfile(pwd, p); end
c = char(java.io.File(p).getCanonicalPath());
end

% ------------------------------------------------------------------------
function out = quiet_call(varargin)
% QUIET_CALL  feval(varargin{:}) with its printed output (fsolve tables)
%   captured by evalc.  INPUTS: varargin {fhandle, args...}.
%   OUTPUTS: out (first output).
[~, out] = evalc('feval(varargin{:})');
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

% ------------------------------------------------------------------------
function v = fdefs(s, f, v0)
% FDEFS  s.(f) if present and non-empty else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
