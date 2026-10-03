%% REFILL_24X48_JOB  (2026-10-02)  Re-run the hole filler on the 24 x 48
% build's round 1 after the pool fix (FINDINGS 92: 47 of the filler's 58
% failures were "The parallel pool has shut down"), then re-package, audit
% and sweep exactly as run_phase_torus step 4 does, and refresh final/.
%
% Why a job and not a rerun of reproduce_library_24x48_job: the campaign
% reached its fixed point (status 'done'), so run_phase_torus returns
% 'nothing to do' and never re-enters the filler. This job replays step 4
% with the driver's OWN helpers (commonOpts, roundExtras, assertPackaged via
% its test seam) and the spec saved in reproduce_result.mat, so the inputs
% are the build's, not a copy.
%
% Not replayed: spine-root registration and discovery (steps 4b, 5). A
% certified spine cell (iD = 1) could in principle seed a new anchor; the
% verdict reports how many spine cells the refill certified so that is
% visible rather than assumed.
%
% Launch with nohup, never in the interactive session:
%   nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%       "run('<this file>')" > ~/refill_24x48.out 2>&1 &
% Watch <round_01>/fill_holes_direct.log and the verdict file. It resumes
% (the filler keeps the cells an earlier run certified): run it again.
%
% INPUTS:  none (edit the paths below)
% OUTPUTS: <verdictF> (one line); round_01 catalog re-packaged; final/ refreshed
verdictF = fullfile(getenv('HOME'), 'REFILL_24x48_VERDICT.txt');
R = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results/reproduce_70mN_24x48';

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
try
    V = fullfile(R, 'round_01');  catMat = fullfile(V, 'costate_catalog_dro_tulip_70mN.mat');
    L = load(fullfile(R, 'reproduce_result.mat'));  spec = L.out.spec;
    S = load(fullfile(R, 'torus_state.mat'));  st = S.st;
    H = run_phase_torus('localfunctions');
    in = H.userInputs(spec);
    common = H.commonOpts(in, st.anchors, fullfile(R, 'arcs'), sprintf('arrival_arc_%s_*.mat', in.tag));
    if in.roundDeadlineHours > 0, common.roundDeadlineHours = in.roundDeadlineHours; end
    seedFile = fullfile(R, 'direct_certified.mat');
    if isfile(seedFile), common.seedFiles = {seedFile}; end
    extra = H.roundExtras(R, 1);

    C0 = load(catMat);  c0 = C0.(char(fieldnames(C0)));  n0 = nnz(c0.sheets(1).has_solution(:, :, 1));
    fh = fill_holes_direct(catMat, struct('logFile', fullfile(V, 'fill_holes_direct.log'), 'improveDays', in.improveDays));
    msg = sprintf('REFILL: %d holes, %d tried, %d certified', fh.nHoles, fh.nTried, fh.nCert);
    if fh.nCert > 0
        args = common;
        args.outDir = V;  args.launch = false;
        args.extraRibFiles = unique([extra, {fh.file}], 'stable');
        args.run = struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true);
        o3 = run_costate_library(args);
        H.assertPackaged(o3, 'round_01', 're-package with the refilled cells');
        C1 = load(catMat);  c1 = C1.(char(fieldnames(C1)));  n1 = nnz(c1.sheets(1).has_solution(:, :, 1));
        A = load(fullfile(V, 'audit_70mN.mat'));  a = A.(char(fieldnames(A)));
        fin = fullfile(R, 'final');
        for f = {'costate_catalog_dro_tulip_70mN.mat', 'catalog_receipt.mat', 'second_order_progress_v3.mat', ...
                 'fine_rib_direct_holes.mat', 'phase_torus_70mN.png', 'phase_torus_findings.png'}
            if isfile(fullfile(V, f{1})), copyfile(fullfile(V, f{1}), fin); end
        end
        nSpine = nnz([fh.cells.ok] & [fh.cells.iD] == 1);
        fid = fopen(fullfile(fin, 'README.txt'), 'a');
        fprintf(fid, 'REFILLED %s by batch/refill_24x48_job.m (pool fix): %d -> %d cells; %d spine cell(s) certified (not registered as roots).\n', ...
                char(datetime('now')), n0, n1, nSpine);
        fclose(fid);
        msg = sprintf('%s | cells %d -> %d of %d | spine cells %d | audit %s', msg, n0, n1, numel(c1.sheets(1).has_solution(:, :, 1)), nSpine, auditLine(a));
    end
    msg = ['REFILL FINISHED: ' msg];
catch ME
    msg = sprintf('REFILL FAILED: %s | %s', ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);

function s = auditLine(a)
% AUDITLINE  The audit's ok/bad counts, whatever its field names.
% INPUTS: a (audit struct).  OUTPUTS: s (char).
if isfield(a, 'nOk'), s = sprintf('%d ok / %d bad', a.nOk, a.nBad);
else, s = sprintf('fields: %s', strjoin(fieldnames(a)', ' ')); end
end
