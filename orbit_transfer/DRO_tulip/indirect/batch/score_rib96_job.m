%% SCORE_RIB96_JOB  (2026-10-03)  The departure-axis score: wrap the finished
% 96-phase departure rib (departure_rib96_job) as a one-column catalog
% (rib_as_catalog) and score it with score_interpolator -- 60 seeded random
% departure phases at the rib's arrival phase, blend vs nearest entry. The
% departure twin of the arrival spine's 52 of 60 (FINDINGS 90).
%
% Run AFTER the rib's verdict file says RIB96 FINISHED. Launch with nohup:
%   nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%       "run('<this file>')" > ~/score_rib96.out 2>&1 &
% About 3 h for 60 queries; it resumes. Verdict:
% <outDir>/interp_score/SCORE_VERDICT.txt, and the one line below.
%
% INPUTS:  none (edit the block below)
% OUTPUTS: <outDir>/rib96_col28_catalog.mat, <outDir>/interp_score/, <verdictF>
verdictF = fullfile(getenv('HOME'), 'SCORE_RIB96_VERDICT.txt');
res      = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results';
sheetMat = fullfile(res, 'sheet96_resolution_test', 'arrival_sheet_70mN_nA96.mat');
refMat   = fullfile(res, 'library_70mN_24x48_merged', 'costate_catalog_dro_tulip_70mN.mat');   % problem identity
outDir   = fullfile(res, 'departure_rib96_test');
ribMat   = fullfile(outDir, 'rib96_col28.mat');
nD       = 96;

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'));
cd(ind);
try
    L = load(ribMat);  R = L.R;  assert(isscalar(R), 'score_rib96_job: one rib expected in %s', ribMat);
    S = load(sheetMat);  S = S.S;
    ref = load(refMat);  ref = ref.(char(fieldnames(ref)));
    catMat = fullfile(outDir, sprintf('rib96_col%02d_catalog.mat', R.j));
    c = rib_as_catalog(R, S, ref, struct('nD', nD, 'outMat', catMat));
    fprintf('%s\n', c.note);
    Sc = score_interpolator(catMat, struct('nQuery', 60));
    m = Sc.summary;
    msg = sprintf('SCORE RIB96 FINISHED: %d entries (%d off-lattice dropped, rib stop: %s) | blend usable %d of %d (%.0f%%) | nearest %.0f%%', ...
                  c.n_entries, c.rib.nOffLattice, c.rib.stop, m.nUsable, m.nQuery, 100*m.hitRate, 100*m.baseHitRate);
catch ME
    msg = sprintf('SCORE RIB96 FAILED: %s | %s', ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);
