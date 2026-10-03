%% MERGE_24X48_JOB  (2026-10-03)  The refilled 24 x 48 library merged with
% the 24 x 24 record (merge_phase_catalogs: the record's entry wherever it
% is faster or the 24 x 48 has none), then AUDITED in full -- every entry
% re-flown, its tfMin witness and its second-order verdicts recomputed, fail
% closed -- and compared with the record on the shared cells.
%
% Writes into results/library_70mN_24x48_merged/ only; final/ of the build
% and the record are read, never written. Adoption as the library of record
% is a separate decision, made on this job's verdict.
%
% Launch with nohup, never in the interactive session:
%   nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%       "run('<this file>')" > ~/merge_24x48.out 2>&1 &
%
% INPUTS:  none (edit the paths below)
% OUTPUTS: <verdictF> (one line); <outDir>/{costate_catalog_dro_tulip_70mN,
%          merge_info, audit_70mN, comparison_with_record}.mat
verdictF = fullfile(getenv('HOME'), 'MERGE_24x48_VERDICT.txt');
R = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results/';
baseMat  = fullfile(R, 'reproduce_70mN_24x48', 'final', 'costate_catalog_dro_tulip_70mN.mat');
donorMat = fullfile(R, 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
outDir   = fullfile(R, 'library_70mN_24x48_merged');

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
try
    if ~isfolder(outDir), mkdir(outDir); end
    catMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
    [M, info] = merge_phase_catalogs(baseMat, donorMat, struct('donorTag', 'record 24x24 (library_70mN_24x24_final)'));
    L = load(baseMat);  vn = char(fieldnames(L));          % keep the catalog's variable name
    S.(vn) = M;  save(catMat, '-struct', 'S');
    save(fullfile(outDir, 'merge_info.mat'), 'info', 'baseMat', 'donorMat');
    msg = sprintf('MERGE: %d record entries faster, %d holes filled -> %d of %d cells', info.nFaster, info.nFilled, ...
                  M.n_entries, numel(M.sheets.has_solution));

    C = compare_phase_catalogs(catMat, donorMat);
    save(fullfile(outDir, 'comparison_with_record.mat'), 'C');
    msg = sprintf('%s | vs record: %d of %d shared cells agree, missing %d, new slower %d, new faster %d', msg, ...
                  C.nAgree, C.nRef, C.nOnlyRef, C.nNewSlower, C.nNewFaster);

    A = audit_phase_catalog(catMat, struct('out', fullfile(outDir, 'audit_70mN.mat'), 'pool', capped_pool()));
    msg = sprintf('%s | audit %d ok / %d bad', msg, A.nOk, A.nBad);
    msg = ['MERGE FINISHED: ' msg];
catch ME
    msg = sprintf('MERGE FAILED: %s | %s', ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);
