function ok = test_sheet_to_catalog_file()
%% Purpose:
%
%   Tests sheet_to_catalog_file -- the converter from an arrival-phase
%   sheet (build_arrival_sheet) plus its departure ribs into the SHEET FILE
%   layout that build_costate_catalog_family consumes (OK / Z8 / TF / sD /
%   sA / rungs / meta). The 70 mN phase sheet needs no schema change: the
%   min-time catalog schema already keys sheets by sD_frac x sA_frac, so
%   this is a layout conversion, not a new format.
%
%   Checks the grid layout, the placement of spine and rib entries, that
%   uncertified grid points stay OK = false with NaN t_f, that the meta
%   block is complete enough for the packager, and that the packager then
%   builds a catalog which the schema validator passes.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
nD = 12;  nA = 12;  sA0 = 0.0754;  sD0 = 0;

mkC = @(sD, sA, tf) struct('ok', true, 'reason', 'certified', 'z', [(1:7)'*tf; tf], ...
                           'tfDays', tf*4.4335, 'sD', sD, 'sA', sA);
S = struct('sA', mod(sA0 + (0:nA-1)/nA, 1), 'TF', nan(1, nA), 'Z8', nan(8, nA), ...
           'cand', {cell(1, nA)}, 'sA0', sA0, 'nA', nA);
for j = [1 3]                                   % two certified spine points
    C = mkC(sD0, S.sA(j), 4 + j/10);
    S.cand{j} = C;  S.TF(j) = C.tfDays;  S.Z8(:, j) = C.z;
end
ribs = {struct('sA', S.sA(1), 'pts', [mkC(11/12, S.sA(1), 4.5), mkC(10/12, S.sA(1), 4.6)])};

out = fullfile(tempdir, 'test_sheet_file.mat');
Q = sheet_to_catalog_file(S, ribs, out, struct('nD', nD, 'sD0', sD0));

% MATLAB drops a trailing singleton, so the one-rung grid is [nD nA] with
% size(.,3) = 1 -- which is what the packager's [nD,nA,nR] = size(OK) reads
ok = chk(ok, isequal(size(Q.OK, 1, 2, 3), [nD nA 1]) && isequal(size(Q.Z8, 1, 2, 3, 4), [8 nD nA 1]), ...
         sprintf('grid shapes: OK %s, Z8 %s', mat2str(size(Q.OK)), mat2str(size(Q.Z8))));
ok = chk(ok, nnz(Q.OK) == 4, sprintf('four certified entries placed (%d)', nnz(Q.OK)));
ok = chk(ok, Q.OK(1,1) && Q.OK(1,3) && ~Q.OK(1,2), 'spine points land on the sD0 row, gaps stay false');
ok = chk(ok, Q.OK(12,1) && Q.OK(11,1) && abs(Q.TF(12,1) - 4.5) < 1e-12, ...
         sprintf('rib points land on their sD rows: TF(12,1) = %.3f ND', Q.TF(12,1)));
ok = chk(ok, all(isnan(Q.TF(~Q.OK))), 'uncertified grid points carry NaN t_f');
ok = chk(ok, abs(Q.Z8(8,1,3,1) - 4.3) < 1e-12, sprintf('z8 placed by (iD,iA): tf = %.3f ND', Q.Z8(8,1,3,1)));
ok = chk(ok, isscalar(Q.rungs) && abs(Q.rungs - 0.070) < 1e-12, 'single 70 mN rung');
need = {'muStar','lStar','tStar','ispS','m0kg','tauDRO','depFamily','depParams', ...
        'arrFamily','arrParams','NpTulip','pmTulip','periodTulip'};
miss = need(~isfield(Q.meta, need));
ok = chk(ok, isempty(miss), sprintf('meta complete for the packager (missing: %s)', strjoin(miss, ',')));
ok = chk(ok, Q.meta.ispS == 900 && abs(Q.meta.tauDRO - 1) < 1e-12 && Q.meta.NpTulip == 7, ...
         sprintf('operating point recorded: Isp %g s, tau_dep %g, Np %d', Q.meta.ispS, Q.meta.tauDRO, Q.meta.NpTulip));

% the packager must accept it and the schema must validate the catalog
d_ = fullfile(tempdir, 'test_sheet_dir');  if ~isfolder(d_), mkdir(d_); end
copyfile(out, fullfile(d_, 'dro_tulip_70mN_tau1_Np7.mat'));
cat_ = build_costate_catalog_family(d_, fullfile(tempdir, 'test_cat.mat'), struct( ...
    'glob', 'dro_tulip_70mN_*.mat', 'name', 'test_catalog', ...
    'description', 'test', 'provenance', 'test', 'depReconstruction', 'test'));
probs = catalog_schema('validate', cat_);
ok = chk(ok, isempty(probs), sprintf('packaged catalog validates (%s)', strjoin(probs, '; ')));
ok = chk(ok, cat_.n_entries == 4 && abs(cat_.thruster.isp_s - 900) < 1e-12, ...
         sprintf('catalog: %d entries, Isp %g s', cat_.n_entries, cat_.thruster.isp_s));

if ok, fprintf('TEST_SHEET_TO_CATALOG_FILE: ALL PASS\n'); else, fprintf('TEST_SHEET_TO_CATALOG_FILE: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
