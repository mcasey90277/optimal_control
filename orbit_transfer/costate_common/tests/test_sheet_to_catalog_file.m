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

% every certified candidate carries its VERDICTS: the conjugate test and the
% three sufficiency-hypothesis gates. Packaging used to discard them, so the
% catalog reported a conjugate census of 0/0/0 on entries that were all
% conjugate-certified.
mkC = @(sD, sA, tf) struct('ok', true, 'reason', 'certified', 'z', [(1:7)'*tf; tf], ...
                           'tfDays', tf*4.4335, 'sD', sD, 'sA', sA, 'conj', 1, ...
                           'g', struct('minLamV', 3.2 + tf, 'minQmt', 3.8 + tf, 'dimS', 1));
S = struct('sA', mod(sA0 + (0:nA-1)/nA, 1), 'TF', nan(1, nA), 'Z8', nan(8, nA), ...
           'cand', {cell(1, nA)}, 'sA0', sA0, 'nA', nA);
for j = [1 3]                                   % two certified spine points
    C = mkC(sD0, S.sA(j), 4 + j/10);
    S.cand{j} = C;  S.TF(j) = C.tfDays;  S.Z8(:, j) = C.z;
end
ribs = {struct('sA', S.sA(1), 'pts', [mkC(11/12, S.sA(1), 4.5), mkC(10/12, S.sA(1), 4.6)])};

% THE PROBLEM IDENTITY travels with the sheet. Packaging must take the
% engine, the orbits and the departure phase from what was CERTIFIED, not
% from fresh option defaults -- otherwise a catalog can confidently mislabel
% the physics of a real trajectory. (Astra chain review 2026-09-10.)
S.problem = struct('version', 1, 'thrustN', 0.070, 'ispS', 900, 'm0kg', 150, ...
                   'tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1, 'sD', sD0, ...
                   'muStar', 0.012150585609624, 'lStar', 389703.264829278, ...
                   'tStar', 382981.289129055, 'Tnd', 1.234e-4, 'cnd', 6.789);

out = fullfile(tempdir, 'test_sheet_file.mat');
Q = sheet_to_catalog_file(S, ribs, out, struct('nD', nD));

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
ok = chk(ok, isequal(size(Q.CONJ, 1, 2, 3), [nD nA 1]) && all(Q.CONJ(Q.OK) == 1) && all(Q.CONJ(~Q.OK) == -1), ...
         'conjugate verdicts carried, -1 where there is no entry');
ok = chk(ok, abs(Q.MINLV(1,3,1) - (3.2 + 4.3)) < 1e-12 && Q.DIMS(1,3,1) == 1 && ...
             all(isnan(Q.MINQ(~Q.OK))), 'hypothesis gates carried per cell, NaN where empty');
need = {'muStar','lStar','tStar','ispS','m0kg','tauDRO','depFamily','depParams', ...
        'arrFamily','arrParams','NpTulip','pmTulip','periodTulip'};
miss = need(~isfield(Q.meta, need));
ok = chk(ok, isempty(miss), sprintf('meta complete for the packager (missing: %s)', strjoin(miss, ',')));
ok = chk(ok, Q.meta.ispS == 900 && abs(Q.meta.tauDRO - 1) < 1e-12 && Q.meta.NpTulip == 7, ...
         sprintf('operating point recorded: Isp %g s, tau_dep %g, Np %d', Q.meta.ispS, Q.meta.tauDRO, Q.meta.NpTulip));

% (a) the spine goes on the CERTIFIED departure row, not row 1 by default
S2 = S;  S2.problem.sD = 1/12;
Q2 = sheet_to_catalog_file(S2, {}, '', struct('nD', nD));
ok = chk(ok, Q2.OK(2,1) && ~Q2.OK(1,1), 'a sheet certified at sD = 1/12 lands on row 2, not row 1');

% (b) options may ASSERT the identity, never replace it
threw = false;
try, sheet_to_catalog_file(S, {}, '', struct('nD', nD, 'ispS', 1710)); catch, threw = true; end
ok = chk(ok, threw, 'a conflicting engine in opts is refused, not silently applied');

% (c) an entry needs a located CERTIFICATE, not just a finite time
S3 = S;  S3.cand{1} = struct([]);          % summary says solved, no certificate
Q3 = sheet_to_catalog_file(S3, {}, '', struct('nD', nD));
ok = chk(ok, ~Q3.OK(1,1), 'a finite S.TF with no certificate does not export');

% (d) a rib point carrying ok = true but unusable numbers is refused, and
%     must not displace a good entry
badRib = {struct('sA', S.sA(1), 'pts', struct('ok', true, 'z', [(1:7)'; NaN], ...
                 'tfDays', NaN, 'sD', 11/12, 'sA', S.sA(1), 'conj', 1, 'g', []))};
Q4 = sheet_to_catalog_file(S, badRib, '', struct('nD', nD));
ok = chk(ok, ~Q4.OK(12,1), 'a rib point with a NaN final time is refused');

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
s1 = cat_.sheets(1);
ok = chk(ok, isfield(s1, 'conj_pass') && isequal(size(s1.conj_pass), size(s1.has_solution)) && ...
             nnz(s1.conj_pass == 1) == 4, ...
         sprintf('catalog sheet carries conj_pass (%d passes)', nnz(s1.conj_pass == 1)));
ok = chk(ok, isfield(s1, 'gate_dimS') && all(s1.gate_dimS(s1.has_solution) == 1) && ...
             isfield(cat_, 'conj_test') && isfield(cat_, 'hyp_gates'), ...
         'catalog carries the gate grids AND the provenance the schema demands');

if ok, fprintf('TEST_SHEET_TO_CATALOG_FILE: ALL PASS\n'); else, fprintf('TEST_SHEET_TO_CATALOG_FILE: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
