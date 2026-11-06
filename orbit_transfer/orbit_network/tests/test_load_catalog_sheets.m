function ok = test_load_catalog_sheets()
%% Purpose:
%
%   The sheet table: 56 rows in catalog order, v1 rows mapped onto the v2
%   names, every grid 6 x 6 x 11, unique (dep, arr) keys, arrival periods
%   that match what the catalogs store, and the refusal on a bad root.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

T = load_catalog_sheets([]);
ok = check(ok, numel(T) == 56, sprintf('56 sheets (got %d)', numel(T)));
cats = {T.catalog};
ok = check(ok, isequal(unique(cats, 'stable'), {'dro_tulip', 'halo_tulip', 'dpo_tulip', 'halo_halo'}), 'catalog order DRO, HALO, DPO, HALO_HALO');
ok = check(ok, isequal(cellfun(@(c) nnz(strcmp(cats, c)), {'dro_tulip', 'halo_tulip', 'dpo_tulip', 'halo_halo'}), [16 16 16 8]), 'sheet counts 16/16/16/8');

%% v1 rows mapped onto v2 names
d = T(strcmp(cats, 'dro_tulip'));
ok = check(ok, all(strcmp({d.dep_family}, 'dro')) && all(arrayfun(@(r) isfield(r.dep_params, 'tau'), d)), 'DRO rows: dep_family dro with a tau');
ok = check(ok, all(strcmp({d.arr_family}, 'tulip')) && all(arrayfun(@(r) isfield(r.arr_params, 'Np') && isfield(r.arr_params, 'pm'), d)), 'DRO rows: arrival tulip with Np, pm');
ok = check(ok, all(arrayfun(@(r) r.T_dep == r.dep_params.tau, d)), 'DRO rows: T_dep is the DRO period');
ok = check(ok, all(arrayfun(@(r) abs(r.T_arr - 2*pi*(r.arr_params.Np - 2)/(r.arr_params.Np - 1)) < 1e-9, d)), 'DRO rows: T_arr is the locked tulip period 2 pi (Np-2)/(Np-1)');

%% shapes, rungs, keys
ok = check(ok, all(arrayfun(@(r) isequal(size(r.has_solution), [6 6 11]) && isequal(size(r.tf_nd), [6 6 11]) && isequal(size(r.conj_pass), [6 6 11]), T)), 'every grid is 6 x 6 x 11');
ok = check(ok, all(arrayfun(@(r) islogical(r.has_solution) && islogical(r.conj_pass) && islogical(r.accepted), T)), 'has_solution, conj_pass and accepted are logical');
ok = check(ok, all(arrayfun(@(r) ~any(r.accepted(:) & ~r.has_solution(:)), T)), 'accepted implies has_solution');
r1 = find(T(1).rungs_N == 1);
nAcc = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.accepted, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
nSol = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.has_solution, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
ok = check(ok, all(nAcc <= nSol) && all(nAcc(r1+1:end) < nSol(r1+1:end)) && sum(nSol - nAcc) > 50, sprintf('accepted <= solved at every rung, strictly smaller at 0.75 and 0.5 N (%d refuted or gated cells in all)', sum(nSol - nAcc)));
%% the raw conj_pass holds -1 for "not run": it must not be a pass
S = load(fullfile(fileparts(fileparts(here)), 'HALO_HALO', 'direct', 'results', 'costate_catalog_halo_halo.mat'));
raw = S.costate_catalog_halo_halo.sheets;
nMinus = sum(arrayfun(@(s) nnz(s.conj_pass == -1), raw));
hh = T(strcmp({T.catalog}, 'halo_halo'));
nPassT = sum(arrayfun(@(s) nnz(s.conj_pass), hh));
nPassRaw = sum(arrayfun(@(s) nnz(s.conj_pass == 1), raw));
ok = check(ok, nMinus >= 1 && nPassT == nPassRaw, sprintf('conj_pass == 1 only (%d raw -1 entries are not passes)', nMinus));
ok = check(ok, all(arrayfun(@(r) isequal(r.rungs_N, [15 12 10 7 5 3 2 1.5 1 0.75 0.5]), T)), 'rungs identical and as documented');
keys = arrayfun(@(r) [r.dep_key ' -> ' r.arr_key], T, 'UniformOutput', false);
ok = check(ok, numel(unique(keys)) == 56, 'the 56 (dep, arr) keys are unique');
ok = check(ok, any(strcmp({T.dep_key}, 'halo_L1_tau1.8037')) && any(strcmp({T.arr_key}, 'halo_L2_tau2.2')) && any(strcmp({T.arr_key}, 'tulip_Np12')), 'key spelling: halo_L1_tau1.8037, halo_L2_tau2.2, tulip_Np12');
ok = check(ok, all(arrayfun(@(r) all(isfield(r.constants, {'muStar', 'lStar_km', 'tStar_s'})), T)), 'constants carried');

%% refusal
try
    load_catalog_sheets(fullfile(tempdir, 'no_such_orbit_transfer'));
    ok = check(ok, false, 'bad root refused');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:load_catalog_sheets:missingCatalog'), 'bad root refused by name');
end

fprintf('test_load_catalog_sheets: %s\n', verdict(ok));
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
function s = verdict(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
