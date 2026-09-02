function ok = test_catalog_schema_v3()
%% Purpose:
%
%   Tests the schema-v3 (objective/gamma axis) branch of catalog_schema:
%
%     1. COMPAT: all five shipped v1/v2 min-time catalogs still validate
%        clean (byte-identical behavior on old versions).
%     2. A synthetic valid v3 min-fuel catalog validates clean.
%     3. Mutations are rejected: missing .axis3, missing .objective,
%        missing junctions (.Yj) on a minfuel catalog, mf_frac shape
%        mismatch, lam0 column-count mismatch.
%     4. The new 'deltav_from_mf' derivation matches c*ln(1/mf)*lStar/tStar.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All cases passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/02/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
ot = fileparts(fileparts(fileparts(mfilename('fullpath'))));

%% 1. compat: shipped catalogs stay clean:
mats = { fullfile(ot,'DRO_tulip','direct','results','costate_catalog_dro_tulip.mat'); ...
         fullfile(ot,'HALO_tulip','direct','results','costate_catalog_halo_tulip.mat'); ...
         fullfile(ot,'DPO_tulip','direct','results','costate_catalog_dpo_tulip.mat'); ...
         fullfile(ot,'HALO_HALO','direct','results','costate_catalog_halo_halo.mat'); ...
         fullfile(ot,'HALO_HALO','direct','results','costate_catalog_halo_halo_B.mat') };
for k = 1:numel(mats)
    L = load(mats{k});  fn = fieldnames(L);
    pb = catalog_schema('validate', L.(fn{1}));
    [~, nm] = fileparts(mats{k});
    ok = chk(ok, isempty(pb), sprintf('compat: %s validates clean', nm));
end

%% 2. synthetic v3 min-fuel catalog:
c3 = mkV3();
pb = catalog_schema('validate', c3);
ok = chk(ok, isempty(pb), sprintf('synthetic v3 minfuel clean (%s)', strjoin(pb, '; ')));

%% 3. mutations rejected:
muts = { @(c) rmfield(c, 'axis3'),                 'missing axis3'; ...
         @(c) rmfield(c, 'objective'),             'missing objective'; ...
         @(c) rmSheet(c, 'Yj'),                    'minfuel without junctions'; ...
         @(c) resizeSheet(c, 'mf_frac'),           'mf_frac shape mismatch'; ...
         @(c) trimLam0(c),                         'lam0 column mismatch' };
for k = 1:size(muts, 1)
    pb = catalog_schema('validate', muts{k,1}(c3));
    ok = chk(ok, ~isempty(pb), sprintf('mutation rejected: %s', muts{k,2}));
end

%% 4. deltav_from_mf derivation:
dv = catalog_schema('derive', c3, 'deltav_from_mf', struct('mf_frac', 0.9));
ref = c3.thruster.c_nd*log(1/0.9)*c3.constants.lStar_km/c3.constants.tStar_s;
ok = chk(ok, abs(dv - ref) < 1e-14, 'deltav_from_mf matches c*ln(1/mf)');

if ok, fprintf('TEST_CATALOG_SCHEMA_V3: ALL PASS\n');
else,  fprintf('TEST_CATALOG_SCHEMA_V3: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function c = mkV3()
% MKV3  Minimal valid v3 min-fuel catalog (1 sheet, 2x2x2 grid, 3 entries).
% INPUTS: none.  OUTPUTS: c catalog struct.
hs = false(2,2,2);  hs(1,1,1) = true;  hs(2,2,1) = true;  hs(1,2,2) = true;
ei = zeros(2,2,2);  ei(hs) = 1:3;
sh = struct('tauDRO', 1.0, 'tau_dep', 1.0, 'tau_arr', 5.236, 'Np', 7, ...
    'pm', -1, 'dep_family', 'dro', 'dep_params', struct('tau', 1.0), ...
    'arr_family', 'tulip', 'arr_params', struct('Np', 7, 'pm', -1), ...
    'period_tulip_nd', 5.236, 'sD_frac', [0 .5], 'sA_frac', [0 .5], ...
    'has_solution', hs, 'tf_nd', 4.4*hs, 'entry_index', ei, ...
    'tfmin_nd', 3.66*ones(2,2), 'p_floor', 0.002*hs, ...
    'mf_frac', 0.94*hs, 'lam0', rand(7,3), 'Yj', rand(14, 13, 3));
c = struct('name', 'test_v3_minfuel', 'description', 'synthetic', ...
    'schema', 3, 'objective', 'minfuel', ...
    'axis3', struct('name', 'gamma', 'values', [1.1 1.2]), ...
    'constants', struct('muStar', 0.01215, 'lStar_km', 389703.264829278, ...
                        'tStar_s', 382981.289129055), ...
    'thruster', struct('isp_s', 1710, 'm0_kg', 150, 'c_nd', 16.48), ...
    'smoothing', struct('family', 'eps'), ...
    'sheets', sh, 'n_entries', 3, 'env', struct('matlab', 'test'), ...
    'derive', struct('note', 'test'));
end

function c = rmSheet(c, f)
% RMSHEET  Remove a per-sheet field.  INPUTS: c; f.  OUTPUTS: c.
c.sheets = rmfield(c.sheets, f);
end

function c = resizeSheet(c, f)
% RESIZESHEET  Break a per-sheet grid's shape.  INPUTS: c; f.  OUTPUTS: c.
c.sheets(1).(f) = zeros(3, 3, 2);
end

function c = trimLam0(c)
% TRIMLAM0  Drop one lam0 column.  INPUTS: c.  OUTPUTS: c.
c.sheets(1).lam0 = c.sheets(1).lam0(:, 1:2);
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
