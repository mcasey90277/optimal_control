function out = catalog_schema(action, varargin)
%% Purpose:
%
%   THE versioned schema authority for compact costate catalogs -- the
%   normative field list, the validator, and the NAMED FORMULA REGISTRY
%   (retiring the accepted debt of `cat.derive` free-form strings being
%   the only statement of the derivations). One home: packagers stamp,
%   pickers and consumers validate, everyone derives through here.
%
%   VERSIONS:
%     1  (implicit; catalogs without a .schema field) tulip-arrival
%        catalogs: sheets keyed by (tauDRO = departure period, Np).
%        dep_family/dep_params present from the halo catalog on.
%     2  arrival-period axis: sheets additionally carry tau_dep, tau_arr,
%        arr_family, arr_params. tauDRO is kept equal to tau_dep and Np is
%        kept (NaN for non-tulip arrivals) so every v1 consumer keeps
%        working on v2 tulip catalogs. Catalog carries .schema = 2 and
%        .env (environment pinning).
%     3  OBJECTIVE/GAMMA axis (fixed-final-time catalogs, 2026-09-02): one
%        catalog per objective (.objective 'minenergy'|'minfuel'; min-TIME
%        catalogs stay v1/v2, untouched). The per-sheet third array axis is
%        NAMED: top-level .axis3 = struct('name','gamma','values',[...])
%        replaces .rungs_N (tf = gamma * sheet.tfmin_nd). Sheets carry
%        .tfmin_nd [nD x nA], .p_floor and .mf_frac grids (mf is STORED --
%        the all-burn m_final derivation is wrong once coasts exist),
%        .lam0 [7 x n] (no free tf), and -- REQUIRED for minfuel, the
%        identifiability rule -- .Yj [14 x K+1 x n] ms junction states.
%        Top level adds .smoothing (family + notes). Derive registry gains
%        'deltav_from_mf'.
%
%% Inputs:
%
%  action                   char                    'version'  current
%                                                             schema number
%                                                   'validate' check a
%                                                             catalog
%                                                             struct
%                                                   'derive'   named
%                                                             formula
%                                                             evaluation
%
%  varargin                                         validate: (cat_)
%                                                   derive:   (cat_, name,
%                                                             args struct)
%
%% Outputs:
%
%  out                      -                       version: double;
%                                                   validate: cellstr of
%                                                   problems ({} = clean);
%                                                   derive: double
%
%   Derive registry (name -> args -> value):
%     'days'       .t_nd                     -> days
%     'thrust_nd'  .thrustN                  -> ND thrust acceleration
%     'm_final'    .tf_nd .thrustN           -> final mass fraction
%     'deltaV_kms' .tf_nd .thrustN           -> km/s (all-burn min-time)
%     'phase_days' .frac .period_nd          -> days
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

switch lower(action)

case 'version'
    out = 3;

case 'validate'
    cat_ = varargin{1};
    p = {};
    v = 1;
    if isfield(cat_, 'schema'), v = cat_.schema; end
    need = {'name','description','constants','thruster', ...
            'sheets','n_entries','derive'};
    if v >= 3
        need = [need, {'axis3','objective','smoothing'}];
    else
        need = [need, {'rungs_N'}];
    end
    for k = 1:numel(need)
        if ~isfield(cat_, need{k})
            p{end+1} = sprintf('missing top-level field .%s', need{k});
        end
    end
    if ~isempty(p), out = p; return, end
    needS = {'tauDRO','Np','pm','sD_frac','sA_frac','has_solution', ...
             'tf_nd','entry_index'};
    if v >= 3
        needS = [needS, {'lam0','tfmin_nd','p_floor','mf_frac'}];
    else
        needS = [needS, {'z8'}];
    end
    if v >= 2
        needS = [needS, {'tau_dep','tau_arr','dep_family','dep_params', ...
                         'arr_family','arr_params','period_tulip_nd'}];
        if ~isfield(cat_, 'env')
            p{end+1} = 'v2 catalog missing .env (environment pinning)';
        end
        af = arrayfun(@(sh) string(sh.arr_family), cat_.sheets);
        if numel(unique(af)) > 1
            p{end+1} = 'mixed arrival families (one catalog per arrival family)';
        end
    end
    for ks = 1:numel(cat_.sheets)
        sh = cat_.sheets(ks);
        for k = 1:numel(needS)
            if ~isfield(sh, needS{k})
                p{end+1} = sprintf('sheet %d: missing .%s', ks, needS{k});
            end
        end
        if v < 3 && isfield(sh,'z8') && isfield(sh,'has_solution') && ...
           size(sh.z8,2) ~= nnz(sh.has_solution)
            p{end+1} = sprintf('sheet %d: z8 columns (%d) ~= solved entries (%d)', ...
                ks, size(sh.z8,2), nnz(sh.has_solution));
        end
        if v >= 3 && isfield(sh,'has_solution')
            nSol = nnz(sh.has_solution);
            if isfield(sh,'lam0') && size(sh.lam0,2) ~= nSol
                p{end+1} = sprintf('sheet %d: lam0 columns (%d) ~= solved entries (%d)', ...
                    ks, size(sh.lam0,2), nSol);
            end
            for gf = {'p_floor','mf_frac'}
                if isfield(sh, gf{1}) && ...
                   ~isequal(size(sh.(gf{1})), size(sh.has_solution))
                    p{end+1} = sprintf('sheet %d: %s shape ~= has_solution', ks, gf{1});
                end
            end
            if isfield(cat_,'axis3') && isfield(cat_.axis3,'values') && ...
               size(sh.has_solution,3) ~= numel(cat_.axis3.values)
                p{end+1} = sprintf('sheet %d: axis-3 length ~= axis3.values', ks);
            end
            if isfield(cat_,'objective') && strcmpi(cat_.objective,'minfuel')
                % identifiability rule: minfuel entries ship junction states
                if ~isfield(sh,'Yj') || size(sh.Yj,3) ~= nSol || size(sh.Yj,1) ~= 14
                    p{end+1} = sprintf(['sheet %d: minfuel catalog requires ' ...
                        '.Yj [14 x K+1 x nSolved] junction states'], ks);
                end
            end
        end
        % conjugate-point verdicts (OPTIONAL, added 2026-08-23): when a
        % sheet carries conj_pass it must be an int8 grid shaped like
        % has_solution, and the catalog must carry the .conj_test
        % provenance block (conj_catalog_pass writes both or neither):
        if isfield(sh, 'conj_pass') && ~isempty(sh.conj_pass)
            if ~isequal(size(sh.conj_pass), size(sh.has_solution))
                p{end+1} = sprintf('sheet %d: conj_pass shape ~= has_solution', ks);
            end
            if ~isfield(cat_, 'conj_test') && ...
               ~any(contains(p, 'conj_test provenance'))
                p{end+1} = 'sheets carry conj_pass but catalog missing .conj_test provenance';
            end
        end
        if v >= 2 && isfield(sh,'tau_dep') && isfield(sh,'tauDRO') && ...
           abs(sh.tau_dep - sh.tauDRO) > 1e-12
            p{end+1} = sprintf('sheet %d: tau_dep ~= tauDRO alias', ks);
        end
    end
    % Propulsion self-consistency (review 2026-09-05, P0.1): the stored
    % exhaust velocity must be the one implied by the stored Isp, else a
    % consumer recomputing c from isp_s gets a different dV than
    % 'deltav_from_mf' (the shipped min-fuel catalog carried Isp 1710 s
    % beside c_nd for 900 s: dV off 1.90x).
    th = cat_.thruster;  cn = cat_.constants;
    if all(isfield(th, {'isp_s','c_nd'})) && all(isfield(cn, {'lStar_km','tStar_s'}))
        cExp = th.isp_s * 9.80665 * cn.tStar_s / (1000 * cn.lStar_km);
        if abs(th.c_nd - cExp) > 1e-6 * cExp
            p{end+1} = sprintf(['thruster.c_nd (%.6f) inconsistent with isp_s = %g s ' ...
                '(implies %.6f, i.e. Isp %.1f s)'], th.c_nd, th.isp_s, cExp, ...
                th.c_nd * 1000 * cn.lStar_km / (9.80665 * cn.tStar_s));
        end
    end
    nTot = 0;
    for ks = 1:numel(cat_.sheets)
        nTot = nTot + nnz(cat_.sheets(ks).has_solution);
    end
    if nTot ~= cat_.n_entries
        p{end+1} = sprintf('n_entries (%d) ~= sum of solved entries (%d)', ...
            cat_.n_entries, nTot);
    end
    out = p;

case 'derive'
    cat_ = varargin{1};  name = varargin{2};  a = varargin{3};
    tStar = cat_.constants.tStar_s;
    lStar = cat_.constants.lStar_km;
    cnd   = cat_.thruster.c_nd;
    m0    = cat_.thruster.m0_kg;
    switch lower(name)
        case 'days'
            out = a.t_nd * tStar/86400;
        case 'thrust_nd'
            out = (a.thrustN/m0) * tStar^2 / (lStar*1000);
        case 'm_final'
            Tnd = catalog_schema('derive', cat_, 'thrust_nd', a);
            out = 1 - Tnd*a.tf_nd/cnd;
        case 'deltav_kms'
            mf = catalog_schema('derive', cat_, 'm_final', a);
            out = cnd*log(1/mf) * lStar/tStar;
        case 'phase_days'
            out = a.frac * a.period_nd * tStar/86400;
        case 'deltav_from_mf'
            % v3 objectives: mf is STORED (coasts break the all-burn
            % identity); dV follows from the rocket equation alone.
            out = cnd*log(1/a.mf_frac) * lStar/tStar;
        otherwise
            error('catalog_schema:derive', 'unknown derivation ''%s''', name);
    end

otherwise
    error('catalog_schema:action', ...
          'unknown action ''%s'' (version | validate | derive)', action);
end
end
