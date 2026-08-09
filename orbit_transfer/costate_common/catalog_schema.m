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
    out = 2;

case 'validate'
    cat_ = varargin{1};
    p = {};
    need = {'name','description','constants','thruster','rungs_N', ...
            'sheets','n_entries','derive'};
    for k = 1:numel(need)
        if ~isfield(cat_, need{k})
            p{end+1} = sprintf('missing top-level field .%s', need{k});
        end
    end
    if ~isempty(p), out = p; return, end
    v = 1;
    if isfield(cat_, 'schema'), v = cat_.schema; end
    needS = {'tauDRO','Np','pm','sD_frac','sA_frac','has_solution', ...
             'tf_nd','entry_index','z8'};
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
        if isfield(sh,'z8') && isfield(sh,'has_solution') && ...
           size(sh.z8,2) ~= nnz(sh.has_solution)
            p{end+1} = sprintf('sheet %d: z8 columns (%d) ~= solved entries (%d)', ...
                ks, size(sh.z8,2), nnz(sh.has_solution));
        end
        if v >= 2 && isfield(sh,'tau_dep') && isfield(sh,'tauDRO') && ...
           abs(sh.tau_dep - sh.tauDRO) > 1e-12
            p{end+1} = sprintf('sheet %d: tau_dep ~= tauDRO alias', ks);
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
        otherwise
            error('catalog_schema:derive', 'unknown derivation ''%s''', name);
    end

otherwise
    error('catalog_schema:action', ...
          'unknown action ''%s'' (version | validate | derive)', action);
end
end
