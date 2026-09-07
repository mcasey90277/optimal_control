function S = gates_catalog_pass(catMat, opts)
% GATES_CATALOG_PASS  Run the min-time sufficiency-hypothesis gates
%   (mintime_hypothesis_gates: strong Legendre min|lam_v|, all-burn
%   min Q_mt, abnormal-lift dim S) over every entry of a compact costate
%   catalog and record them -- the catalog-scale form of the audit
%   (doc/mintime_second_order_audit.tex, section 7).
%   Same campaign contract as conj_catalog_pass: endpoints rebuilt from the
%   sheet recipes, one tfMinProp flight + one 7x7 adjoint integration per
%   entry, sidecar progress .mat after EVERY entry, attempt counter before
%   each flight, clean batch-budget exit, resume for free; writeback into
%   the catalog only on an explicit call after a complete census.
% INPUTS:
%   catMat - path to a catalog .mat (single variable, schema v1/v2) [char]
%   opts   - (optional) struct: .logFile [''], .batchSec [inf],
%            .maxEntries [inf], .maxAtt [2], .nSamp [200], .rankTol [1e-8],
%            .sideMat [<catMat minus .mat>_gatesprog.mat], .writeback [false]
% OUTPUTS:
%   S - struct: .done, .nDone, .nTodo, .nH2fail (min|lam_v| <= lamVTol),
%       .nH3fail (min Q_mt <= 0), .nAbnormal (dim S ~= 1), .sideMat
% REFERENCES:
%   [1] costate_common/mintime_hypothesis_gates.m (the instrument)
%   [2] costate_common/conj_catalog_pass.m (the campaign skeleton)
%   [3] doc/mintime_second_order_audit.tex (why these three)

if nargin < 2, opts = struct(); end
batchSec = fieldd(opts, 'batchSec', inf);
maxEnt   = fieldd(opts, 'maxEntries', inf);
maxAtt   = fieldd(opts, 'maxAtt', 2);
nSamp    = fieldd(opts, 'nSamp', 200);
rankTol  = fieldd(opts, 'rankTol', 1e-8);
lamVTol  = fieldd(opts, 'lamVTol', 1e-6);
logFile  = fieldd(opts, 'logFile', '');
[pth, base] = fileparts(catMat);
sideMat  = fieldd(opts, 'sideMat', fullfile(pth, [base '_gatesprog.mat']));
writeback = fieldd(opts, 'writeback', false);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

L = load(catMat);  fn = fieldnames(L);
assert(numel(fn) == 1, 'expected one variable in %s', catMat);
cat_ = L.(fn{1});  nS = numel(cat_.sheets);

%% Sidecar grids:
if exist(sideMat, 'file')
    P = load(sideMat);
    assert(numel(P.MINLV) == nS, 'sidecar %s incompatible (sheet count changed)', sideMat);
else
    P = struct('catMat', catMat, 'created', datestr(now), 'nSamp', nSamp, 'rankTol', rankTol, ...
               'MINLV', {cell(1,nS)}, 'MINQ', {cell(1,nS)}, 'DIMS', {cell(1,nS)}, ...
               'SVR', {cell(1,nS)}, 'NULLR', {cell(1,nS)}, 'HRES', {cell(1,nS)}, 'ATT', {cell(1,nS)});
    for ks = 1:nS
        sz = size(cat_.sheets(ks).has_solution);
        P.MINLV{ks} = nan(sz);  P.MINQ{ks} = nan(sz);  P.DIMS{ks} = -ones(sz, 'int8');
        P.SVR{ks} = nan(sz);    P.NULLR{ks} = nan(sz); P.HRES{ks} = nan(sz);
        P.ATT{ks} = zeros(sz, 'uint8');
    end
    save(sideMat, '-struct', 'P');
end

mu = cat_.constants.muStar;  lStar = cat_.constants.lStar_km;  tStar = cat_.constants.tStar_s;
cnd = cat_.thruster.c_nd;    m0 = cat_.thruster.m0_kg;
ndT = @(TN) (TN/m0)*tStar^2/(lStar*1000);

tAll = tic;  nDone = 0;
for ks = 1:nS
    sh = cat_.sheets(ks);
    todo = sh.has_solution & isnan(P.MINLV{ks}) & (P.ATT{ks} < maxAtt);
    if ~any(todo(:)), continue, end
    % only the DEPARTURE state is needed (the gates fly the stored z8):
    if isfield(sh, 'dep_family'), depFam = sh.dep_family; depPar = sh.dep_params;
    else,                         depFam = 'dro';         depPar = struct('tau', sh.tauDRO); end
    [tD, rvD] = get_family_orbit(depFam, depPar);
    lg('[sheet %d/%d] %s(%g): %d entries to gate', ks, nS, depFam, sh.tauDRO, nnz(todo));
    [nD, nA, nR] = size(sh.has_solution);
    for iD = 1:nD
     for iA = 1:nA
      for kr = 1:nR
        if ~todo(iD,iA,kr), continue, end
        if toc(tAll) > batchSec || nDone >= maxEnt
            save(sideMat, '-struct', 'P');
            lg('[batch] clean exit: %d entries this call', nDone);
            S = census(P, cat_, sideMat, maxAtt, lamVTol); return
        end
        z8 = sh.z8(:, sh.entry_index(iD,iA,kr));
        if ~all(isfinite(z8)) || z8(8) <= 0
            P.ATT{ks}(iD,iA,kr) = maxAtt;  save(sideMat, '-struct', 'P');
            lg('[s%d %d,%d,r%d] SKIP: bad z8', ks, iD, iA, kr);  continue
        end
        P.ATT{ks}(iD,iA,kr) = P.ATT{ks}(iD,iA,kr) + 1;
        save(sideMat, '-struct', 'P');
        rv0 = interp1(tD, rvD, mod(sh.sD_frac(iD),1)*tD(end), 'spline');
        Tnd = ndT(cat_.rungs_N(kr));
        try
            g = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, mu, struct('nSamp', nSamp, 'rankTol', rankTol));
            P.MINLV{ks}(iD,iA,kr) = g.minLamV;  P.MINQ{ks}(iD,iA,kr) = g.minQmt;
            P.DIMS{ks}(iD,iA,kr)  = int8(min(g.dimS, 127));
            P.SVR{ks}(iD,iA,kr)   = g.svRatio;  P.NULLR{ks}(iD,iA,kr) = g.nullResid;
            P.HRES{ks}(iD,iA,kr)  = g.Hresid;
            if g.minLamV <= lamVTol || g.minQmt <= 0 || g.dimS ~= 1 || g.nullResid > 1e-4
                lg('[s%d %d,%d,r%d] FLAG: min|lv|=%.2e minQ=%.2e dimS=%d svr=%.1e nullR=%.1e Hres=%.1e', ...
                   ks, iD, iA, kr, g.minLamV, g.minQmt, g.dimS, g.svRatio, g.nullResid, g.Hresid);
            end
        catch err
            lg('[s%d %d,%d,r%d] ERROR: %s (att %d/%d)', ks, iD, iA, kr, err.message, P.ATT{ks}(iD,iA,kr), maxAtt);
        end
        save(sideMat, '-struct', 'P');
        nDone = nDone + 1;
        if mod(nDone, 50) == 0, lg('[progress] %d entries this call, %.1f s elapsed', nDone, toc(tAll)); end
      end
     end
    end
end

S = census(P, cat_, sideMat, maxAtt, lamVTol);
lg('[census] gated %d / todo %d / H2 fail %d / H3 fail %d / abnormal %d  (done=%d)', ...
   S.nDone, S.nTodo, S.nH2fail, S.nH3fail, S.nAbnormal, S.done);

if writeback
    assert(S.done, 'census incomplete (%d todo) -- not writing back', S.nTodo);
    bak = [catMat '.bak_gates'];
    if ~exist(bak, 'file'), copyfile(catMat, bak); end
    for ks = 1:nS
        cat_.sheets(ks).gate_min_lamv = P.MINLV{ks};
        cat_.sheets(ks).gate_min_qmt  = P.MINQ{ks};
        cat_.sheets(ks).gate_dimS     = P.DIMS{ks};
    end
    cat_.hyp_gates = struct('date', datestr(now, 'yyyy-mm-dd'), ...
        'instrument', 'costate_common/mintime_hypothesis_gates', 'nSamp', nSamp, 'rankTol', rankTol, ...
        'meaning', ['gate_min_lamv: min_t |lam_v| on the flown arc (strong Legendre, need > 0); ' ...
                    'gate_min_qmt: min_t (|lam_v|/m + lam_m/c) (all-burn is the PMP control, need > 0); ' ...
                    'gate_dimS: dimension of the lift space (1 = no abnormal lift). ' ...
                    'With conj_pass = 1 these are the BCT sufficiency hypotheses, sampled -- ' ...
                    'see doc/mintime_second_order_audit.tex section 6']);
    Lout = struct(fn{1}, cat_);  save(catMat, '-struct', 'Lout');
    lg('[writeback] gates stored in %s (backup: %s)', catMat, bak);
end
end

% ------------------------------------------------------------------------
function S = census(P, cat_, sideMat, maxAtt, lamVTol)
% CENSUS  Count gate results from the sidecar.  INPUTS: P; cat_; sideMat;
% maxAtt; lamVTol.  OUTPUTS: S struct.
nDone = 0; nTodo = 0; nH2 = 0; nH3 = 0; nAb = 0;
for ks = 1:numel(cat_.sheets)
    ok = cat_.sheets(ks).has_solution;
    d = ok & ~isnan(P.MINLV{ks});
    nDone = nDone + nnz(d);
    nTodo = nTodo + nnz(ok & isnan(P.MINLV{ks}) & (P.ATT{ks} < maxAtt));
    nH2 = nH2 + nnz(d & P.MINLV{ks} <= lamVTol);
    nH3 = nH3 + nnz(d & P.MINQ{ks} <= 0);
    nAb = nAb + nnz(d & P.DIMS{ks} ~= 1);
end
S = struct('done', nTodo == 0, 'nDone', nDone, 'nTodo', nTodo, ...
           'nH2fail', nH2, 'nH3fail', nH3, 'nAbnormal', nAb, 'sideMat', sideMat);
end

function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end

function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
