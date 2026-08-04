function P = polish_catalog(sweepMat, outMat)
% POLISH_CATALOG  Wave 2: pumpkyn indirect refinement of every green map cell.
%
% For each cell the sweep proved green, seed pumpkyn.cr3bp.tfMin with the
% harvested costates [lambda0; tf], verify the converged solution by
% re-propagation, and store the REFINED catalog entry. Saves after EVERY
% square -- the raw refined solution is the data product, and nothing here is
% ever the only copy of anything for more than one solve.
%
% INPUTS:
%   sweepMat - path to a sweep result .mat holding S (map, LAM0 catalog, meta
%              with the ORBIT DEFINITION -- the orbits are rebuilt from meta,
%              never assumed)
%   outMat   - output .mat path
%
% OUTPUTS:
%   P - struct: per-cell [nD x nA] .TFI (indirect tf) .RES (terminal residual)
%       .DTF (|tf_direct - tf_indirect|) .OKI (verified) .WALL, refined
%       costates .LAMI [8 x nD x nA], and .meta copied from the sweep so the
%       file stands alone.
%
% REFERENCES:
%   [1] sweep_phasing_direct.m -- the map and catalog this refines.
%   [2] certify/costate_compare.m -- validation of the covector mapping.

M = load(sweepMat);  S = M.S;
ob = S.meta.orbit;
muStar = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
g0   = 9.80665*tStar^2/(1000*lStar);
Tmax = (S.meta.thrustN/S.meta.m0kg)*tStar^2/(lStar*1000);
c    = (S.meta.ispS/tStar)*g0;
[~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);
dep = @(f) interp1(tD, rvD, mod(f,1)*tD(end), 'spline');
arr = @(f) interp1(tT, rvT, mod(f,1)*tT(end), 'spline');

nD = numel(S.sD);  nA = numel(S.sA);
TFI = nan(nD,nA); RES = TFI; DTF = TFI; WALL = TFI; OKI = false(nD,nA);
LAMI = nan(8,nD,nA);
meta = S.meta;  sD = S.sD;  sA = S.sA; %#ok<NASGU>
fprintf('POLISH: %d green cells\n', nnz(S.PASS));
nDone = 0;
for iD = 1:nD
    for iA = 1:nA
        if ~S.PASS(iD,iA), continue, end
        rv0 = dep(S.sD(iD));  rvf = arr(S.sA(iA));
        z0 = S.LAM0(:,iD,iA);
        t0 = tic;
        try
            evalc('z = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z0, Tmax, c, muStar);');
            [~, rv] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6), 1, z(1:7)'], Tmax, c, muStar);
            res = norm(rv(end,1:6) - rvf(1:6));
            ok = res < 1e-7 && all(isfinite(z)) && z(8) > 0.5 && z(8) < 20;
            TFI(iD,iA) = z(8);  RES(iD,iA) = res;
            DTF(iD,iA) = abs(z(8) - S.TF(iD,iA));  WALL(iD,iA) = toc(t0);
            OKI(iD,iA) = ok;
            if ok, LAMI(:,iD,iA) = z; end
            nDone = nDone + 1;
            fprintf('  (%2d,%2d) tfD=%.5f tfI=%.5f |d|=%.2e res=%.1e %s (%.0fs) [%d]\n', ...
                iD, iA, S.TF(iD,iA), z(8), DTF(iD,iA), res, ...
                local_ok(ok), WALL(iD,iA), nDone);
        catch ME
            fprintf('  (%2d,%2d) ERROR %.0fs: %s\n', iD, iA, toc(t0), ME.message);
        end
        save(outMat, 'TFI','RES','DTF','OKI','WALL','LAMI','sD','sA','meta');  % EVERY square
    end
end
P = struct('TFI',TFI,'RES',RES,'DTF',DTF,'OKI',OKI,'WALL',WALL,'LAMI',LAMI, ...
           'sD',S.sD,'sA',S.sA,'meta',S.meta);
fprintf('POLISH DONE: %d/%d verified. median |tfD-tfI| = %.3e ; fastest refined tf = %.7f\n', ...
    nnz(OKI), nnz(S.PASS), median(DTF(OKI)), min(TFI(OKI)));
end

% ---------------------------------------------------------------------------
function s = local_ok(t)
% LOCAL_OK  OK/fail marker. INPUTS: t (logical)  OUTPUTS: s [char]
if t, s = 'OK'; else, s = 'fail'; end
end
