% Certify the 70 mN transfer: harvest costates from the converged direct
% solution, refine by multiple shooting (FREE t_f), then hand z8 to the
% independent single-shooting solver. This is the catalog pipeline's own
% route (direct -> covector harvest -> ms_bvp -> tfMin acceptance).
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect');
if isempty(which('casadi.Opti')), addpath(fullfile(getenv('HOME'),'casadi-3.7.0')); end
S='/private/tmp/claude-501/-Users-msc-Desktop-optimal-control/0a56dbd5-11b8-4df7-b2f2-e22812d68425/scratchpad/';
L=load([S 'direct70.mat']); o=L.o;
tStar=382981.289129055; lStar=389703.264829278;
fprintf('direct 70 mN: tf=%.5f ND (%.3f d), success=%d\n', o.tf, o.tf*tStar/86400, o.success);
fprintf('lamDef size: %s\n', mat2str(size(o.lamDef)));
best=[];
for K=[24 48 96]
    try
        [seed, hd] = harvest_ms_seed(o, K);
    catch ME
        fprintf('K=%d harvest FAILED: %s\n', K, ME.message); continue
    end
    fprintf('K=%2d harvest: sign=%+d voteMargin=%.2f', K, hd.sign, hd.voteMargin);
    if isfield(hd,'lamTf'), fprintf(' lamTf=%.4f', median(hd.lamTf)); end
    fprintf('\n');
    [z,it] = ms_tfmin(L.rv0(1:6), L.rvf(1:6), seed, L.Tnd, L.cnd, L.muStar, ...
                      struct('wallSec',900,'conjTest',true));
    fprintf('   ms: converged=%d normR=%.2e tf=%.5f ND (%.3f d) condJ=%.1e conj=%s\n', ...
        it.converged, it.normR, z(8), z(8)*tStar/86400, it.condJ, it.conj.verdict);
    if it.converged && (isempty(best) || z(8)<best.z(8)), best=struct('z',z,'it',it,'K',K); end
end
if ~isempty(best)
    z=best.z;
    zt = pumpkyn.cr3bp.tfMin(L.rv0(1:6)', L.rvf(1:6)', z(:), L.Tnd, L.cnd, L.muStar);
    dz = norm(zt(:)-z(:));
    [~,Y] = pumpkyn.cr3bp.tfMinProp(z(8), [L.rv0(1:6);1;z(1:7)], L.Tnd, L.cnd, L.muStar);
    flyKm = sqrt(sum((Y(end,1:3)-L.rvf(1:3)').^2))*lStar;
    mf = Y(end,7);
    g = mintime_hypothesis_gates(z, L.rv0(1:6), L.Tnd, L.cnd, L.muStar);
    fprintf('\n===== 70 mN CERTIFIED CANDIDATE (K=%d) =====\n', best.K);
    fprintf('  tf      = %.5f ND = %.3f days\n', z(8), z(8)*tStar/86400);
    fprintf('  dV      = %.4f km/s\n', L.cnd*log(1/mf)*lStar/tStar);
    fprintf('  prop    = %.2f kg of 150 (%.1f%%)\n', 150*(1-mf), 100*(1-mf));
    fprintf('  ms normR= %.2e ;  flown arrival = %.4f km\n', best.it.normR, flyKm);
    fprintf('  tfMin   : |dz| = %.2e  -> %s\n', dz, string(dz<1e-6));
    fprintf('  conj    : %s ;  gates: min|lam_v|=%.3e minQmt=%.3e dimS=%d\n', ...
        best.it.conj.verdict, g.minLamV, g.minQmt, g.dimS);
    save([S 'certified70.mat'],'z','best','dz','flyKm','mf','g');
end
