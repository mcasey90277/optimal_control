% DIAG_PRIMER  Primer-angle distribution on burn nodes, M0 vs M1.
cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo');
for f = {'M0_coplanar','M1_3d_fixedLf'}
    S = load(['results/' f{1} '.mat']);  o = S.res.fuel;
    N = size(o.lamDef,2);  ss = o.U(4,1:N);  burn = ss > 0.5;
    lamV = o.lamDef(4:6,:);  lVn = sqrt(sum(lamV.^2,1));
    cosb = zeros(1,N);
    for k = 1:N, cosb(k) = dot(o.U(1:3,k), -lamV(:,k))/max(lVn(k),1e-30); end
    if mean(cosb(burn)) < 0, cosb = -cosb; end
    ang = real(acosd(max(-1,min(1,cosb(burn)))));
    sa  = sort(ang);  n = numel(sa);
    fprintf('%s: solver=%.2f mean=%.2f med=%.4f p90=%.2f p99=%.2f max=%.1f nBurn=%d\n', ...
        f{1}, o.primerAlignDeg, mean(ang), median(ang), sa(round(0.9*n)), sa(round(0.99*n)), max(ang), n);
end
