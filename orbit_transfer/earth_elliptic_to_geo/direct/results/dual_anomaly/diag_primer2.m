% DIAG_PRIMER2  Node-centered (adjacent-interval-averaged) primer alignment.
% KKT stationarity for alpha_k involves dtau_{k-1}*lam^{(k-1)} + dtau_k*lam^{(k)}
% (uniform dtau here), so the node costate = average of adjacent interval duals.
cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo');
for f = {'M0_coplanar','M1_3d_fixedLf'}
    S = load(['results/' f{1} '.mat']);  o = S.res.fuel;
    N = size(o.lamDef,2);  ss = o.U(4,:);
    lamV = o.lamDef(4:6,:);
    % node-centered: interior nodes k=2..N average intervals (k-1,k)
    lamN = zeros(3, N+1);
    lamN(:,1) = lamV(:,1);  lamN(:,end) = lamV(:,end);
    lamN(:,2:N) = 0.5*(lamV(:,1:N-1) + lamV(:,2:N));
    burn = ss > 0.5;
    cosb = zeros(1, N+1);
    for k = 1:N+1
        cosb(k) = dot(o.U(1:3,k), -lamN(:,k)) / max(norm(lamN(:,k)),1e-30);
    end
    if mean(cosb(burn)) < 0, cosb = -cosb; end
    ang = real(acosd(max(-1,min(1,cosb(burn)))));
    sa = sort(ang);  n = numel(sa);
    fprintf('%s (node-centered): mean=%.4f med=%.6f p90=%.4f p99=%.2f max=%.2f nBurn=%d\n', ...
        f{1}, mean(ang), median(ang), sa(round(0.9*n)), sa(round(0.99*n)), max(ang), n);
end
