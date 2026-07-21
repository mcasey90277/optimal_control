% PROBE_PERIGEE  Does perigee-concentrated placement flatten the per-arc
% amplification (the linearity indicator that limits the trust radius)?
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
T = load('results/z0_accept2_trace.mat');
A = load('results/p0i_fd_finish.mat');  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();  Tmax = 3*P0.Tmax25;  tfL = a.tf;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
lam0 = T.lam(:);  M = 104;

methods = {'uniform', 'perigee1', 'perigee1.5', 'perigee2', 'perigee3'};
fprintf('=== per-arc ||Phi|| profile vs placement (M=%d) ===\n', M);
for im = 1:numel(methods)
    tN = ztl_ms_nodes(lam0, rv0, tfL, P, M, methods{im});
    Y = [rv0(:); 1; lam0];  amp = zeros(1, M);
    for k = 1:M
        o = ztl_flow(Y, [tN(k) tN(k+1)], P, true);  amp(k) = norm(o.PHI);  Y = o.yf;
    end
    fprintf('  %-11s: max||Phi||=%.2e  med=%.2e  max/med=%.0f  (worst 3 arcs: %s)\n', ...
        methods{im}, max(amp), median(amp), max(amp)/median(amp), ...
        mat2str(round(maxk(amp,3)),3));
end
