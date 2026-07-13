% PROBE_COND_AMP  cond(J_MS) and per-arc amplification ratio under
% AMPLIFICATION-equidistributing node placement vs uniform.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

T = load(fullfile(here, 'results', 'z0_accept2_trace.mat'));
A = load(fullfile(here, 'results', 'p0i_fd_finish.mat'));  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;  tfL = a.tf;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
lam0 = T.lam(:);

fprintf('=== cond(J_MS): amplification vs uniform placement ===\n');
for M = [26 52 104 156]
    [zu, pu] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, M, 'uniform');
    [~, Ju]  = ztl_ms_residual(zu, pu, true);
    [za, pa] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, M, 'amplification');
    [~, Ja]  = ztl_ms_residual(za, pa, true);
    amp = arc_amp(lam0, rv0, pa.tNodes, P);
    fprintf('  M=%3d: cond uniform=%.2e  cond amp=%.2e  | amp-placed ||Phi|| max/med=%.1f\n', ...
            M, cond(Ju), cond(Ja), max(amp)/median(amp));
end

% ---------------------------------------------------------------------------
function amp = arc_amp(lam0, rv0, tN, P)
M = numel(tN) - 1;  Y = [rv0(:); 1; lam0(:)];  amp = zeros(1, M);
for k = 1:M
    o = ztl_flow(Y, [tN(k) tN(k+1)], P, true);  amp(k) = norm(o.PHI);  Y = o.yf;
end
end
