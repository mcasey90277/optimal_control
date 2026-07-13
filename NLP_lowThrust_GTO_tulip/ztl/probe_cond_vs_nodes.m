% PROBE_COND_VS_NODES  How fast does cond(J_MS) fall with node count, and
% where does the amplification concentrate?
%
% Builds the MS seed at several M from the banked Z0 iterate and reports
% cond(J). Also reports the per-arc STM norm profile at a fixed M to reveal
% whether amplification concentrates at perigee arcs (-> perigee-aware nodes)
% or is uniform (-> just add uniform nodes).

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

Mlist = [26 52 78 104 156 208];
fprintf('=== cond(J_MS) vs node count ===\n');
for M = Mlist
    [z, prob] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, M);
    [~, J] = ztl_ms_residual(z, prob, true);
    fprintf('  M=%3d (%.3f rev/arc): unknowns=%4d  cond(J)=%.3e\n', ...
            M, 13/M, numel(z), cond(J));
end

% per-arc STM norm profile at M=52 (perigee concentration?)
fprintf('\n=== per-arc STM norm profile (M=52) ===\n');
M = 52;  tN = linspace(0, tfL, M+1);
Y = [rv0(:); 1; lam0];  amp = zeros(1, M);
for k = 1:M
    o = ztl_flow(Y, [tN(k) tN(k+1)], P, true);
    amp(k) = norm(o.PHI);
    Y = o.yf;
end
fprintf('  max arc ||Phi|| = %.2e at arc %d/%d;  median = %.2e;  ratio = %.1f\n', ...
        max(amp), find(amp==max(amp),1), M, median(amp), max(amp)/median(amp));
fprintf('  arc ||Phi|| deciles: %s\n', mat2str(round(quantile(amp,0:0.1:1)), 3));
save(fullfile(here, 'results', 'probe_cond_vs_nodes.mat'), 'Mlist', 'amp');
