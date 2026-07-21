% PROBE_EQUIL  Does Ruiz two-sided equilibration lower cond(J_MS)?
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
T = load('results/z0_accept2_trace.mat');
A = load('results/p0i_fd_finish.mat');  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();  Tmax = 3*P0.Tmax25;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
for M = [52 104]
    [z, prob] = ztl_ms_seed(T.lam, rv0, rvf, a.tf, P, M);
    [~, J] = ztl_ms_residual(z, prob, true);
    m = size(J,1);  n = size(J,2);  dr = ones(m,1);  dc = ones(n,1);
    for it = 1:5
        B = (spdiags(dr,0,m,m)*J)*spdiags(dc,0,n,n);
        rr = sqrt(max(abs(B),[],2));  rr(rr==0)=1;
        cc = sqrt(max(abs(B),[],1)).'; cc(cc==0)=1;
        dr = dr./rr;  dc = dc./cc;
    end
    Je = (spdiags(dr,0,m,m)*J)*spdiags(dc,0,n,n);
    % column-only (objective-preserving reparametrization dz = Dc w)
    cn = sqrt(sum(J.^2,1)).';  cn(cn==0)=1;
    Jcol = J * spdiags(1./cn, 0, n, n);
    fprintf('M=%3d: cond(J)=%.2e  cond(two-sided)=%.2e  cond(col-only)=%.2e\n', ...
            M, cond(J), cond(full(Je)), cond(Jcol));
end
