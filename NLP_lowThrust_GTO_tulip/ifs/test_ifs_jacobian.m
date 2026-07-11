function test_ifs_jacobian()
% TEST_IFS_JACOBIAN  CS Jacobian matches a real finite-difference Jacobian.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
p = cr3bp_lt_params(0.025,15,2100);  pSund = 1.5;
odeOpts = odeset('RelTol',1e-13,'AbsTol',1e-15);

% self-consistent k=1 fixedState problem (as in test_ifs_residual)
Y0 = [0.9;0.0;-0.02; 0.05;1.3;-0.2; 1.0;0.0; -0.3;0.1;0.05; -0.7;0.2;0.03; -1.25;-0.4];
tau1 = 0.03;  tauf = 0.06;
[~,Ya] = ode113(@(s,y) ifs_eom(s,y,p.Tmax,p.c,p.muStar,pSund,1),[0 tau1],Y0,odeOpts);
N1 = Ya(end,:).';
[~,Yb] = ode113(@(s,y) ifs_eom(s,y,p.Tmax,p.c,p.muStar,pSund,0),[tau1 tauf],N1,odeOpts);
prob = struct('rv0',Y0(1:6),'m0',Y0(7),'t0',Y0(8),'tau0',0,'Tmax',p.Tmax,'c',p.c, ...
    'muStar',p.muStar,'pSund',pSund,'tauf',tauf,'k',1,'uArc',[1 0], ...
    'termMode','fixedState','termTarget',Yb(end,1:8).','odeOpts',odeOpts);

% perturb off the self-consistent point so the Jacobian is exercised generally
Z = ifs_pack(Y0(9:16), N1, tau1) + 1e-3*randn_fixed(25);
[R0, J] = ifs_residual(Z, prob);
Jfd = zeros(numel(R0));
for c = 1:numel(Z)
    h = 1e-6*max(1, abs(Z(c)));
    Zp = Z;  Zp(c) = Zp(c) + h;
    Rp = ifs_residual(Zp, prob);
    Jfd(:, c) = (Rp - R0)/h;
end
rel = full(max(abs(J - Jfd), [], 'all')) / max(1, full(max(abs(Jfd), [], 'all')));
assert(rel < 1e-4, 'CS Jacobian must match FD, rel err %.2e', rel);
assert(issparse(J), 'J should be sparse');
fprintf('ALL PASS (rel err %.2e, nnz %d/%d)\n', rel, nnz(J), numel(J));
end

function v = randn_fixed(n)
% deterministic pseudo-random perturbation (no rng dependence for reproducibility)
v = sin((1:n).' * 1.7) .* cos((1:n).' * 0.9);
end
