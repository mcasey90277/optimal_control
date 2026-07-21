function [z, prob, info] = ztl_ms_seed_sun(lam0, rv0, rvf, tf, P, M)
% ZTL_MS_SEED_SUN  Dynamically-consistent Sundman MS seed by chopping the
% single-shooting trajectory in sigma (continuity ~0 by construction).
%
% Finds tauF (t(tau_f)=tf) by integrating dY/dtau, then integrates the flow
% arc-by-arc in sigma (uniform sigma_k) recording each node state.
%
% INPUTS:  lam0[7], rv0[1x6], rvf[1x6], tf, P (.pSund etc), M
% OUTPUTS: z[15M-7], prob (.rv0 .rvf .tf .sig .M .P), info (.maxContSeed
%          .termErrSeed .tauF)

Y0 = [rv0(:); 1; lam0(:); 0];
oo = odeset('RelTol',getdef(P,'odeRelTol',1e-13),'AbsTol',getdef(P,'odeAbsTol',1e-15), ...
            'Events',@(tau,Y) tf_event(tau,Y,tf));
[~,~,tauE] = ode89(@(tau,Y) ztl_eom_sun(Y,P,'medium'), [0 1e7], Y0, oo);
assert(~isempty(tauE), 'ztl_ms_seed_sun: t never reached tf');
tauF = tauE(1);

sN = linspace(0, 1, M+1);
prob = struct('rv0',rv0(:).','rvf',rvf(:).','tf',tf,'sig',sN,'M',M,'P',P);

Y = cell(1,M);  Y{1} = Y0;
for k = 1:M-1
    o = ztl_flow_sun(Y{k}, tauF, [sN(k) sN(k+1)], P, false);
    Y{k+1} = o.Yf;
end

z = zeros(15*M-7, 1);
z(1:7) = lam0(:);
for k = 2:M
    z(7 + 15*(k-2) + (1:15)) = Y{k};
end
z(end) = tauF;

[~, ~, ri] = ztl_ms_residual_sun(z, prob, false);
info = struct('maxContSeed',ri.maxCont,'termErrSeed',ri.termErr,'tauF',tauF);
end

% ---------------------------------------------------------------------------
function [v,isterm,dir] = tf_event(~, Y, tf)
v = Y(15) - tf;  isterm = 1;  dir = 1;
end
function v = getdef(s,f,d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
