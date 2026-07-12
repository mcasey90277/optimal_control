function [Z, prob, meta] = ifs_seed_mintime(factor, odeOpts)
% IFS_SEED_MINTIME  k=0 (all-burn) IFS anchor from the converged min-time costates.
%
% Builds the continuation anchor for t_f-marching: at t_f = factor*tfMin with
% factor ~ 1, the min-FUEL optimum is all-burn (no coast is beneficial), so the
% IFS structure is k=0 -- a clean 8x8 rendezvous shooting problem. The seed uses
% the min-time indirect solution's initial costates (run_gto_tulip_indirect ->
% zSol(1:7)) plus lamT0 = -Ht(0) (extended-Hamiltonian normalization, H_sigma=0
% at sigma=0), and the total Sundman length tauf is found by integrating the
% hard-burn EOM to t = tf. At factor = 1 the min-time trajectory IS the solution,
% so the seed residual is tiny; ifs_solve2 closes it and continuation steps up.
%
% INPUTS:
%   factor  - t_f as a multiple of tfMin (6.290694 ND); use ~1.00-1.02 [scalar]
%   odeOpts - (optional) integrator options [odeset], default RelTol 1e-12
% OUTPUTS:
%   Z    - anchor unknown vector = lam0 [8x1] (k=0: no nodes, no switch times)
%   prob - IFS problem struct (k=0, uArc=[1], rendezvous terminal, tauParam
%          'direct')
%   meta - struct: k=0, tauSwitch=[], uArc=[1], seedResNorm, lam0
%
% REFERENCES:
%   [1] ms_band/sms_seed_mintime.m (the time-domain analogue).
%   [2] ifs/PLAN_OF_ATTACK.md (Rung 2 t_f-continuation, min-time anchor).

if nargin < 2 || isempty(odeOpts)
    odeOpts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
end

P    = sms_problem(factor, 1);              % endpoints + params + tf=factor*tfMin
ref  = run_gto_tulip_indirect(false);       % converged min-time indirect solution
lam7 = ref.zSol(1:7);
tf   = P.tf;

% lamT0 = -Ht(0) with the hard burn (uArc=1)
y0c = [P.rv0(:); P.m0; 0; lam7(:); 0];
[~, ~, Ht0] = ifs_eom(0, y0c, P.Tmax, P.c, P.muStar, P.pSund, 1);
lam0 = [lam7(:); -Ht0];                     % [lamR; lamV; lamM; lamT] (8x1)

% total Sundman length: integrate the hard-burn EOM to the event t = tf
y0 = [P.rv0(:); P.m0; 0; lam0];
ev = odeset(odeOpts, 'Events', @(s, y) tf_event(y, tf));
sol = ode113(@(s, y) ifs_eom(s, y, P.Tmax, P.c, P.muStar, P.pSund, 1), ...
             [0, 400], y0, ev);
assert(~isempty(sol.xe), 'ifs_seed_mintime: t never reached tf within sigma<=400');
tauf = sol.xe(end);

prob = struct('rv0',P.rv0(:), 'm0',1, 't0',0, 'tau0',0, 'tauf',tauf, ...
    'Tmax',P.Tmax, 'c',P.c, 'muStar',P.muStar, 'pSund',P.pSund, ...
    'k',0, 'uArc',1, 'termMode','rendezvous', 'rvf',P.rvf(:), 'tf',tf, ...
    'odeOpts',odeOpts, 'tauParam','direct');

Z = lam0;
meta = struct('k',0, 'tauSwitch',[], 'uArc',1, 'lam0',lam0);
meta.seedResNorm = norm(ifs_residual(Z, prob));
end

% -------------------------------------------------------------------------
function [val, isterm, dir] = tf_event(y, tf)
% Terminal event: carried physical time state reaches tf (upward crossing).
val = y(8) - tf;  isterm = 1;  dir = 1;
end
