function R = mee_residual(o, par, sigma, opts)
% MEE_RESIDUAL  True continuous-time (continuous-longitude) local error of a
%   direct MEE solution -- the shared G1 gate for BOTH MEE campaigns
%   (earth 2-body and CR3BP-GEO), routed through the library engine
%   oc.local_residual. Moved here from earth_elliptic_to_geo/direct/verify
%   2026-08-26 when the CR3BP campaign became the second consumer
%   (migration rule). The lunar third-body term rides in transparently:
%   pass a par carrying .pert (the verify_cr3bp_pmp fingerprint pattern)
%   and lt_mee_rhs's opt-in branch handles it -- this gate never looks.
%
%   For each mesh interval: restart from the transcription's own left node,
%   integrate the TRUE Gauss equations (lt_mee_rhs, d/dL) with the
%   trapezoid-consistent linearly-interpolated control, and measure the miss
%   at the right node. A machine-tight NLP defect says the nodes satisfy the
%   quadrature rule; THIS number says whether they satisfy the dynamics
%   (measured gap on DRO->tulip: 1e7).
%
%   Components are reported SEPARATELY with their own scales -- never norm
%   P (ND length) against t (ND time) and call it kilometres (the exact
%   dimensional bug external review caught in the original DRO gate).
%
% INPUTS:
%   o     - casadi_lt_mee result struct; uses .X [7x(N+1)] = [P;ex;ey;hx;
%           hy;m;t], .U [4x(N+1)] = [beta(3);thr], .dL (total longitude
%           span) [struct]
%   par   - kepler_lt_params struct (lt_mee_rhs parameter block); include
%           .pert (muM/DM/nM/phi0/gain) for a CR3BP lunar-aware row [struct]
%   sigma - node grid in longitude fraction, [(N+1)x1] or [1x(N+1)]
%   opts  - (optional) .relTol [1e-11] .absTol [1e-13] integrator tols
%           (must sit well below the residuals being measured)
%
% OUTPUTS:
%   R - struct:
%       .RP/.RPMax/.RPMed      semi-latus P error [ND; x par.LU for km]
%       .Re/.ReMax             shape-element error norm(ex,ey,hx,hy) [-]
%       .Rm/.RmMax             mass-fraction error [-]
%       .Rt/.RtMax             time error [ND; x par.TU for s]
%       .kWorst .LMid [1xN]    worst interval; interval midpoints [rad]
%       .relTol .absTol .scheme
%
% REFERENCES:
%   [1] oclib/+oc/local_residual.m (the engine; restart-from-node contract)
%   [2] DRO_tulip/direct/certify/dro_residual.m (the reference gate this
%       generalizes; FINDINGS.md carries the 1e7 measurement)
%   [3] core/lt_mee_rhs.m (the dynamics; par.L set per evaluation here)

if nargin < 4, opts = struct(); end
relTol = fieldd(opts, 'relTol', 1e-11);
absTol = fieldd(opts, 'absTol', 1e-13);

% Node longitudes: L_k = pi + sigma_k * DeltaL (casadi_lt_mee convention):
Lg = pi + sigma(:).' * o.dL;
N  = numel(Lg) - 1;

dX = oc.local_residual(o.X, Lg, @(LL, z) rhsL(LL, z, o, par, Lg), ...
        struct('solver', @ode113, 'RelTol', relTol, 'AbsTol', absTol));

RP = abs(dX(1,:));
Re = sqrt(sum(dX(2:5,:).^2, 1));
Rm = abs(dX(6,:));
Rt = abs(dX(7,:));
mix = sqrt(sum(dX.^2, 1));                 % ranking only, dimensionless
[~, kW] = max(mix);

R = struct('RP',RP, 'RPMax',max(RP), 'RPMed',median(RP), ...
           'Re',Re, 'ReMax',max(Re), 'Rm',Rm, 'RmMax',max(Rm), ...
           'Rt',Rt, 'RtMax',max(Rt), 'kWorst',kW, ...
           'LMid',0.5*(Lg(1:N)+Lg(2:N+1)), ...
           'relTol',relTol, 'absTol',absTol, 'scheme','trapezoid');
end

% ------------------------------------------------------------------------
function dz = rhsL(LL, z, o, par, Lg)
% RHSL  d/dL dynamics at longitude LL with the trapezoid-consistent
% linearly-interpolated control. par.L is set per evaluation (lt_mee_rhs
% contract). INPUTS: LL; z [7x1]; o; par; Lg [1xN+1]. OUTPUTS: dz [7x1].
k = find(LL >= Lg, 1, 'last');
k = min(max(k, 1), numel(Lg)-1);
w = (LL - Lg(k)) / max(Lg(k+1)-Lg(k), eps);
u = o.U(:,k) + w*(o.U(:,k+1) - o.U(:,k));
u(1:3) = u(1:3) / max(norm(u(1:3)), eps);      % direction stays unit
u(4)   = min(max(u(4), 0), 1);                 % throttle stays admissible
par.L  = LL;
dz = lt_mee_rhs(z, u, par);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
