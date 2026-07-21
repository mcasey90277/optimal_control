function [dY, aux] = ztl_eom_sun(Y, P, regime)
% ZTL_EOM_SUN  Sundman-regularized augmented PMP dynamics (per Sundman-time
% tau), for the ramp-family min-fuel CR3BP problem.
%
% Reparametrizes the physical-time EOM by dt/dtau = kappa = r1^pSund
% (r1 = distance to Earth), stretching the perigee passage so the costate
% rates (dlambda/dtau ~ r1^(pSund-3)) are gentle instead of the physical-time
% 1/r^3 blow-up. The costates are UNCHANGED (physical-time PMP costates); only
% the integration variable is scaled. Physical time t is carried as a passive
% 15th state (for the terminal t=tf constraint). Returns the UN-normalized
% Sundman-time RHS dY/dtau; the flow applies the tauF (total length) scaling.
%
% INPUTS:
%   Y      - augmented state [15x1]: [r(3); v(3); m; lam_r(3); lam_v(3);
%            lam_m; t]; MAY be complex (CS probing)
%   P      - struct: .muStar .c .Tmax .eps .pSund
%   regime - 'on' | 'medium' | 'off'
%
% OUTPUTS:
%   dY  - dY/dtau [15x1] = [kappa*f(y); kappa]
%   aux - (optional) struct: .S .u .kappa .r1  (diagnostic; real Y)
%
% REFERENCES: Sundman regularization (PSR pSund=1.5); ztl_eom.m; SUN_BUILD.md.

y = Y(1:14);
yDot = ztl_eom(y, P, regime);            % physical RHS (14), CS-safe (nargout=1)

r  = y(1:3);
dd = [r(1) + P.muStar; r(2); r(3)];
r1 = sqrt(sum(dd.^2));                    % distance to Earth (CS-safe)
kappa = r1^P.pSund;

dY = [kappa*yDot; kappa];                 % [dy/dtau (14); dt/dtau (1)]

if nargout > 1
    [~, auxy] = ztl_eom(y, P, regime);    % aux uses real() -- Y real here
    aux = struct('S', auxy.S, 'u', auxy.u, 'kappa', real(kappa), 'r1', real(r1));
end
end
