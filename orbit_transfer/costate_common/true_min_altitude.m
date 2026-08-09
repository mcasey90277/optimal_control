function amin = true_min_altitude(o, muStar, Tmax, c, lStar, rMoonKm)
%% Purpose:
%
%   Minimum lunar altitude of the PROPAGATED trajectory, not of the nodes.
%   A collocation altitude floor binds at nodes only; this checks it
%   BETWEEN nodes, where periselene actually happens. Extracted verbatim
%   from certify_dro_mintime/local_true_min_alt (migration #4).
%
%  ASSUMPTIONS / NOTES:
%
% • 65 samples per interval, not 9: the floor is active at periselene,
%   where the trajectory turns fastest, and a coarse sample can step
%   straight over the minimum.
%
%% Inputs:
%
%  o                        struct                  Direct solution (same
%                                                   fields as
%                                                   flown_control_error)
%
%  muStar                   double                  CR3BP mass ratio
%
%  Tmax                     double                  ND thrust acceleration
%
%  c                        double                  ND exhaust velocity
%
%  lStar                    double                  Length unit (km)
%
%  rMoonKm                  double                  Lunar radius (km)
%
%% Outputs:
%
%  amin                     double                  Minimum altitude (km)
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

mu = muStar;
if isfield(o,'tNodes') && ~isempty(o.tNodes), t = o.tNodes(:).';
else, t = o.s(:).' * o.tf; end
odeo = odeset('RelTol',1e-11,'AbsTol',1e-13);
amin = inf;
hasMid = isfield(o,'Um') && ~isempty(o.Um);
for k = 1:numel(t)-1
    a = t(k);  b = t(k+1);  h = max(b-a,eps);
    ua = o.U(:,k);  ub = o.U(:,k+1);
    if hasMid
        um = o.Um(:,k);
        uOf = @(tt) ctrl_quad(ua, um, ub, (tt-a)/h);
    else
        uOf = @(tt) ua + ((tt-a)/h)*(ub-ua);
    end
    [~,Z] = ode113(@(tt,z) cr3bp_thrust_rhs(z, uOf(tt), mu, Tmax, c), ...
                   linspace(a,b,65), o.X(:,k), odeo);
    rr = vecnorm(Z(:,1:3) - [1-mu 0 0], 2, 2);
    amin = min(amin, min(rr)*lStar - rMoonKm);
end
end
