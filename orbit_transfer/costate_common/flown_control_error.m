function [erNd, evNd] = flown_control_error(o, muStar, Tmax, c)
%% Purpose:
%
%   THE flown-control verifier (migration #4): flies a direct solution's
%   RECONSTRUCTED CONTROL once, end to end, and reports where the
%   spacecraft actually arrives relative to the solution's own terminal
%   state. This is the physically meaningful accuracy number -- 'if you
%   flew this control, where would you arrive?' -- and it is the G1b gate
%   of every campaign's certification. Extracted verbatim from
%   certify_dro_mintime (the math is family-free: CR3BP + thrust +
%   quadratic control reconstruction).
%
%   Why it matters: every solve in the first campaign reported a defect of
%   1e-14 and every one of them was wrong. The defect is a statement about
%   the discretization; THIS is a statement about the trajectory.
%
%  ASSUMPTIONS / NOTES:
%
% • Control reconstruction matches the transcription: quadratic through
%   (U_k, Um_k, U_{k+1}) when midpoint controls are present (Hermite-
%   Simpson), linear otherwise.
% • Integration at RelTol 1e-12 / AbsTol 1e-14 (ode113), restarting the
%   CONTROL parametrization at each interval but never the state.
%
%% Inputs:
%
%  o                        struct                  Direct solution:
%   .X                      [7 x N+1]               States at nodes
%                                                   [r; v; m]
%   .U                      [4 x N+1]               Controls at nodes
%                                                   [dir(3); throttle]
%   .Um                     [4 x N] or []           Midpoint controls
%   .tNodes                 [1 x N+1] or []         Node times; [] uses
%                                                   .s * .tf
%   .s                      [1 x N+1]               Normalized nodes (only
%                                                   read if tNodes empty)
%   .tf                     double                  Final time (ND)
%
%  muStar                   double                  CR3BP mass ratio
%
%  Tmax                     double                  ND thrust acceleration
%                                                   at unit mass fraction
%
%  c                        double                  ND exhaust velocity
%
%% Outputs:
%
%  erNd                     double                  Terminal POSITION error
%                                                   (ND norm)
%
%  evNd                     double                  Terminal VELOCITY error
%                                                   (ND norm)
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isfield(o,'tNodes') && ~isempty(o.tNodes), t = o.tNodes(:).';
else, t = o.s(:).' * o.tf; end
odeo = odeset('RelTol',1e-12,'AbsTol',1e-14);
hasMid = isfield(o,'Um') && ~isempty(o.Um);
z = o.X(:,1);
for k = 1:numel(t)-1
    a = t(k);  b = t(k+1);  h = max(b-a,eps);
    ua = o.U(:,k);  ub = o.U(:,k+1);
    if hasMid
        um = o.Um(:,k);  uOf = @(tt) ctrl_quad(ua, um, ub, (tt-a)/h);
    else
        uOf = @(tt) ua + ((tt-a)/h)*(ub-ua);
    end
    [~,Z] = ode113(@(tt,zz) cr3bp_thrust_rhs(zz, uOf(tt), muStar, Tmax, c), ...
                   [a b], z, odeo);
    z = Z(end,:).';
end
erNd = norm(z(1:3) - o.X(1:3,end));
evNd = norm(z(4:6) - o.X(4:6,end));
end
