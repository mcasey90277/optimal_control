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

% Since oclib move 2 the integration structure lives in the shared engine
% oc.fly_control (perInterval mode = the G1b convention: integrator
% restarts at every node); this wrapper owns the CR3BP dynamics and the
% scheme-matched control reconstruction closure.
if isfield(o,'tNodes') && ~isempty(o.tNodes), t = o.tNodes(:).';
else, t = o.s(:).' * o.tf; end
if isempty(which('oc.fly_control'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'oclib'));
end
hasMid = isfield(o,'Um') && ~isempty(o.Um);
    function u = uOf(tt)
    % Global-time control lookup: locate the interval, evaluate its own
    % reconstruction. Continuous at nodes (both neighbors interpolate the
    % node control exactly), so interval choice at a shared node is
    % value-identical to the old per-interval closures.
    k = min(max(1, find(tt >= t, 1, 'last')), numel(t)-1);
    a = t(k);  h = max(t(k+1)-a, eps);
    if hasMid
        u = ctrl_quad(o.U(:,k), o.Um(:,k), o.U(:,k+1), (tt-a)/h);
    else
        u = o.U(:,k) + ((tt-a)/h)*(o.U(:,k+1)-o.U(:,k));
    end
    end
zEnd = oc.fly_control(o.X(:,1), t, ...
    @(tt,zz) cr3bp_thrust_rhs(zz, uOf(tt), muStar, Tmax, c), ...
    struct('mode','perInterval', 'solver',@ode113, ...
           'RelTol',1e-12, 'AbsTol',1e-14));
erNd = norm(zEnd(1:3) - o.X(1:3,end));
evNd = norm(zEnd(4:6) - o.X(4:6,end));
end
