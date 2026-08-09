function [zEnd, out] = fly_control(z0, tGrid, rhs, opts)
%% Purpose:
%
%   THE flown-control engine (oclib move 2): integrate a solution's
%   reconstructed control through the true continuous dynamics, end to
%   end, and return where the trajectory actually arrives. This is the
%   physically meaningful accuracy check behind orbit_transfer's G1b gate
%   and booster_landing's G2 gate -- "defect is not accuracy" made
%   executable. The DYNAMICS AND CONTROL RECONSTRUCTION belong to the
%   caller (one rhs closure): reconstruction is domain policy (the
%   booster's annulus-feasible direction/magnitude split is provably NOT
%   a plain quadratic -- see hs_quad_ctrl's ADAPTATION note), so the
%   engine owns only the integration structure.
%
%  ASSUMPTIONS / NOTES:
%
% • Mode 'perInterval' integrates each [t_k, t_{k+1}] separately, carrying
%   the state but restarting the integrator at every node (the orbit
%   G1b convention -- restarts affect step selection, so the mode is part
%   of the measured number, not a detail).
% • Mode 'span' integrates [t_1, t_end] in one call (the booster G2
%   convention).
% • The caller's rhs must evaluate its control from GLOBAL time t; inside
%   a per-interval call t stays within that interval, so an
%   interval-lookup closure reproduces per-interval-built controls
%   exactly (the reconstruction is continuous at the shared node).
%
%% Inputs:
%
%  z0                       [nx x 1]                Initial state
%
%  tGrid                    [1 x N+1]               Node times; only the
%                                                   endpoints are used in
%                                                   'span' mode
%
%  rhs                      fhandle                 dz = rhs(t, z): true
%                                                   dynamics with the
%                                                   reconstructed control
%                                                   inside
%
%  opts                     struct (optional)
%   .mode                   char                    'perInterval' (default)
%                                                   | 'span'
%   .solver                 fhandle                 @ode113 (default) |
%                                                   @ode45 | ...
%   .RelTol                 double                  [1e-12]
%   .AbsTol                 double                  [1e-14]
%
%% Outputs:
%
%  zEnd                     [nx x 1]                Terminal state of the
%                                                   flown trajectory
%
%  out                      struct                  .t, .Z: the flown
%                                                   trajectory samples
%                                                   (concatenated across
%                                                   intervals in
%                                                   'perInterval' mode)
%
%% Revision History:
%  M. Casey                                                   (c) 08/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
mode   = d('mode', 'perInterval');
solver = d('solver', @ode113);
oo = odeset('RelTol', d('RelTol', 1e-12), 'AbsTol', d('AbsTol', 1e-14));

t = tGrid(:).';
wantTraj = nargout > 1;
switch lower(mode)
    case 'perinterval'
        z = z0(:);
        tAll = [];  Zall = [];
        for k = 1:numel(t)-1
            [tk, Z] = solver(rhs, [t(k) t(k+1)], z, oo);
            z = Z(end,:).';
            if wantTraj, tAll = [tAll; tk]; Zall = [Zall; Z]; end
        end
        zEnd = z;
    case 'span'
        [tAll, Zall] = solver(rhs, [t(1) t(end)], z0(:), oo);
        zEnd = Zall(end,:).';
    otherwise
        error('oc:fly_control:mode', ...
              'unknown mode ''%s'' (perInterval | span)', mode);
end
if wantTraj, out = struct('t', tAll, 'Z', Zall); end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
