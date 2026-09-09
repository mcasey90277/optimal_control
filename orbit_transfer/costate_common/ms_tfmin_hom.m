function [z, info] = ms_tfmin_hom(rv0, rvf, seed, Tmax, c, muStar, opts)
%% Purpose:
%
%   HOMOGENEOUS minimum-time multiple shooting: the same problem as
%   ms_tfmin with the objective multiplier rho FREE and the multipliers
%   confined to the unit sphere,
%
%       H(t_f) = rho + lam'f = 0,        rho^2 + |lam_0|^2 = 1.
%
%   Why: in the normal chart (rho = 1) the DRO -> tulip fast family runs to
%   infinity as thrust falls toward 72 mN -- |lam_0| 46 -> 1449 along the
%   pseudo-arclength arc with every step pure costate motion (FINDINGS 33
%   correction). That is not a fold; it is the chart degenerating as the
%   branch approaches an ABNORMAL configuration. On the sphere the
%   multipliers stay bounded and rho -> 0 becomes a visible, finite event.
%
%   The costate ODE and the control law are positively homogeneous in lam,
%   so (a*lam, a*rho) flies the SAME trajectory: rho enters only the
%   terminal H and the normalisation. pumpkyn's tfMinEoM returns H with
%   the normal "1 +" built in, so H_hom = (H_pumpkyn - 1) + rho and its
%   state Jacobian is unchanged.
%
%   A rho = 1 root (lam_0, t_f) maps to the sphere as
%       a = 1/sqrt(1 + |lam_0|^2),  lam -> a*lam (everywhere), rho = a.
%   Continuing algebraically into rho < 0 does NOT give minimum-time PMP
%   candidates; rho -> 0 is where the family LOSES normality.
%
%% Inputs:
%
%  rv0, rvf                 [6 x 1]                 departure / arrival states
%
%  seed                     struct                  .tf, .tGrid [1 x K+1],
%                                                   .Y [14 x K+1]; .extra =
%                                                   rho seed (if absent the
%                                                   costates are treated as
%                                                   a rho = 1 root and
%                                                   re-normalised)
%
%  Tmax, c, muStar          double                  as ms_tfmin
%
%  opts                     struct (optional)       ms_bvp options
%
%% Outputs:
%
%  z                        [8 x 1]                 [lam_0 (on the sphere);
%                                                   t_f]
%
%  info                     struct                  ms_bvp info plus .rho,
%                                                   .zNormal = [lam_0/rho;
%                                                   t_f] (the rho = 1 chart,
%                                                   valid while rho > 0)
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 7, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here);

if ~isfield(seed, 'extra') || isempty(seed.extra)
    a = 1/sqrt(1 + norm(seed.Y(8:14, 1))^2);
    seed.Y(8:14, :) = a*seed.Y(8:14, :);
    seed.extra = a;
end

prob = struct('ny', 14, 'freeIdx0', 8:14, 'nExtra', 1, ...
    'prop',     @(dt, y0, needSTM) mintime_prop_seg(dt, y0, needSTM, Tmax, c, muStar), ...
    'rhs',      @(y) mintime_rhs_point(y, Tmax, c, muStar), ...
    'terminal', @(y, needJ, x) terminalHom(y, rvf, Tmax, c, muStar, x), ...
    'extraEq',  @(p1, x) sphereEq(p1, x));

seed.Y(1:7, 1) = [rv0(:); 1];
[p, info] = ms_bvp(prob, seed, opts);

z = [p(1:7); p(end-1)];
info.rho = p(end);
if abs(info.rho) > 0, info.zNormal = [p(1:7)/info.rho; p(end-1)]; else, info.zNormal = nan(8,1); end
end

% ------------------------------------------------------------------------
function [g, dgdy, dgdx] = terminalHom(y, rvf, Tmax, c, muStar, x)
% TERMINALHOM  Homogeneous free-tf min-time terminal conditions: r, v
% matched; lam_m(tf) = 0; rho + lam'f = 0.  INPUTS: y [14x1]; rvf; Tmax; c;
% muStar; x = rho.  OUTPUTS: g [8x1]; dgdy [8x14]; dgdx [8x1].
[~, Hn, dHdy] = pumpkyn.cr3bp.tfMinEoM(0, [y; reshape(eye(14), [], 1)], ...
                                       Tmax, c, muStar);
g = [y(1:6) - rvf(:); y(14); (Hn - 1) + x];       % Hn = 1 + lam'f
dgdy = [eye(6), zeros(6, 8);
        zeros(1, 13), 1;
        dHdy];
dgdx = [zeros(7, 1); 1];
end

function [e, dedp1, dedx] = sphereEq(p1, x)
% SPHEREEQ  The normalisation rho^2 + |lam_0|^2 - 1 = 0.
% INPUTS: p1 = lam_0 [7x1]; x = rho.  OUTPUTS: e; dedp1 [1x7]; dedx.
e = x^2 + p1(:)'*p1(:) - 1;
dedp1 = 2*p1(:)';
dedx  = 2*x;
end
