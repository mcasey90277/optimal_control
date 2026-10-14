function P = entry_thrust_program(z8, junctions, rv0, rvf, phys, opts)
%% Purpose:
%
%   A library entry's TRANSFER AND THRUST PROGRAM, flown segment by segment
%   from its multiple-shooting junction states -- accurate where a single
%   flight from z8 amplifies the costates' error over 17-26 days (FINDINGS
%   97). Reports the junction defects and the arrival miss so a recipient
%   sees how well the stored junctions represent the transfer.
%
%  ASSUMPTIONS / NOTES:
%
% • Junctions sit on the uniform grid linspace(0, t_f, K+1)
%   (flight_to_junctions); column k is the augmented state
%   [r; v; m; lam_r; lam_v; lam_m] at the start of segment k.
% • Minimum-time, all-burn: thrust direction u = -lam_v/|lam_v| (checked in
%   the test against the equations of motion), throttle 1.
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lambda0(7); t_f] (ND)
%  junctions                [14 x K]                junction start states
%  rv0, rvf                 [6 x 1]                 departure / arrival states (ND)
%  phys                     struct                  .Tnd .cnd .mu .lStar .tStar .m0kg
%  opts                     struct (optional)       .nSample [50] per segment
%
%% Outputs:
%
%  P                        struct                  .t .tDays .Y .u .massKg
%                                                   .defectKm .defectVms
%                                                   .startErr .flyKm .flyVms .K
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
nS = 50;  if isfield(opts, 'nSample') && ~isempty(opts.nSample), nS = opts.nSample; end
K = size(junctions, 2);
assert(size(junctions, 1) == 14 && K >= 1 && all(isfinite(junctions(:))), 'entry_thrust_program:junctions', ...
       'junctions must be a finite 14 x K array');
tf = z8(8);  h = tf/K;
y1 = [rv0(:); 1; z8(1:7)];
P = struct('K', K, 'startErr', sqrt(sum((junctions(:, 1) - y1).^2))/sqrt(sum(y1.^2)));
t = [];  Y = [];  dK = zeros(1, max(K - 1, 0));  dV = dK;
for k = 1:K
    [tk, Yk] = pumpkyn.cr3bp.tfMinProp(h, junctions(:, k), phys.Tnd, phys.cnd, phys.mu);
    tq = linspace(0, h, nS).';
    Yq = interp1(tk, Yk, tq, 'pchip');
    Yq(end, :) = Yk(end, :);
    if k < K
        dK(k) = sqrt(sum((Yk(end, 1:3) - junctions(1:3, k+1).').^2))*phys.lStar;
        dV(k) = sqrt(sum((Yk(end, 4:6) - junctions(4:6, k+1).').^2))*phys.lStar/phys.tStar*1000;
    end
    if k > 1, tq = tq(2:end);  Yq = Yq(2:end, :); end
    t = [t; (k - 1)*h + tq];  Y = [Y; Yq];
end
lv = Y(:, 11:13);
P.t = t;  P.tDays = t*phys.tStar/86400;  P.Y = Y;
P.u = -lv ./ sqrt(sum(lv.^2, 2));
P.massKg = Y(:, 7)*phys.m0kg;
P.defectKm = dK;  P.defectVms = dV;
P.flyKm = sqrt(sum((Y(end, 1:3) - rvf(1:3).').^2))*phys.lStar;
P.flyVms = sqrt(sum((Y(end, 4:6) - rvf(4:6).').^2))*phys.lStar/phys.tStar*1000;
end
