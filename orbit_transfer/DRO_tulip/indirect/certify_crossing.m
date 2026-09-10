function C = certify_crossing(p, sA, B, anc, opts)
%% Purpose:
%
%   The gate stack for ONE arrival-phase candidate before it may enter the
%   sheet. A candidate is a homogeneous-chart multiple-shooting root
%   (arclength_ms crossing, p packed as [lam0(7); Y_2..Y_K; t_f; rho]) at a
%   grid arrival phase sA. This function converts the chart and hands the
%   seed to certify_root, which owns the gates (and is reused by the
%   departure ribs, whose seeds are already normal-chart). It is certified
%   only if ALL of:
%
%     1. normal-chart polish converges (ms_tfmin, rho = 1, tolR 3e-11);
%     2. the flight from z8 alone lands within gateKm AND gateVms of the
%        target -- position and velocity (a sheet harness that gated on
%        position only accepted arrivals with the wrong velocity, Astra
%        review 2026-09-09);
%     3. pumpkyn tfMin, from the polished z8, returns within tolDz of it
%        (foreign-solver witness, MUST converge -- an exception is a FAIL);
%     4. the free-time conjugate test passes (ENDPOINT = inconclusive = FAIL);
%     5. the min-time hypothesis gates hold: min|lam_v| > 0, min Q_mt > 0,
%        dim S = 1 (no abnormal lift).
%
%   Every candidate keeps its numbers whether it passes or not; C.reason
%   names the FIRST gate that failed. The sheet keeps all certified
%   candidates per grid point and takes the minimum t_f; nothing is chosen
%   here.
%
%% Inputs:
%
%  p                        [n x 1]                 homogeneous-chart ms root
%  sA                       double                  its arrival phase (grid)
%  B, anc                   struct                  from arclength_arrival
%                                                   ('setup')
%  opts                     struct (optional)
%   .gateKm [100] .gateVms [10] .tolDz [1e-6] .wallSec [600]
%   .nSamp [200] (gates) .rankTol [1e-8]
%
%% Outputs:
%
%  C                        struct                  .ok .reason .z [8x1]
%                                                   .Y [14 x K] junctions
%                                                   .tfDays .dvKms .mfKg
%                                                   .flyKm .flyVms .dz
%                                                   .conj (1/0/-1) .g (gates
%                                                   struct or []) .sA .sD
%                                                   .rho .normR .wallSec
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end

% ---- unpack the homogeneous chart into a normal-chart seed ---------------
K = anc.K;  ctf = anc.ctf;
% the chart layout is CHECKED, not assumed: p = [lam0(7); Y_2..Y_K; t_f; rho]
nExpect = 7 + 14*(K-1) + 1 + 1;
assert(isnumeric(p) && numel(p) == nExpect && all(isfinite(p)), ...
       'candidate vector is %d long; this chart needs %d (K = %d)', numel(p), nExpect, K);
assert(ctf == nExpect - 1, 'anchor says t_f is at index %d; the chart puts it at %d', ctf, nExpect-1);
rho = p(end);
if ~(rho > 1e-6)
    C = struct('ok', false, 'reason', sprintf('rho = %.1e: abnormal, no normal chart', rho), ...
               'z', nan(8,1), 'Y', [], 'tfDays', NaN, 'dvKms', NaN, 'mfKg', NaN, ...
               'flyKm', NaN, 'flyVms', NaN, 'dz', NaN, 'conj', -1, 'g', [], ...
               'sA', sA, 'sD', anc.sD, 'rho', rho, 'normR', NaN, 'wallSec', 0);
    return
end
lam0 = p(1:7)/rho;
Yj = reshape(p(8:8+14*(K-1)-1), 14, K-1);
Yj(8:14, :) = Yj(8:14, :)/rho;
rv0 = B.rv0(1:6);  rvf = B.stateA(sA);
Y = [[rv0; 1; lam0], Yj];
seed = struct('tf', p(ctf), 'tGrid', linspace(0, p(ctf), K+1), 'Y', [Y, Y(:,end)]);

copts = opts;  copts.sA = sA;  copts.sD = anc.sD;
C = certify_root(seed, rv0, rvf(1:6), B, copts);
C.rho = rho;
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
