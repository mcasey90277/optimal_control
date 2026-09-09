function C = certify_crossing(p, sA, B, anc, opts)
%% Purpose:
%
%   The gate stack for ONE arrival-phase candidate before it may enter the
%   sheet. A candidate is a homogeneous-chart multiple-shooting root
%   (arclength_ms crossing, p packed as [lam0(7); Y_2..Y_K; t_f; rho]) at a
%   grid arrival phase sA. It is certified only if ALL of:
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
d = @(f,v) fieldd(opts, f, v);
gateKm = d('gateKm', 100);  gateVms = d('gateVms', 10);  tolDz = d('tolDz', 1e-6);
wallSec = d('wallSec', 600);
lStar = 389703.264829278;  tStar = 382981.289129055;
m0kg = d('m0kg', 150);
t0 = tic;

C = struct('ok', false, 'reason', '', 'z', nan(8,1), 'Y', [], 'tfDays', NaN, ...
           'dvKms', NaN, 'mfKg', NaN, 'flyKm', NaN, 'flyVms', NaN, 'dz', NaN, ...
           'conj', -1, 'g', [], 'sA', sA, 'sD', anc.sD, 'rho', NaN, 'normR', NaN, ...
           'wallSec', NaN);

% ---- unpack the homogeneous chart into a normal-chart seed ---------------
K = anc.K;  ctf = anc.ctf;
rho = p(end);  C.rho = rho;
if ~(rho > 1e-6), C.reason = sprintf('rho = %.1e: abnormal, no normal chart', rho); C.wallSec = toc(t0); return, end
lam0 = p(1:7)/rho;
Yj = reshape(p(8:8+14*(K-1)-1), 14, K-1);
Yj(8:14, :) = Yj(8:14, :)/rho;
rv0 = B.rv0(1:6);  rvf = B.stateA(sA);
Y = [[rv0; 1; lam0], Yj];
seed = struct('tf', p(ctf), 'tGrid', linspace(0, p(ctf), K+1), 'Y', [Y, Y(:,end)]);

% ---- 1. normal-chart polish + conjugate test ---------------------------
try
    [z, it] = ms_tfmin(rv0, rvf(1:6), seed, B.Tnd, B.cnd, B.mu, ...
                       struct('tolR', 3e-11, 'wallSec', wallSec, 'conjTest', true));
catch ME
    C.reason = ['ms_tfmin threw: ' ME.message];  C.wallSec = toc(t0);  return
end
C.normR = it.normR;
if ~it.converged, C.reason = sprintf('normal-chart polish did not converge (|R| = %.1e)', it.normR); C.wallSec = toc(t0); return, end
C.z = z(:);  C.Y = it.Y;  C.tfDays = z(8)*tStar/86400;

% ---- 2. flown arrival, position AND velocity ---------------------------
[~, Yf] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0; 1; z(1:7)], B.Tnd, B.cnd, B.mu);
C.flyKm  = norm(Yf(end,1:3) - rvf(1:3)')*lStar;
C.flyVms = norm(Yf(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
mf = Yf(end,7);
C.mfKg = (1 - mf)*m0kg;  C.dvKms = B.cnd*log(1/mf)*lStar/tStar;
if ~(C.flyKm < gateKm),   C.reason = sprintf('flown position miss %.1f km > %g', C.flyKm, gateKm);  C.wallSec = toc(t0); return, end
if ~(C.flyVms < gateVms), C.reason = sprintf('flown velocity miss %.2f m/s > %g', C.flyVms, gateVms); C.wallSec = toc(t0); return, end

% ---- 3. foreign witness ---------------------------------------------------
try
    za = pumpkyn.cr3bp.tfMin(rv0', rvf(1:6)', z(:), B.Tnd, B.cnd, B.mu);
    C.dz = norm(za(:) - z(:));
catch ME
    C.reason = ['tfMin witness threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~(isfinite(C.dz) && C.dz <= tolDz), C.reason = sprintf('tfMin witness |dz| = %.2e > %g', C.dz, tolDz); C.wallSec = toc(t0); return, end

% ---- 4. conjugate test ----------------------------------------------------
if isfield(it, 'conj') && isfield(it.conj, 'pass'), C.conj = double(it.conj.pass); end
if C.conj ~= 1, C.reason = sprintf('conjugate test verdict %d', C.conj); C.wallSec = toc(t0); return, end

% ---- 5. hypothesis gates --------------------------------------------------
try
    g = mintime_hypothesis_gates(z(:), rv0, B.Tnd, B.cnd, B.mu, ...
                                 struct('nSamp', d('nSamp', 200), 'rankTol', d('rankTol', 1e-8)));
catch ME
    C.reason = ['gates threw: ' ME.message];  C.wallSec = toc(t0);  return
end
C.g = g;
if ~(g.minLamV > 0), C.reason = sprintf('min|lam_v| = %.2e not > 0', g.minLamV); C.wallSec = toc(t0); return, end
if ~(g.minQmt > 0),  C.reason = sprintf('min Q_mt = %.2e not > 0', g.minQmt);   C.wallSec = toc(t0); return, end
if g.dimS ~= 1,      C.reason = sprintf('dim S = %d (abnormal lift)', g.dimS);   C.wallSec = toc(t0); return, end

C.ok = true;  C.reason = 'certified';  C.wallSec = toc(t0);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
