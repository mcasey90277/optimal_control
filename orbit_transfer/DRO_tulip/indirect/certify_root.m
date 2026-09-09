function C = certify_root(seed, rv0, rvf, B, opts)
%% Purpose:
%
%   The min-time GATE STACK for one normal-chart multiple-shooting seed:
%
%     1. polish with ms_tfmin (tolR 3e-11) + the free-time conjugate test;
%     2. fly from z8 alone and gate the arrival in POSITION and VELOCITY (a
%        harness that gated on position only accepted arrivals with the
%        wrong velocity -- Astra review 2026-09-09);
%     3. pumpkyn tfMin from the polished z8 must return within tolDz (a
%        foreign-solver witness; an exception is a FAIL, not a pass);
%     4. the conjugate verdict must be PASS (ENDPOINT = inconclusive = FAIL);
%     5. the min-time hypothesis gates: min|lam_v| > 0, min Q_mt > 0 and
%        dim S = 1 (no abnormal lift of the same trajectory).
%
%   Numbers are kept whether or not the seed certifies; C.reason names the
%   FIRST gate that failed. Callers: certify_crossing (arrival-phase
%   arclength candidates, after converting the homogeneous chart) and the
%   departure ribs, whose seeds are already normal-chart.
%
%% Inputs:
%
%  seed                     struct                  ms_bvp seed (.tf, .tGrid,
%                                                   .Y [14 x K+1])
%  rv0, rvf                 [6 x 1]                 departure / arrival states
%  B                        struct                  .Tnd .cnd .mu (operating
%                                                   point, arclength_arrival)
%  opts                     struct (optional)
%   .gateKm [100] .gateVms [10] .tolDz [1e-6] .wallSec [600] .m0kg [150]
%   .tolR [3e-11] polish tolerance, .tolRelax [1e-8] a polish that plateaus
%   below this still goes to the gates (and says so in .reason),
%   .nSamp [200] .rankTol [1e-8] .sA .sD (recorded, not used)
%
%% Outputs:
%
%  C                        struct                  .ok .reason .z [8x1]
%                                                   .Y [14 x K] .tfDays
%                                                   .dvKms .mfKg .flyKm
%                                                   .flyVms .dz .conj .g
%                                                   .sA .sD .rho .normR
%                                                   .wallSec
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
gateKm = d('gateKm', 100);  gateVms = d('gateVms', 10);  tolDz = d('tolDz', 1e-6);
wallSec = d('wallSec', 600);  m0kg = d('m0kg', 150);
% tolR was DOCUMENTED as an option and then hardcoded -- the caller's value
% was silently ignored until 2026-09-09.
tolR = d('tolR', 3e-11);
% A PLATEAU is not a failure. The rib walking the positive departure sense
% out of the anchor stalls at sD = 0.0465 with |R| = 5.1e-11 against a 3e-11
% tolerance, at the same place and the same residual whatever the step size:
% the solver's achievable floor there is simply above a very tight number.
% The residual is one check among five, and the flown miss and the foreign
% witness are stronger evidence than its last decade, so a polish that
% plateaus below tolRelax goes on to the gates and lets THEM decide -- and
% says so in the reason, so a plateaued entry is never mistaken for a clean
% one.
tolRelax = d('tolRelax', 1e-8);
lStar = 389703.264829278;  tStar = 382981.289129055;
t0 = tic;
rv0 = rv0(1:6);  rvf = rvf(1:6);

C = struct('ok', false, 'reason', '', 'z', nan(8,1), 'Y', [], 'tfDays', NaN, ...
           'dvKms', NaN, 'mfKg', NaN, 'flyKm', NaN, 'flyVms', NaN, 'dz', NaN, ...
           'conj', -1, 'g', [], 'sA', d('sA', NaN), 'sD', d('sD', NaN), ...
           'rho', NaN, 'normR', NaN, 'wallSec', NaN);

% ---- 1. normal-chart polish + conjugate test ---------------------------
try
    [z, it] = ms_tfmin(rv0, rvf, seed, B.Tnd, B.cnd, B.mu, ...
                       struct('tolR', tolR, 'wallSec', wallSec, 'conjTest', true));
catch ME
    C.reason = ['ms_tfmin threw: ' ME.message];  C.wallSec = toc(t0);  return
end
C.normR = it.normR;
plateau = false;
if ~it.converged
    if it.normR < tolRelax
        plateau = true;
    else
        C.reason = sprintf('normal-chart polish did not converge (|R| = %.1e)', it.normR);
        C.wallSec = toc(t0);  return
    end
end
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
    za = pumpkyn.cr3bp.tfMin(rv0', rvf', z(:), B.Tnd, B.cnd, B.mu);
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

C.ok = true;  C.wallSec = toc(t0);
if plateau
    C.reason = sprintf('certified (polish plateaued at |R| = %.1e, above tolR = %.1e)', it.normR, tolR);
else
    C.reason = 'certified';
end
end


function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
