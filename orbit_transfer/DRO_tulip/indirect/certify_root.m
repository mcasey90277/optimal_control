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
%   .pool [gcp('nocreate')] the fence's worker pool -- WITHOUT one the
%   external calls are unfenced and can hang for hours,
%   .capPolishSec [900] .capFlySec [300] .capWitnessSec [300]
%   .capGatesSec [900] hard caps; exceeding one is a NAMED FAILURE,
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
% HARD WALL-CLOCK FENCE on every external call (2026-09-10). `wallSec`
% below is an IN-PROCESS check: it fires only between solver iterations, so
% it cannot bound a single grinding evaluation -- and a CR3BP segment
% propagation whose iterate parks a junction near a primary can grind for
% hours inside ONE evaluation (run_capped's header, after it happened twice
% on the 0.5 N campaign; then again on 2026-09-09, when the sheet assembly
% and a departure rib each burned 12-16 h at 100% CPU with no progress).
% run_capped cancels the worker, so it is the only real fence. A capped
% call FAILS with a named reason: silently accepting a candidate whose
% witness or gates could not be run would put an unverified entry in the
% catalog, which is the one thing this stack exists to prevent.
pool = d('pool', gcp('nocreate'));
capPolish  = d('capPolishSec',  900);
capFly     = d('capFlySec',     300);
capWitness = d('capWitnessSec', 300);
capGates   = d('capGatesSec',   900);
if isempty(pool)
    warning('certify_root:unfenced', ...
        'no parallel pool: external calls run UNFENCED and can hang indefinitely (see capped_pool)');
end
lStar = 389703.264829278;  tStar = 382981.289129055;
t0 = tic;
rv0 = rv0(1:6);  rvf = rvf(1:6);

C = struct('ok', false, 'reason', '', 'z', nan(8,1), 'Y', [], 'tfDays', NaN, ...
           'dvKms', NaN, 'mfKg', NaN, 'flyKm', NaN, 'flyVms', NaN, 'dz', NaN, ...
           'conj', -1, 'g', [], 'sA', d('sA', NaN), 'sD', d('sD', NaN), ...
           'rho', NaN, 'normR', NaN, 'wallSec', NaN);

% ---- 1. normal-chart polish + conjugate test ---------------------------
try
    [okC, z, it] = fenced(pool, capPolish, @ms_tfmin, 2, rv0, rvf, seed, B.Tnd, B.cnd, B.mu, ...
                          struct('tolR', tolR, 'wallSec', wallSec, 'conjTest', true));
catch ME
    C.reason = ['ms_tfmin threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okC
    C.reason = sprintf('polish exceeded its %g s cap (or the worker errored)', capPolish);
    C.wallSec = toc(t0);  return
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
[okF, ~, Yf] = fenced(pool, capFly, @pumpkyn.cr3bp.tfMinProp, 2, z(8), [rv0; 1; z(1:7)], B.Tnd, B.cnd, B.mu);
if ~okF
    C.reason = sprintf('the flight exceeded its %g s cap', capFly);  C.wallSec = toc(t0);  return
end
C.flyKm  = norm(Yf(end,1:3) - rvf(1:3)')*lStar;
C.flyVms = norm(Yf(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
mf = Yf(end,7);
C.mfKg = (1 - mf)*m0kg;  C.dvKms = B.cnd*log(1/mf)*lStar/tStar;
if ~(C.flyKm < gateKm),   C.reason = sprintf('flown position miss %.1f km > %g', C.flyKm, gateKm);  C.wallSec = toc(t0); return, end
if ~(C.flyVms < gateVms), C.reason = sprintf('flown velocity miss %.2f m/s > %g', C.flyVms, gateVms); C.wallSec = toc(t0); return, end

% ---- 3. foreign witness ---------------------------------------------------
try
    [okW, za] = fenced(pool, capWitness, @pumpkyn.cr3bp.tfMin, 1, rv0', rvf', z(:), B.Tnd, B.cnd, B.mu);
catch ME
    C.reason = ['tfMin witness threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okW
    C.reason = sprintf('tfMin witness exceeded its %g s cap', capWitness);  C.wallSec = toc(t0);  return
end
C.dz = norm(za(:) - z(:));
if ~(isfinite(C.dz) && C.dz <= tolDz), C.reason = sprintf('tfMin witness |dz| = %.2e > %g', C.dz, tolDz); C.wallSec = toc(t0); return, end

% ---- 4. conjugate test ----------------------------------------------------
if isfield(it, 'conj') && isfield(it.conj, 'pass'), C.conj = double(it.conj.pass); end
if C.conj ~= 1, C.reason = sprintf('conjugate test verdict %d', C.conj); C.wallSec = toc(t0); return, end

% ---- 5. hypothesis gates --------------------------------------------------
try
    [okG, g] = fenced(pool, capGates, @mintime_hypothesis_gates, 1, z(:), rv0, B.Tnd, B.cnd, B.mu, ...
                      struct('nSamp', d('nSamp', 200), 'rankTol', d('rankTol', 1e-8)));
catch ME
    C.reason = ['gates threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okG
    C.reason = sprintf('hypothesis gates exceeded their %g s cap', capGates);  C.wallSec = toc(t0);  return
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


function varargout = fenced(pool, capSec, fh, nout, varargin)
% FENCED  One external call under a HARD wall-clock cap when a pool is
% available, direct otherwise. Returns [ok, outputs...].
% INPUTS: pool; capSec; fh; nout; varargin.  OUTPUTS: ok; nout outputs.
if isempty(pool)
    varargout = cell(1, nout + 1);
    [varargout{2:nout+1}] = feval(fh, varargin{:});
    varargout{1} = true;
    return
end
varargout = cell(1, nout + 1);
[varargout{1}, varargout{2:nout+1}] = run_capped(pool, fh, nout, capSec, varargin{:});
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
