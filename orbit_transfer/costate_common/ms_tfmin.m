function [z, info] = ms_tfmin(rv0, rvf, seed, Tmax, c, muStar, opts)
% MS_TFMIN  Multiple-shooting solve of the CR3BP minimum-time PMP problem.
%
% Solves the same problem as pumpkyn.cr3bp.tfMin (terminal conditions
% r(tf)=rf, v(tf)=vf, lambda_m(tf)=0, H(tf)=0 with H = 1 + lambda'f), but
% with the arc split into K short segments whose junction states are extra
% unknowns. Short segments kill the ~1e3x Lyapunov amplification that makes
% single shooting from collocation-derived costates intractable (measured:
% catalog seeds miss by 36,000-560,000 km when single-shot). The seed is a
% full state+costate TRAJECTORY, e.g. from a direct collocation solution.
%
% Since migration #3 this is a thin PROBLEM DEFINITION bound to the generic
% engine costate_common/ms_bvp: this file owns only the CR3BP min-time
% closures (propagation via tfMinProp, dynamics and terminal conditions via
% tfMinEoM) and the z8 packing; structure, Jacobian assembly, guards, and
% the trust-region solve live in the engine. All dynamics, STMs, and the
% Hamiltonian come from pumpkyn calls so every convention matches tfMin
% exactly; a converged z here is interchangeable with a tfMin solution.
%
% INPUTS:
%   rv0    - initial rotating-frame pos/vel, ND [6x1 or 1x6]
%   rvf    - final   rotating-frame pos/vel, ND [6x1 or 1x6]
%   seed   - struct:
%              .tf     initial time-of-flight guess, ND [scalar]
%              .tGrid  segment-boundary times, 0..tf [1 x K+1]
%              .Y      augmented states [r;v;m;lam_r;lam_v;lam_m] at the
%                      boundary times [14 x K+1] (col 1 state part is
%                      overwritten with [rv0; 1])
%   Tmax   - max thrust accel, ND [scalar]
%   c      - exhaust velocity, ND [scalar]
%   muStar - CR3BP mass ratio [scalar]
%   opts   - (optional) struct: .maxIter [100], .tolR [1e-10], .wallSec
%            [300], .verbose [false], .keepSTMs [false] (adds info.PHI for
%            ms_conjugate_test), .conjTest [false] (implies keepSTMs; runs
%            the free-time quotiented Jacobi test and attaches info.conj)
%
% OUTPUTS:
%   z    - [lambda0(7); tf] in tfMin's convention [8x1] (best iterate)
%   info - struct: .converged, .normR (final inf-norm), .iters, .wall (s),
%          .Y (14 x K junction START states; the final propagated endpoint
%          is not included), .tGrid, .PHI {1 x K} (only if keepSTMs),
%          .conj (ms_conjugate_test output, only if conjTest)
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMin / tfMinProp / tfMinEoM (D. Koblick, Coorbital) -
%       single-shooting original; conventions and dynamics reused verbatim.
%   [2] Betts, "Practical Methods for Optimal Control", ch. 3 (multiple
%       shooting structure).

if nargin < 7, opts = struct(); end

% Self-resolve the shared-library dependency: many existing callers add
% only indirect/ to the path (review finding, GPT 2026-08-08).
if isempty(which('ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'costate_common'));
end

rv0 = rv0(:);  rvf = rvf(:);

prob = struct('ny', 14, 'freeIdx0', 8:14, ...
    'prop',     @(dt, y0, needSTM) propSeg(dt, y0, needSTM, Tmax, c, muStar), ...
    'rhs',      @(y) rhsPoint(y, Tmax, c, muStar), ...
    'terminal', @(y, needJ) terminalMinTime(y, rvf, Tmax, c, muStar));

conjTest = isfield(opts, 'conjTest') && opts.conjTest;
if conjTest, opts.keepSTMs = true; end

seed.Y(1:7,1) = [rv0; 1];                      % state part of col 1 is fixed
[p, info] = ms_bvp(prob, seed, opts);

z = [p(1:7); p(end)];

if conjTest
    % Free-time quotiented Jacobi test: columns are (lam_r, lam_v)
    % variations modulo the scaling invariance; the flow column closes the
    % square. lam_m is excluded (control never depends on it).
    info.conj = ms_conjugate_test(info, ...
        struct('flow', @(y) flow6(y, Tmax, c, muStar)));
end
end

% ------------------------------------------------------------------------
function [yh, PHI] = propSeg(dt, y0, needSTM, Tmax, c, muStar)
% PROPSEG  One-segment propagation via pumpkyn tfMinProp, with STM.
% INPUTS: dt; y0 [14x1]; needSTM logical; Tmax; c; muStar.
% OUTPUTS: yh [14x1]; PHI [14x14] or [].
if needSTM, y0 = [y0; reshape(eye(14), [], 1)]; end
[~, Yout] = pumpkyn.cr3bp.tfMinProp(dt, y0, Tmax, c, muStar);
yh = Yout(end, 1:14)';
if needSTM, PHI = reshape(Yout(end, 15:210), 14, 14); else, PHI = []; end
end

% ------------------------------------------------------------------------
function F = rhsPoint(y, Tmax, c, muStar)
% RHSPOINT  Min-time dynamics at a point via pumpkyn tfMinEoM.
% INPUTS: y [14x1]; Tmax; c; muStar.  OUTPUTS: F [14x1].
F = pumpkyn.cr3bp.tfMinEoM(0, [y; reshape(eye(14), [], 1)], Tmax, c, muStar);
F = F(1:14);
end

% ------------------------------------------------------------------------
function f = flow6(y, Tmax, c, muStar)
% FLOW6  Position/velocity rows of the min-time dynamics at a point.
% INPUTS: y [14x1]; Tmax; c; muStar.  OUTPUTS: f [6x1].
F = rhsPoint(y, Tmax, c, muStar);
f = F(1:6);
end

% ------------------------------------------------------------------------
function [g, dgdy] = terminalMinTime(y, rvf, Tmax, c, muStar)
% TERMINALMINTIME  Free-tf min-time terminal conditions and Jacobian:
% r,v matched; lambda_m(tf) = 0; H(tf) = 0.
% INPUTS: y [14x1]; rvf [6x1]; Tmax; c; muStar.
% OUTPUTS: g [8x1]; dgdy [8x14].
[~, H, dHdy] = pumpkyn.cr3bp.tfMinEoM(0, [y; reshape(eye(14), [], 1)], ...
                                      Tmax, c, muStar);
g = [y(1:6) - rvf; y(14); H];
dgdy = [eye(6), zeros(6, 8);
        zeros(1, 13), 1;
        dHdy];
end
