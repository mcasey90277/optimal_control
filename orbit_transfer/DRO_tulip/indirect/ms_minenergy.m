function [z, info] = ms_minenergy(rv0, rvf, tf, seed, Tmax, c, muStar, opts)
% MS_MINENERGY  Multiple-shooting solve of the CR3BP fixed-tf MINIMUM-ENERGY
% PMP problem -- the min-energy sibling of ms_tfmin.
%
% Cost J = Int s^2 dt (Bertrand-Epenoy energy endpoint, GTO_tulip's eps = 1
% convention), final time FIXED at tf, final mass free. Terminal conditions:
% r(tf) = rf, v(tf) = vf, lambda_m(tf) = 0 (seven, for the seven unknown
% initial costates); there is no H(tf) = 0 condition -- instead H is a first
% integral along the arc and is reported for checking (info.H).
%
% Thin PROBLEM DEFINITION bound to costate_common/ms_bvp with opts.fixedTf:
% this file owns only the three closures (propagation + STM via
% cr3bp_minenergy_prop, point dynamics + Hamiltonian via
% cr3bp_minenergy_pmp, the terminal conditions) and the z7 packing.
% Optionally runs the generic single-shooting acceptance gate
% (costate_common/ss_bvp_accept) on the same closures -- the role pumpkyn
% tfMin plays for min-time entries, for which no pumpkyn twin exists here.
%
% INPUTS:
%   rv0    - initial rotating-frame pos/vel, ND [6x1 or 1x6]
%   rvf    - final   rotating-frame pos/vel, ND [6x1 or 1x6]
%   tf     - FIXED time of flight, ND [scalar]
%   seed   - struct:
%              .tGrid  segment-boundary times, 0..tf [1 x K+1]
%              .Y      augmented states [r;v;m;lam_r;lam_v;lam_m] at the
%                      boundary times [14 x K+1] (col 1 state part is
%                      overwritten with [rv0; 1]); .tf is ignored/overwritten
%   Tmax   - max thrust accel, ND [scalar]
%   c      - exhaust velocity, ND [scalar]
%   muStar - CR3BP mass ratio [scalar]
%   opts   - (optional) struct: ms_bvp options (.maxIter [100], .tolR
%            [1e-10], .wallSec [300], .verbose [false], .keepSTMs [false])
%            plus .accept [false] (run ss_bvp_accept on the result; attaches
%            info.accept), .tolDz [1e-6] (acceptance move tolerance),
%            .tolRss [1e-6] (single-shooting residual floor for the gate;
%            NOT .tolR -- see ss_bvp_accept)
%
% OUTPUTS:
%   z    - lambda0 = [lam_r; lam_v; lam_m] [7x1] (best iterate)
%   info - struct: .converged, .normR, .iters, .wall, .Y (14 x K junction
%          START states), .tGrid, .H [1 x K] Hamiltonian at the junctions
%          (with the s^2 running cost; constant on a converged extremal),
%          .Hdrift = max |H - H(1)|, .s [1 x K] throttle at the junctions,
%          .PHI (keepSTMs), .accept (opts.accept: ss_bvp_accept info)
%
% REFERENCES:
%   [1] costate_common/cr3bp_minenergy_pmp (the field; unit-tested against
%       pumpkyn.cr3bp.tfMinEoM at saturation).
%   [2] Bertrand & Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," OCAM 23(4), 2002 (the energy endpoint).
%   [3] Betts, "Practical Methods for Optimal Control", ch. 3.

if nargin < 8, opts = struct(); end

% Self-resolve the shared-library dependency (as ms_tfmin does).
if isempty(which('ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'costate_common'));
end

rv0 = rv0(:);  rvf = rvf(:);

prob = struct('ny', 14, 'freeIdx0', 8:14, ...
    'prop',     @(dt, y0, needSTM) cr3bp_minenergy_prop(dt, y0, needSTM, Tmax, c, muStar), ...
    'rhs',      @(y) cr3bp_minenergy_pmp(y, Tmax, c, muStar), ...
    'terminal', @(y, needJ) terminalMinEnergy(y, rvf));

doAccept = isfield(opts, 'accept') && opts.accept;
tolDz = 1e-6;   tolRss = 1e-6;
if isfield(opts, 'tolDz'),  tolDz  = opts.tolDz;  end
if isfield(opts, 'tolRss'), tolRss = opts.tolRss; end
mo = opts;
for f = {'accept', 'tolDz', 'tolRss'}
    if isfield(mo, f{1}), mo = rmfield(mo, f{1}); end
end
mo.fixedTf = true;

conjTest = isfield(opts, 'conjTest') && opts.conjTest;
if conjTest, mo.keepSTMs = true; end
if isfield(mo, 'conjTest'), mo = rmfield(mo, 'conjTest'); end

seed.tf = tf;
seed.Y(1:7,1) = [rv0; 1];                      % state part of col 1 is fixed
[p, info] = ms_bvp(prob, seed, mo);
z = p(1:7);

if conjTest
    % FIXED-tf Jacobi test (validated on the analytic LQ pi-conjugate case,
    % tests/test_conj_fixedtf): no flow column, no scaling quotient (L
    % breaks the invariance), rows = the components vanishing under the
    % terminal conditions (r, v fixed; lam_m(tf) = 0), cols = ALL seven
    % initial costates (the throttle law depends on lam_m here).
    info.conj = ms_conjugate_test(info, struct( ...
        'stateRows', [1:6 14], 'costateCols', 8:14, ...
        'quotientDir', [], 'freeTime', false));
end

% Hamiltonian and throttle at the junction starts (first-integral check).
K = size(info.Y, 2);
info.H = zeros(1, K);  info.s = zeros(1, K);
for k = 1:K
    [~, ~, ax] = cr3bp_minenergy_pmp(info.Y(:,k), Tmax, c, muStar);
    info.H(k) = ax.H;  info.s(k) = ax.s;
end
info.Hdrift = max(abs(info.H - info.H(1)));

if doAccept
    y1 = info.Y(:,1);
    ao = struct('fixedTf', true, 'tolDz', tolDz, 'tolR', tolRss);
    [~, info.accept] = ss_bvp_accept(prob, y1, tf, ao);
end
end

% ------------------------------------------------------------------------
function [g, dgdy] = terminalMinEnergy(y, rvf)
% TERMINALMINENERGY  Fixed-tf min-energy terminal conditions and Jacobian:
% r, v matched; lambda_m(tf) = 0 (final mass free).
% INPUTS: y [14x1]; rvf [6x1].  OUTPUTS: g [7x1]; dgdy [7x14].
g = [y(1:6) - rvf; y(14)];
dgdy = [eye(6), zeros(6, 8);
        zeros(1, 13), 1];
end
