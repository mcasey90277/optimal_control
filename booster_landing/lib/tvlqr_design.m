function ctrl = tvlqr_design(sol, P, opts)
% TVLQR_DESIGN  Time-varying LQR about a PDG guidance trajectory.
%
% Linearizes pdg_dynamics along (x*(t), T*(t)) and integrates the
% differential Riccati equation backward:
%   -Pdot = A'P + PA - PBR^{-1}B'P + Q,   P(tf) = Qf
% Gains K(t) = R^{-1} B(t)' P(t), stored on a dense grid; controller is
%   T_cmd(t) = T*(t) - K(t) (x - x*(t))  (saturation applied by the sim).
%
% INPUTS:
%   sol  - collocation solution (Task 3 interface)
%   P    - booster_params
%   opts - (optional) .Q .R .Qf weight overrides (7x7, 3x3, 7x7)
% OUTPUTS:
%   ctrl - .tgrid(1xM) .K(3x7xM) .Pt(7x7xM) .xnom(@t->7x1)
%          .Tnom(@t_scalar->3x1) -- SCALAR-t ONLY (see note below); every
%          caller in this repo (ricrhs below, sim_closed_loop's
%          control_law) already calls it with a scalar t.
%
% ADAPTATION (task-7 fix report, 2026-08-08): ctrl.Tnom was originally a
% global pchip spline threaded through all nodes+midpoints. certify_pdg.m's
% G2 gate already disproved that exact representation for reconstructing a
% Hermite-Simpson control (task-5 fix report): a global spline is not the
% per-segment quadratic Simpson's rule is actually built against, and
% measured a non-shrinking ~0.6-0.8 kg mass "floor" that the per-segment
% reconstruction collapsed to ~0.0001-0.0005 kg. ctrl.Tnom now calls the
% shared lib/hs_quad_ctrl.m (promoted out of certify_pdg.m, byte-identical
% logic, regression-checked by test_certify_nominal.m) instead, so the
% TVLQR feedforward flies the same control G2 already certified accurate
% to ~1.6e-4 m, not a spline known to smear the bang-bang switch and
% front-load braking near touchdown. hs_quad_ctrl is scalar-t only (see
% its own header) -- ctrl.Tnom inherits that restriction; ctrl.xnom stays
% a pchip state interpolant (states are smooth between HS nodes, no
% representation mismatch there).
%
% ADAPTATION 2 (task-7 fix report round 2, 2026-08-08): default R raised
% from 1e-10 to 1.5e-9 at the time, under a v(tf)=0 terminal BC. Root
% cause: under the deep terminal-arc saturation this trajectory forces
% (P.Tmin exceeds vehicle weight at every mass -- see sim_closed_loop.m),
% a saturated control law only has DIRECTION left to work with, not
% magnitude; at R=1e-10 the direction is chosen almost entirely to null a
% large dispersed horizontal error, leaving little of the saturated
% thrust vector's direction "budget" for arresting vertical speed.
%
% ADAPTATION 3 (task-7 fix report round 3, 2026-08-08): default R lowered
% again, to 7e-10, after P.vf changed the guidance terminal velocity from
% 0 to [0;0;-1.5] m/s (booster_params.m; the v=0 endpoint was singular
% under P.Tmin>weight). A fresh mini-sweep under the new BC (R from 1e-12
% to 2e-9, plus a Qf(4:6) probe at R=7e-10 -- task-7 fix report round 3
% for the full table) found the dispersed-50m-safe window narrowed to
% R in [6.5e-10, 7.0e-10] (R<=6e-10 fails vtd>2; R>=7.5e-10 reverts to an
% arrest, td.landed=false, even though its vtd number looks fine).
% R=7e-10 sits at the good end of that window with markedly more margin
% than round 2's R=1.5e-9 gave under the old BC (dispersed miss=9.64 m
% vs pad_radius=15, vtd=0.88 vs vtd_max=2.0 -- round 2's setting was
% razor-thin at miss=14.77/vtd=1.377).
%
% The zero-dispersion case STILL does not reach a genuine touchdown at
% ANY R tested in the dispersed-safe window (arrest ~0.65 m up, vtd~0.64
% m/s, far from the P.vf=-1.5 m/s target) -- removing the v(tf)=0
% singularity substantially improved nominal tracking (a LOW R, ~1e-11 to
% 7e-11, now DOES achieve a genuine landing near the -1.5 m/s target) but
% did not eliminate the arrest possibility under closed-loop tracking
% error, and that low-R band is disjoint from the dispersed-safe band by
% over an order of magnitude. See the task-7 fix report round 3 for the
% full sweep and the resulting honest test status.
%
% REFERENCES:
%   [1] Anderson & Moore, "Optimal Control: Linear Quadratic Methods."
if nargin < 3, opts = struct(); end
if ~isfield(opts,'Q'),  opts.Q  = diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0]); end
if ~isfield(opts,'R'),  opts.R  = 7e-10*eye(3);                            end
if ~isfield(opts,'Qf'), opts.Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);          end

%% Nominal interpolants (state: pchip over nodes; thrust: shared
%% per-segment HS quadratic reconstruction, scalar t only -- see the
%% ADAPTATION note above):
Nn = size(sol.X,2) - 1;  h = sol.tf/Nn;
ctrl.xnom = @(t) interp1(sol.t.', sol.X.', clampt(t, sol.tf), 'pchip').';
ctrl.Tnom = @(t) hs_quad_ctrl(clampt(t, sol.tf), sol.U, sol.Um, h, Nn);

%% Backward Riccati on vec(P), dense output grid:
M     = 4*(Nn+1);
ctrl.tgrid = linspace(0, sol.tf, M);
Rinv  = inv(opts.R);
ric   = @(t, pv) ricrhs(t, pv, ctrl, P, opts.Q, Rinv);
oo    = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, PV] = ode45(ric, fliplr(ctrl.tgrid), opts.Qf(:), oo);   % integrates tf->0
PV    = flipud(PV);                                          % re-order 0->tf
ctrl.Pt = zeros(7,7,M);  ctrl.K = zeros(3,7,M);
for k = 1:M
    Pk = reshape(PV(k,:), 7, 7);  Pk = (Pk + Pk.')/2;
    [~, ~, Bk] = pdg_dynamics(ctrl.xnom(ctrl.tgrid(k)), ...
                              ctrl.Tnom(ctrl.tgrid(k)), P);
    ctrl.Pt(:,:,k) = Pk;
    ctrl.K(:,:,k)  = Rinv * (Bk.' * Pk);
end
end

function t = clampt(t, tf), t = min(max(t, 0), tf); end

function pdot = ricrhs(t, pv, ctrl, P, Q, Rinv)
Pk = reshape(pv, 7, 7);
[~, A, B] = pdg_dynamics(ctrl.xnom(t), ctrl.Tnom(t), P);
Pd   = -(A.'*Pk + Pk*A - Pk*B*Rinv*B.'*Pk + Q);
pdot = Pd(:);
end
