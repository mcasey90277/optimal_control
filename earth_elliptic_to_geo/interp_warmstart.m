function W = interp_warmstart(Xsrc, Usrc, dLsrc, sigmaSrc, sigmaDst)
% INTERP_WARMSTART  Pure mesh-refine helper: interpolate a converged MEE
% trajectory (Xsrc/Usrc, sampled on sigmaSrc in [0,1]) onto a new sigma grid
% (sigmaDst), producing a warm-start point for casadi_lt_mee.m at a different
% node density/revolution count. Factored out (review finding, Fix 2) from
% THREE previously-inline copies of the exact same interp1 pattern --
% run_transfer_mee.m's cfg.warmStart path, run_mintime_mee.m's Stage-B
% cfg.warmStartAnchor path, and nodestudy_mee.m>solve_warm_node -- so the
% logic is unit-testable without a solve (test_warmstart_mee.m) and the two
% ladder-critical call sites can no longer drift apart.
%
% Interpolation convention (unchanged from the pre-refactor inline code):
% LINEAR for the continuous state X and the RTN thrust-direction rows of U
% (rows 1-3, "beta"), NEAREST for the throttle row (row 4, "thr") to keep
% bang-bang switch edges crisp instead of blurring them into intermediate
% throttle values. dL is a PASSTHROUGH, not interpolated -- it is a single
% scalar total-longitude-span value, invariant of node count.
%
% LATENT BUG FIX (review finding, Fix 2): linearly interpolating two unit
% RTN thrust-direction vectors (beta, constrained to |beta|=1 at every node
% in the actual casadi_lt_mee.m optimization, see its
% "beta(1,k)^2+beta(2,k)^2+beta(3,k)^2==1" constraint) does NOT in general
% produce a unit vector at the interpolated point -- especially near a
% throttle switch, where beta can rotate quickly between adjacent source
% nodes. The pre-refactor code fed this sub-unit Ubeta straight into U0
% without renormalizing, silently handing casadi_lt_mee.m a warm start that
% violates its own unit-norm constraint at every interpolated node (a primal
% infeasibility at iteration 0, the exact class of problem the STAGE-B SEED
% THROTTLE FIX and ANCHOR WARM-START FIX comments elsewhere in this
% directory went to some lengths to avoid). FIX: renormalize each column of
% the interpolated beta rows to unit norm before returning.
%
% INPUTS:  Xsrc     - source state trajectory [7 x Msrc]
%          Usrc     - source control trajectory [4 x Msrc], row 4 = throttle
%          dLsrc    - source total longitude span [scalar, ND] (passthrough)
%          sigmaSrc - source sigma grid, [0,1] [Msrc x 1]
%          sigmaDst - destination sigma grid, [0,1] [Mdst x 1]
% OUTPUTS: W - struct: .X [7 x Mdst] (linear-interpolated), .U [4 x Mdst]
%          (rows 1-3 linear-interpolated then renormalized to unit norm,
%          row 4 nearest-interpolated), .dL [scalar, = dLsrc unchanged]
%
% REFERENCES: [1] run_transfer_mee.m (cfg.warmStart path, one of two
%   ladder-critical callers). [2] run_mintime_mee.m (cfg.warmStartAnchor
%   Stage-B path, the other ladder-critical caller). [3] nodestudy_mee.m
%   >solve_warm_node (the third, pre-existing copy of this pattern; not
%   itself refactored to call this helper -- out of this fix's scope, noted
%   for a future pass). [4] casadi_lt_mee.m (the unit-norm beta constraint
%   this fix restores compatibility with). [5] test_warmstart_mee.m (unit
%   tests: sizes, endpoint preservation, nearest-throttle bang-bang
%   preservation, unit-norm renormalization, dL passthrough).
W.X = interp1(sigmaSrc, Xsrc.', sigmaDst, 'linear').';

Ubeta = interp1(sigmaSrc, Usrc(1:3,:).', sigmaDst, 'linear').';
betaNorm = sqrt(sum(Ubeta.^2, 1));
Ubeta = Ubeta ./ betaNorm;   % restore |beta|=1 (see LATENT BUG FIX above)

Uthr = interp1(sigmaSrc, Usrc(4,:).', sigmaDst, 'nearest').';

W.U  = [Ubeta; Uthr];
W.dL = dLsrc;
end
