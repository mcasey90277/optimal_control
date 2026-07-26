function p = cr3bp_ipopt_opts(maxIter, warmTight)
% CR3BP_IPOPT_OPTS  The CR3BP-family IPOPT option set (Sundman solvers).
%
% Single source for the solver options shared, verbatim, by the three
% CR3BP-family Sundman transcriptions:
%   GTO_tulip/direct/sundman_minfuel/casadi_minfuel_sundman.m
%   GTO_ELFO/direct/elfo/casadi_energy_freetf.m
%   GTO_ELFO/direct/elfo/casadi_mintime_freetf.m
%
% Measured 2026-07-26 before extraction: the three carried 20 identical option
% assignments; the only difference was one REDUNDANT `mu_strategy='monotone'`
% in the tulip copy, set before the warmTight branch and immediately overridden
% by it in both arms. That redundancy is dropped here (no behavioural change --
% both arms set mu_strategy themselves).
%
% DELIBERATELY NOT SHARED WITH THE EARTH CAMPAIGN. earth's casadi_lt_mee.m has
% ZERO option lines in common with these (measured): it uses a different option
% style and carries settings this family does not, notably
% mumps_pivot_order = 0 -- the AMD-ordering workaround for a hard METIS crash
% at its problem sizes. Forcing one helper across both would either impose that
% workaround where it is not wanted or dilute it into a flag. earth also cannot
% see cr3bp_common, so sharing would require a NEW cross-cutting path
% dependency; see orbit_transfer/CODE_STRUCTURE.md on why that trade is refused.
% Likewise PSR/lib keeps its own frozen copy by design (guarded instead by
% test_psr_vendor_drift.m).
%
% THE TWO WARM-START REGIMES (the substance of this file):
%   warmTight = true   re-solving AT an already-converged point (an eps
%                      sharpening step). Hug the bounds, start from a small
%                      barrier: mu_init 1e-4 and 1e-9 pushes.
%   warmTight = false  a genuine continuation MOVE (a muGain step, or a
%                      t_f-floating re-solve). Default bound push and a larger
%                      initial barrier so IPOPT has room -- pinning tightly on
%                      a real step makes inf_du blow up.
%
% INPUTS:
%   maxIter   - IPOPT iteration cap [scalar]
%   warmTight - true for a sharpening re-solve, false for a continuation move
%               [logical scalar]
%
% OUTPUTS:
%   p - options struct for opti.solver('ipopt', p) [struct], with fields
%       .print_time [logical] and .ipopt.* as described above
%
% REFERENCES:
%   [1] orbit_transfer/CODE_STRUCTURE.md (Tier-0 extraction; the measurement
%       that these three blocks were identical, and why earth/PSR are excluded).
%   [2] cr3bp_common/tests/test_cr3bp_ipopt_opts.m (the byte-identity gate:
%       asserts this helper reproduces each solver's former inline struct
%       exactly, so the extraction is provably behaviour-preserving).
p = struct();
p.print_time      = true;
p.ipopt.max_iter  = maxIter;
p.ipopt.tol       = 1e-7;
p.ipopt.nlp_scaling_method = 'gradient-based';
p.ipopt.linear_solver      = 'mumps';
p.ipopt.print_level        = 5;

if warmTight
    % TIGHT: re-solve AT a converged point (sharpen). Hug bounds, small barrier.
    p.ipopt.mu_strategy                 = 'monotone';
    p.ipopt.warm_start_init_point       = 'yes';
    p.ipopt.mu_init                     = 1e-4;
    p.ipopt.warm_start_bound_push       = 1e-9;
    p.ipopt.warm_start_bound_frac       = 1e-9;
    p.ipopt.warm_start_slack_bound_push = 1e-9;
    p.ipopt.warm_start_slack_bound_frac = 1e-9;
    p.ipopt.warm_start_mult_bound_push  = 1e-9;
else
    % LOOSE: a genuine continuation move (a muGain step, or a t_f-floating
    % re-solve). Monotone barrier, default bound push, larger initial barrier so
    % IPOPT has room to move (tight pinning makes inf_du blow up on a real step).
    p.ipopt.mu_strategy           = 'monotone';
    p.ipopt.warm_start_init_point = 'yes';
    p.ipopt.mu_init               = 0.1;
end
end
