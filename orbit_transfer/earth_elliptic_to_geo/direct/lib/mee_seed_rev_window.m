function w = mee_seed_rev_window()
% MEE_SEED_REV_WINDOW  Admissible revolution count for the two-pass MEE seed.
%
% The two-pass seed protocol runs a cheap N=50 mee_seed probe, reads its
% revolution count, and sizes the real seed as N = round(nodesPerRev * nRev).
% This is the window the probe's nRev must land in for that seed to be usable.
%
% WHY IT IS A WINDOW, AND WHY THESE NUMBERS. At full throttle the MEE seed ODE
% reaches P = 1 in ~3.06 revs and goes coordinate-singular past ~rev 4 -- it
% chases an e -> 1 singularity once P overshoots the target. Throttling to
% ~0.4 with tangential steering stretches the crossing to ~7.5 revs, which both
% matches the paper's optimal revolution count and gives a well-conditioned
% collocation seed. Landing outside [6.5, 9] means seedThr is wrong for this
% thrust, and the right response is to fix seedThr, not to widen this window.
%
% WHY IT IS SHARED RATHER THAN INLINED. Two drivers implement the two-pass
% protocol -- run_transfer_mee.m (earth) and bridge_mu_continuation.m (CR3BP) --
% and the protocol itself is deliberately NOT factored out: the four lines they
% share are surrounded by caching and knob-sourcing that each driver legitimately
% owns and that have already diverged (registry lookup vs cfg, resume-gated vs
% unconditional cache reads). Extracting those four lines would need a
% seven-argument helper -- a worse trade than the duplication.
%
% What must NOT diverge is the physics: this window and the N formula. A change
% applied to one driver and not the other would silently give the two campaigns
% different admissible seeds. So the window lives here, in the one place both
% drivers already reach (CR3BP delegates its path to the earth campaign), while
% the orchestration stays where it belongs.
%
% INPUTS:  (none)
% OUTPUTS: w - [nRevMin nRevMax] admissible revolution count [1x2]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/drivers/run_transfer_mee.m (stage 1).
%   [2] earth_elliptic_to_geo_CR3BP/direct/bridge_mu_continuation.m (stage 1).
%   [3] orbit_transfer/CODE_STRUCTURE.md Tier 2 (why the protocol is not
%       extracted but this constant is).
w = [6.5, 9];
end
