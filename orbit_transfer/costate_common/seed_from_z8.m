function seed = seed_from_z8(z8, rv0, K, Tmax, c, muStar)
%% Purpose:
%
%   Build a multiple-shooting SEED from a converged (or trusted) z8 =
%   [lam0(7); tf]: fly it end-to-end with pumpkyn's min-time propagator and
%   cut the flight into K+1 junction states. Seeded AT a root, ms_bvp
%   converges in 1-2 Newton iterations -- the pattern behind the conjugate
%   catalog sweep, the golden-cells regression, and the GTO flagship
%   min-time probe, extracted to one home on its third appearance
%   (migration rule).
%
%  ASSUMPTIONS / NOTES:
%
% • Min-time all-burn flight (tfMinProp): mass fraction starts at 1, the
%   augmented state is [rv; m; lam] (14).
% • Junction states are pchip-interpolated from the flight samples on a
%   UNIFORM tGrid -- the convention every consumer already used; duplicate
%   propagator time samples are removed before interpolation.
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf], tfMin
%                                                   convention
%
%  rv0                      [6 x 1] or [1 x 6]      Departure state,
%                                                   rotating frame ND
%
%  K                        double                  Segment count (K+1
%                                                   junctions)
%
%  Tmax                     double                  ND thrust accel at unit
%                                                   mass fraction
%
%  c                        double                  ND exhaust velocity
%
%  muStar                   double                  CR3BP mass ratio
%
%% Outputs:
%
%  seed                     struct                  .tf, .tGrid [1 x K+1],
%                                                   .Y [14 x K+1] -- the
%                                                   ms_bvp / ms_tfmin seed
%                                                   contract
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

tG = linspace(0, z8(8), K+1);
[tj, yj] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(:); 1; z8(1:7)], Tmax, c, muStar);
[tu, iu] = unique(tj);
Y = interp1(tu, yj(iu, 1:14), tG, 'pchip')';
seed = struct('tf', z8(8), 'tGrid', tG, 'Y', Y);
end
