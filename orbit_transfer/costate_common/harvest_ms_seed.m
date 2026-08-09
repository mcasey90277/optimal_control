function [seed, diag_] = harvest_ms_seed(o, K)
%% Purpose:
%
%   Builds a multiple-shooting SEED from a direct collocation solution --
%   the harvest path, made a single-home library function (migration #3;
%   the sign-vote + midpoint-association rules previously lived inline in
%   thrust_ladder_library, the exact one-home-per-rule violation that let
%   the Hermite-Simpson midpoint bug exist in two places).
%
%   The covector mapping itself (station association, sign vote, lambda_t
%   check) is delegated to duals_to_costates; this function owns only the
%   interpolation of states and mapped costates onto the K+1 segment
%   boundaries.
%
%  ASSUMPTIONS / NOTES:
%
% • Hermite-Simpson defect multipliers belong at interval MIDPOINTS
%   (review finding, Gemini 2026-08-04) -- enforced by duals_to_costates.
% • With Sundman regularization off the mesh is uniform in time and the
%   solver leaves tNodes empty; it is reconstructed here.
% • Costate interpolation uses pchip with EXTRAP: midpoint stations do not
%   reach t = 0 or t = tf, so the boundary values are extrapolated.
%
%% Inputs:
%
%  o                        struct                  Direct-solve output:
%   .X                      [>=7 x N+1]             States at nodes
%   .lamDef                 [>=7 x N]               Defect multipliers per
%                                                   interval (rows 1:7 map
%                                                   to costates)
%   .Um                     [>=3 x N]               Solved thrust directions
%   .tNodes                 [1 x N+1] or []         Node times ([] =
%                                                   uniform over [0, tf])
%   .tf                     double                  Final time (ND)
%
%  K                        double                  Number of shooting
%                                                   segments
%
%% Outputs:
%
%  seed                     struct                  ms_tfmin seed: .tf,
%                                                   .tGrid [1 x K+1], .Y
%                                                   [14 x K+1]
%
%  diag_                    struct                  duals_to_costates
%                                                   diagnostics (.sign,
%                                                   .voteMargin, .lamT,
%                                                   .lamTOK, .notes)
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

lamD = o.lamDef(1:7,:);
if isempty(o.tNodes)
    tN = linspace(0, o.tf, size(o.X,2));
else
    tN = o.tNodes(:)';
end

spec = struct('scheme', 'hermite-simpson', 'mu', lamD, 'tNodes', tN, ...
              'uDir', o.Um(1:3, 1:size(lamD,2)));
if size(o.lamDef, 1) >= 8
    spec.lamTf = o.lamDef(8,:);                % enables the lambda_t = +1 check
end
[lam, tMid, diag_] = duals_to_costates(spec);

tG = linspace(0, o.tf, K+1);
Xg = interp1(tN, o.X', tG, 'pchip')';
Lg = interp1(tMid, lam', tG, 'pchip', 'extrap')';
seed = struct('tf', o.tf, 'tGrid', tG, 'Y', [Xg(1:7,:); Lg]);
end
