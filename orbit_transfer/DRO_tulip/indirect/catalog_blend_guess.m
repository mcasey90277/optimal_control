function [zGuess, blend] = catalog_blend_guess(K, sD, sA, opts)
%% Purpose:
%
%   The costate guess for an off-grid phase pair, as a consumer wants it: the
%   interpolator's blend (phase_catalog_interp, with the time-consistency
%   rule switched on from the setup's sensitivity maps), and the guess's t_f
%   taken from what the blended corners' own costates predict,
%
%       t_f = sum_k  w_k ( T_k + dT/ds_D|_k (s_D - s_D,k) + dT/ds_A|_k (s_A - s_A,k) ),
%
%   rather than from the blend of their t_f. The costates carry dT/ds, so
%   this uses information the plain blend throws away; where t_f bends
%   sharply between grid lines it is closer by hours (FINDINGS 89).
%
%% Inputs:
%
%  K                        struct                  catalog_blend_setup output
%  sD, sA                   double                  the query phases
%  opts                     struct (optional)       .tfFromSensitivities [true]
%                                                   .maxJump, .nearestOnly
%                                                   (passed to the
%                                                   interpolator)
%
%% Outputs:
%
%  zGuess                   [8 x 1]                 NaN when tier 'none'
%  blend                    struct                  phase_catalog_interp's
%                                                   info, plus .tfBlend (the
%                                                   plain blend's t_f) and
%                                                   .tfCorner [1 x n] each
%                                                   corner's prediction
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
fromG = true;  if isfield(opts, 'tfFromSensitivities'), fromG = opts.tfFromSensitivities; end
io = struct('dTdsD', K.GD, 'dTdsA', K.GA, 'maxEdgeResidual', K.edgeCap);
if isfield(opts, 'maxJump'),     io.maxJump = opts.maxJump;         end
if isfield(opts, 'nearestOnly'), io.nearestOnly = opts.nearestOnly; end
[zGuess, blend] = phase_catalog_interp(K.sheet, sD, sA, io);
blend.tfBlend = zGuess(8);  blend.tfCorner = [];
if strcmp(blend.tier, 'none'), return, end

shortWay = @(ds) mod(ds + 0.5, 1) - 0.5;
tfCorner = NaN(1, numel(blend.weights));
for kc = 1:numel(blend.weights)
    iD = blend.iD(kc);  iA = blend.iA(kc);
    tfCorner(kc) = K.sheet.tf_nd(iD, iA) + K.GD(iD, iA)*shortWay(sD - K.sheet.sD_frac(iD)) ...
                                         + K.GA(iD, iA)*shortWay(sA - K.sheet.sA_frac(iA));
end
blend.tfCorner = tfCorner;
if fromG && all(isfinite(tfCorner)), zGuess(8) = blend.weights(:).'*tfCorner(:); end
end
