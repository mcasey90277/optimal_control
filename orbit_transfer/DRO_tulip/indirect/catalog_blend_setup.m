function K = catalog_blend_setup(cat_)
%% Purpose:
%
%   Everything an interpolation needs from a phase catalog, made ONCE: the
%   sheet, the engine in ND units, the two endpoint closures rebuilt from the
%   catalog's own keys, the phase sensitivities dT/ds_D and dT/ds_A of every
%   entry (one flight each), the grid steps, and the time-consistency cap
%   that goes with those steps.
%
%   interp_study does these steps in the open, for one query. A consumer that
%   asks many queries (score_interpolator) does them here.
%
%  ASSUMPTIONS / NOTES:
%
% • A single-sheet DRO -> tulip phase catalog (ladder_endpoints' orbit keys).
% • The cap is 300 min at a grid step of 1/24, scaled as step^1.5 per axis:
%   between the residual it must admit (third order in the step) and the
%   jump it must refuse (the same at any step). FINDINGS 89.
% • A phase with ONE grid line has step 0 and cap 0 (nothing to compare).
%
%% Inputs:
%
%  cat_                     char | struct           catalog file, or the
%                                                   catalog struct
%
%% Outputs:
%
%  K                        struct                  .sheet .phys (.Tnd .cnd
%                                                   .muStar .lStar .tStar
%                                                   .m0kg) .stateD .stateA
%                                                   .dstateD .dstateA
%                                                   (closures of a phase)
%                                                   .GD .GA [nD x nA] ND time
%                                                   per unit phase (NaN: no
%                                                   entry) .gapD .gapA
%                                                   .edgeCap [1 x 2] ND time
%                                                   .catMat ('' for a struct)
%                                                   .problem (certify_root's
%                                                   B.problem)
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

catMat = '';
if ischar(cat_) || isstring(cat_)
    catMat = char(cat_);  L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});
end
assert(isscalar(cat_.sheets), 'catalog_blend_setup:sheets', 'a single-sheet phase catalog is expected');
sheet = cat_.sheets(1);
muStar = cat_.constants.muStar;  lStar = cat_.constants.lStar_km;  tStar = cat_.constants.tStar_s;
thrustN = cat_.rungs_N(1);  ispS = cat_.thruster.isp_s;  m0kg = cat_.thruster.m0_kg;
ndp = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
phys = struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', muStar, 'lStar', lStar, 'tStar', tStar, 'm0kg', m0kg);

keys = struct('muStar', muStar, 'lStar', lStar, 'tStar', tStar, 'tauDRO', sheet.tauDRO, 'NpTulip', sheet.Np, ...
              'tauTulip', sheet.period_tulip_nd, 'pmTulip', sheet.pm, 'ispS', ispS, 'm0kg', m0kg);
[tD, rvD, tT, rvT] = ladder_endpoints(keys);
[stateD, ~, dstateD] = phase_state(tD, rvD);
[stateA, ~, dstateA] = phase_state(tT, rvT);

% the sensitivities of every entry: fly it once for lambda(t_f), two dot products
has = logical(sheet.has_solution(:, :, 1));
GD = NaN(size(has));  GA = NaN(size(has));
for iD = find(any(has, 2)).'
    x0 = stateD(sheet.sD_frac(iD));  dxD = dstateD(sheet.sD_frac(iD));
    for iA = find(has(iD, :))
        z8 = sheet.z8(:, sheet.entry_index(iD, iA));
        [~, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [x0(1:6); 1; z8(1:7)], phys.Tnd, phys.cnd, muStar);
        dxA = dstateA(sheet.sA_frac(iA));
        [GD(iD, iA), GA(iD, iA)] = phase_sensitivity(z8(1:6), Y(end, 8:13).', dxD(1:6), dxA(1:6));
    end
end

gapD = typicalStep(sheet.sD_frac);  gapA = typicalStep(sheet.sA_frac);
edgeCap = (300*60/tStar)*(24*[gapD gapA]).^1.5;

K = struct('sheet', sheet, 'phys', phys, 'stateD', stateD, 'stateA', stateA, 'dstateD', dstateD, 'dstateA', dstateA, ...
           'GD', GD, 'GA', GA, 'gapD', gapD, 'gapA', gapA, 'edgeCap', edgeCap, 'catMat', catMat, ...
           'problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, 'thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
                             'tauDRO', sheet.tauDRO, 'NpTulip', sheet.Np, 'pmTulip', sheet.pm));
end

% ==========================================================================
function h = typicalStep(phases)
% TYPICALSTEP  The median gap between neighbouring grid lines round the circle
% (0 for a single line).  INPUTS: phases [1 x n].  OUTPUTS: h (double).
p = sort(mod(phases(:).', 1));
if isscalar(p), h = 0;  return, end
h = median(diff([p, p(1) + 1]));
end
