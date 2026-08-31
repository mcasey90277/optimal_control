function [tD, rvD, tT, rvT] = ladder_endpoints(ob)
%% Purpose:
%
%   Rebuilds a thrust-ladder sheet's departure and arrival orbits from its
%   OWN meta recipe, family-agnostically -- the endpoint provider for
%   extend_thrust_ladder (and any future ladder tool that reads engine-
%   native sheet files). Mirrors densify_ladder's proven meta-dispatch
%   pattern (diagnostic-A, 2026-08-26):
%
%     - a sheet whose meta carries .depFamily/.depParams (every campaign
%       since thrust_ladder_library went family-agnostic: halo, dpo,
%       halo_halo, gto) routes through the shared
%       costate_common/get_family_orbit provider using the sheet's own
%       recipes;
%     - a sheet built before those fields existed (the DRO fine sheet and
%       the DRO catalog sheets) falls back to the ORIGINAL hardcoded
%       DRO/tulip construction verbatim, so old sheets rebuild
%       byte-identically to the code they were solved with.
%
%% Inputs:
%
%  ob                       struct                  A sheet's meta struct
%                                                   (thrust_ladder_library
%                                                   format: muStar, tauDRO,
%                                                   NpTulip, tauTulip,
%                                                   pmTulip, and optionally
%                                                   depFamily/depParams/
%                                                   arrFamily/arrParams)
%
%% Outputs:
%
%  tD                       [N x 1]                 Departure-orbit time grid
%                                                   (ND, one period)
%
%  rvD                      [N x 6]                 Departure-orbit states
%                                                   (ND rotating frame)
%
%  tT                       [M x 1]                 Arrival-orbit time grid
%                                                   (ND, one period)
%
%  rvT                      [M x 6]                 Arrival-orbit states
%                                                   (ND rotating frame)
%
%% Revision History:
%  M. Casey                                                   (c) 08/31/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: rebuild the DRO fine sheet's endpoint pair and plot both orbits:
     here = fileparts(mfilename('fullpath'));
        Q = load(fullfile(here, '..', 'direct', 'results', ...
                          'thrust_ladder_12x12.mat'), 'meta');
[tD, rvD, tT, rvT] = ladder_endpoints(Q.meta);
     figure('color',[1 1 1]);
     plot3(rvD(:,1), rvD(:,2), rvD(:,3), 'k'); hold on;
     plot3(rvT(:,1), rvT(:,2), rvT(:,3), 'r');
     axis equal; grid on; legend('departure','arrival'); view(-35,25);
     return;
end

muStar = ob.muStar;

%% Departure orbit -- sheet recipe if present, legacy DRO rule otherwise:
if isfield(ob, 'depFamily') && isfield(ob, 'depParams')
    depParams = ob.depParams;
    if ~isfield(depParams, 'muStar'), depParams.muStar = muStar; end
    [tD, rvD] = get_family_orbit(ob.depFamily, depParams);
else
    [~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
         rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
    [tD, rvD] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
end

%% Arrival orbit -- same dispatch:
if isfield(ob, 'arrFamily') && isfield(ob, 'arrParams')
    arrParams = ob.arrParams;
    if ~isfield(arrParams, 'muStar'), arrParams.muStar = muStar; end
    [tT, rvT] = get_family_orbit(ob.arrFamily, arrParams);
else
    [~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
         rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
    [tT, rvT] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);
end
end
