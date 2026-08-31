function ok = test_ladder_endpoints()
%% Purpose:
%
%   Tests ladder_endpoints -- the family-aware endpoint-orbit rebuild used
%   by extend_thrust_ladder (mirrors densify_ladder's proven meta-dispatch
%   pattern). Three cases:
%
%     1. LEGACY meta (no depFamily fields; the DRO fine sheet and the DRO
%        catalog sheets): the rebuild must be BITWISE identical to the
%        original hardcoded getDRO/getTulip + cont_np + prop construction.
%     2. FAMILY meta (a real HALO catalog sheet): the rebuild must be
%        bitwise identical to get_family_orbit on the sheet's own recipes,
%        and must NOT equal the legacy DRO fallback (dispatch proof).
%     3. DRO catalog sheet meta: the fallback path must reproduce the
%        legacy construction for that sheet's tau/Np.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All cases passed
%
%% Revision History:
%  M. Casey                                                   (c) 08/31/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
ind  = fileparts(here);                         % DRO_tulip/indirect
ot   = fileparts(fileparts(ind));               % orbit_transfer
ok   = true;

%% Case 1 -- legacy meta (fine sheet): bitwise vs the original rule:
Qf = load(fullfile(ind, '..', 'direct', 'results', 'thrust_ladder_12x12.mat'), 'meta');
ob = Qf.meta;  muStar = ob.muStar;
[tD, rvD, tT, rvT] = ladder_endpoints(ob);

[~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
     rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
[tDr,rvDr] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
     rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
[tTr,rvTr] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);

ok = check(ok, isequal(tD,tDr) && isequal(rvD,rvDr) && ...
               isequal(tT,tTr) && isequal(rvT,rvTr), ...
           'case 1: legacy fine-sheet meta bitwise vs original rule');

%% Case 2 -- family meta (halo catalog sheet): dispatch + bitwise:
Qh = load(fullfile(ot, 'HALO_tulip', 'direct', 'results', 'catalog', ...
                   'halo_tau1p75_Np12.mat'), 'meta');
oh = Qh.meta;
[tD2, rvD2, tT2, rvT2] = ladder_endpoints(oh);

dp = oh.depParams;  if ~isfield(dp,'muStar'), dp.muStar = oh.muStar; end
ap = oh.arrParams;  if ~isfield(ap,'muStar'), ap.muStar = oh.muStar; end
[tDh, rvDh] = get_family_orbit(oh.depFamily, dp);
[tTh, rvTh] = get_family_orbit(oh.arrFamily, ap);

ok = check(ok, isequal(tD2,tDh) && isequal(rvD2,rvDh) && ...
               isequal(tT2,tTh) && isequal(rvT2,rvTh), ...
           'case 2a: halo meta bitwise vs get_family_orbit recipes');

% dispatch proof: the halo departure is NOT the legacy DRO(tau=1.75):
[~, rvDx] = pumpkynPie.cr3bp.getDRO(oh.tauDRO);
     rvDx = pumpkyn.cr3bp.cont_np(rvDx, oh.tauDRO, oh.muStar, 1e-12);
[~, rvDd] = pumpkyn.cr3bp.prop(oh.tauDRO, rvDx, oh.muStar);
ok = check(ok, ~isequal(rvD2, rvDd), ...
           'case 2b: halo departure differs from the DRO fallback');

%% Case 3 -- DRO catalog sheet (legacy names, no family fields):
dC = dir(fullfile(ind, '..', 'direct', 'results', 'catalog', '*.mat'));
assert(~isempty(dC), 'no DRO catalog sheets found');
Qc = load(fullfile(dC(1).folder, dC(1).name), 'meta');
oc_ = Qc.meta;
ok = check(ok, ~isfield(oc_, 'depFamily'), ...
           'case 3a: DRO catalog sheet is indeed legacy-format');
[tD3, rvD3, tT3, rvT3] = ladder_endpoints(oc_);
[~, rvD0] = pumpkynPie.cr3bp.getDRO(oc_.tauDRO);
     rvD0 = pumpkyn.cr3bp.cont_np(rvD0, oc_.tauDRO, oc_.muStar, 1e-12);
[tDc,rvDc] = pumpkyn.cr3bp.prop(oc_.tauDRO, rvD0, oc_.muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(oc_.tauTulip, oc_.NpTulip, oc_.pmTulip);
     rvT0 = pumpkyn.cr3bp.cont_np(rvT0, oc_.tauTulip, oc_.muStar, 1e-12);
[tTc,rvTc] = pumpkyn.cr3bp.prop(oc_.tauTulip, rvT0, oc_.muStar);
ok = check(ok, isequal(tD3,tDc) && isequal(rvD3,rvDc) && ...
               isequal(tT3,tTc) && isequal(rvT3,rvTc), ...
           'case 3b: DRO catalog meta bitwise vs original rule');

if ok, fprintf('TEST_LADDER_ENDPOINTS: ALL PASS\n');
else,  fprintf('TEST_LADDER_ENDPOINTS: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, label)
% CHECK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label.
% OUTPUTS: ok updated.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
