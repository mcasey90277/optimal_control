function [rv0, rvf, p] = dro_tulip_endpoints(opts)
% DRO_TULIP_ENDPOINTS  Departure and arrival states for the DRO->tulip transfer.
%
% Built from pumpkyn's own catalogs rather than hard-coded, per the standing
% rule: getter first, then cont_np to close the periodic orbit. This mirrors
% exactly how pumpkynPie's demos/lowThrustDRO2Tulip.m constructs them, so the
% direct and indirect formulations are solving the SAME problem -- which is the
% whole point of building the twin.
%
% THE ARRIVAL PHASE IS A CHOICE, not a given. The demo selects the point on the
% tulip whose velocity direction is most opposed to the DRO departure velocity
% (max angle). A different choice gives a different transfer. That criterion is
% reproduced here so the two methods are comparable; it is not claimed optimal.
%
% INPUTS:
%   opts - struct (optional):
%          .tauDRO  DRO period, ND                     [default 1.0]
%          .Np      tulip petal count                  [default 7]
%          .tauTul  tulip period, ND                   [default 5*2*pi/6]
%          .pm      tulip hemisphere (-1 southern)     [default -1]
%          .tol     cont_np tolerance                  [default 1e-12]
%
% OUTPUTS:
%   rv0 - departure state on the DRO [1x6, ND rotating barycentric]
%   rvf - arrival state on the tulip [1x6]
%   p   - struct: .muStar .lStar .tStar .tauDRO .tauTul .Np .idxArrival
%
% REFERENCES:
%   [1] pumpkynPie demos/lowThrustDRO2Tulip.m (the construction mirrored here).
%   [2] orbit_transfer/doc/pumpkyn_reference.md (getter + cont_np convention).

if nargin < 1, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
tauDRO = d('tauDRO', 1.0);
Np     = d('Np', 7);
tauTul = d('tauTul', 5*2*pi/6);
pm     = d('pm', -1);
tol    = d('tol', 1e-12);

p.muStar = 0.012150585609624;
p.lStar  = 389703.264829278;
p.tStar  = 382981.289129055;

% departure: DRO, closed by natural-parameter continuation
[~, rv0] = pumpkynPie.cr3bp.getDRO(tauDRO);
rv0 = pumpkyn.cr3bp.cont_np(rv0, tauDRO, p.muStar, tol);
[~, rvDRO] = pumpkyn.cr3bp.prop(tauDRO, rv0, p.muStar);

% arrival: tulip, then pick the phase most opposed to the departure velocity
[~, rvT] = pumpkyn.cr3bp.getTulip(tauTul, Np, pm);
rvT = pumpkyn.cr3bp.cont_np(rvT, tauTul, p.muStar, tol);
[~, rvTul] = pumpkyn.cr3bp.prop(tauTul, rvT, p.muStar);
ang = pumpkyn.util.bsxAng(rvTul(:,4:6), rvDRO(1,4:6), 2);
[~, idx] = max(ang);
rvf = rvTul(idx,:);

p.tauDRO = tauDRO;  p.tauTul = tauTul;  p.Np = Np;  p.idxArrival = idx;
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
