function [rv0, rvf, elfoTrace] = gto_elfo_endpoints(p, opts)
% GTO_ELFO_ENDPOINTS  Boundary states for the GTO -> lunar-ELFO transfer.
%
% Same departure as gto_tulip_endpoints (GTO 350 x 35786 km, argp -25 deg) but
% the arrival state is a rendezvous point on an elliptical lunar frozen orbit
% (ELFO), constructed exactly as proj7's im_elfo_states.m: build the orbit about
% the Moon with fromOrb, propagate under CR3BP into the ND barycentric rotating
% frame, and pick the rendezvous point. Default = APOLUNE (highest, slowest,
% over the south pole -- the natural, gentle low-thrust capture point).
%
% INPUTS:
%   p    - parameter struct from CR3BP_LT_PARAMS (uses .muStar,.lStar,.tStar)
%   opts - (optional) struct:
%          .oe   - [sma_km ecc inc_deg argp_deg]  (default [12000 0.69 56.5 90])
%          .raan - RAAN of the ELFO member, deg    (default 0)
%          .M0   - mean anomaly seed, deg          (default 0)
%          .point- 'apolune' | 'perilune' | 'nearest' (default 'apolune')
%                  'nearest' = ELFO state closest (full 6D) to opts.ref; use
%                  for orbit INSERTION (arrive on the ELFO at the easiest phase,
%                  avoiding a large terminal-velocity swing from the transfer).
%          .ref  - reference state [1x6] for 'nearest' (required if point=nearest)
%
% OUTPUTS:
%   rv0       - departure GTO state [1x6] (ND, rotating frame) [r v]
%   rvf       - ELFO rendezvous state [1x6] (ND, rotating frame) [r v]
%   elfoTrace - one-period ELFO trace [Kx6] for plotting (optional)
%
% REFERENCES:
%   [1] pumpkyn.cr3bp (Koblick) - orb2eci, fromPCI, fromOrb, prop.
%   [2] proj7/pipelines/gdop/im_elfo_states.m, im_elfo_optimum.m (ELFO shape).

if nargin < 2, opts = struct(); end
gd = @(f,d) subsref_default(opts,f,d);
oe    = gd('oe',   [12000, 0.69, 56.5, 90]);
raan  = gd('raan', 0);
M0    = gd('M0',   0);
point = gd('point','apolune');

% --- departure: identical GTO to the tulip problem --------------------------
muEarth = 6.67384e-20*(1-p.muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;   ecc = (35786-350)/(2*sma);
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0, v0], p.muStar, p.tStar, p.lStar, 1);

% --- arrival: ELFO rendezvous point in the ND barycentric rotating frame ----
Mmass = 5.9736e24 + 7.35e22;
nu0   = pumpkyn.util.mean2True(M0*pi/180, oe(2));
oev   = [oe(1), oe(2), oe(3)*pi/180, oe(4)*pi/180, raan*pi/180, nu0];
x0    = pumpkyn.cr3bp.fromOrb(0, oev, Mmass, p.muStar, p.tStar, p.lStar, 2);

tau   = linspace(0, 0.5, 8000)';           % >~ one ELFO period (~0.31 ND)
[~, xs] = pumpkyn.cr3bp.prop(tau, x0, p.muStar);
dMoon = sqrt(sum((xs(:,1:3) - [1-p.muStar, 0, 0]).^2, 2));
switch lower(point)
    case 'perilune', [~, idx] = min(dMoon);
    case 'nearest'
        ref = gd('ref', []);
        assert(~isempty(ref), 'gto_elfo_endpoints:noRef', 'point=nearest requires opts.ref');
        d6 = sqrt(sum((xs(:,1:6) - ref(:).').^2, 2));
        [~, idx] = min(d6);
    otherwise, [~, idx] = max(dMoon);        % apolune
end
rvf = xs(idx, 1:6);
elfoTrace = xs(:, 1:6);
end

% ---------------------------------------------------------------------------
function v = subsref_default(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
