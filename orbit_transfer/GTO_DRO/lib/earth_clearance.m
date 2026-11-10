function [hMinKm, tAtMin, pass, out] = earth_clearance(tau, y, floorKm, p)
%% Purpose:
%
%   Altitude of a flown arc above a primary, in km, and whether it stays
%   above a floor. Pure geometry on the caller's samples: no propagation,
%   no interpolation between samples (the number is the SAMPLE minimum,
%   and says so). Default primary is the EARTH at [-muStar 0 0] with the
%   equatorial radius 6378.137 km -- the check the shared ladder engine
%   lacks for an Earth departure (it enforces a LUNAR floor only). With
%   p.center = [1-muStar 0 0], p.radiusKm = 1737.4 it is the lunar
%   altitude, so a study script measures both clearances with ONE
%   instrument.
%
%   Self-demo (nargin == 0): a circular 6678.137 km Earth orbit converted
%   with pumpkyn fromPCI at t = 0 reads 300 km.
%
%% Inputs:
%
%  tau                      [M x 1]                 sample times, ND
%
%  y                        [M x >= 6]              rotating-frame
%                                                   barycentric ND states,
%                                                   one row per time;
%                                                   columns 1:3 are used
%
%  floorKm                  double                  altitude floor, km
%
%  p                        struct                  .muStar, .lStar
%                                                   (required); .center
%                                                   [1 x 3] ND
%                                                   [-muStar 0 0]; .radiusKm
%                                                   [6378.137]
%
%% Outputs:
%
%  hMinKm                   double                  minimum altitude over the
%                                                   samples, km (negative
%                                                   when the arc crosses the
%                                                   surface)
%
%  tAtMin                   double                  tau at that sample
%
%  pass                     logical                 hMinKm >= floorKm
%
%  out                      struct                  .hKm [M x 1] altitude
%                                                   history, .center,
%                                                   .radiusKm, .iMin
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
    muE = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
    rC  = 6678.137;  vC = sqrt(muE/rC);  th = linspace(0, 2*pi, 64).';
    rvKm = [rC*cos(th), rC*sin(th), zeros(64,1), -vC*sin(th), vC*cos(th), zeros(64,1)];
    y = zeros(64, 6);
    for kk = 1:64
        y(kk,:) = pumpkyn.cr3bp.fromPCI(0, rvKm(kk,:), muStar, tStar, lStar, 1);
    end
    tau = linspace(0, 1, 64).';  floorKm = 200;
    p = struct('muStar', muStar, 'lStar', lStar);
end
if (nargin > 0 && nargin < 4) || ~isfield(p, 'muStar') || ~isfield(p, 'lStar')
    error('earth_clearance:params', 'p.muStar and p.lStar are required');
end
center   = fieldd(p, 'center',   [-p.muStar, 0, 0]);
radiusKm = fieldd(p, 'radiusKm', 6378.137);
tau = tau(:);
if size(y, 2) < 6 || size(y, 1) ~= numel(tau)
    error('earth_clearance:shape', 'y must be [numel(tau) x >= 6]; got [%d x %d] for %d times', ...
          size(y, 1), size(y, 2), numel(tau));
end

dKm  = p.lStar*sqrt(sum((y(:, 1:3) - center).^2, 2));     % distance to the centre, km
hKm  = dKm - radiusKm;
[hMinKm, iMin] = min(hKm);
tAtMin = tau(iMin);
pass   = hMinKm >= floorKm;
out    = struct('hKm', hKm, 'center', center, 'radiusKm', radiusKm, 'iMin', iMin);

if nargin == 0 || nargout == 0
    fprintf('earth_clearance: min altitude %.6f km at tau = %.4f (floor %.0f km) -> %s\n', ...
            hMinKm, tAtMin, floorKm, tern(pass, 'PASS', 'FAIL'));
end
end

function v = fieldd(s, f, d_)
% FIELDD  s.(f) if present else d_.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
