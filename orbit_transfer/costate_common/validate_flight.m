function V = validate_flight(t, Y, tf, Tnd, cnd, mu, lStar, opts)
%% Purpose:
%
%   ONE admissibility check for a flown all-burn trajectory, shared by the
%   certifier and the study script so that "the flight is admissible" means
%   the same thing everywhere. A returned array is not a completed flight:
%   an integrator that stops early without throwing hands back a short,
%   perfectly finite trajectory, and every metric taken from its last row
%   then describes a flight that never happened. This checks
%
%     1. the flight reached t_f (last sample within relTf of t_f);
%     2. every state sample is finite;
%     3. the all-burn mass law m(t) = 1 - (T/c) t holds at the end and the
%        final mass is positive;
%     4. closest approach to the Moon and to the Earth exceeds the stated
%        clearances (a trajectory through a primary is not a solution).
%
%   V.reason names the FIRST failed check.
%
%% Inputs:
%
%  t                        [N x 1]                 sample times [ND]
%  Y                        [N x >=7]               states, columns 1:3 r,
%                                                   4:6 v, 7 m
%  tf                       scalar                  intended final time [ND]
%  Tnd, cnd, mu             scalars                 thrust, exhaust speed,
%                                                   mass ratio [ND]
%  lStar                    scalar                  length unit [km]
%  opts                     struct (optional)
%   .relTf [1e-8] .relMass [1e-6] .moonKmMin [1900] .earthKmMin [6600]
%
%% Outputs:
%
%  V                        struct                  .ok .reason .tEnd
%                                                   .mEnd .mExpect .moonKm
%                                                   .earthKm (closest
%                                                   approaches)
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 8, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
relTf = d('relTf', 1e-8);  relMass = d('relMass', 1e-6);
moonKmMin = d('moonKmMin', 1900);  earthKmMin = d('earthKmMin', 6600);

V = struct('ok', false, 'reason', '', 'tEnd', NaN, 'mEnd', NaN, ...
           'mExpect', NaN, 'moonKm', NaN, 'earthKm', NaN);
if isempty(t) || isempty(Y) || size(Y, 2) < 7
    V.reason = 'flight returned no usable samples';  return
end
[okT, tEnd] = scalar_verdict(t(end));
V.tEnd = tEnd;
if ~okT || ~(tf > 0) || abs(tEnd - tf) > relTf*max(tf, 1)
    V.reason = sprintf('flight did not reach t_f (t_end = %.6g vs t_f = %.6g)', tEnd, tf);
    return
end
if ~all(isfinite(Y(:, 1:7)), 'all')
    V.reason = 'flight returned non-finite states';  return
end
V.mEnd = Y(end, 7);  V.mExpect = 1 - (Tnd/cnd)*tf;
if ~(V.mEnd > 0) || abs(V.mEnd - V.mExpect) > relMass*max(V.mExpect, 1)
    V.reason = sprintf('mass law violated: m(t_f) = %.9g, all-burn expects %.9g', ...
                       V.mEnd, V.mExpect);
    return
end
% closest approaches on the SAMPLES: a coarse grid can miss a graze, so
% the clearances are a screen against gross penetration, not a proof
rMoon  = sqrt(sum((Y(:, 1:3) - [1 - mu, 0, 0]).^2, 2));
rEarth = sqrt(sum((Y(:, 1:3) - [-mu, 0, 0]).^2, 2));
V.moonKm  = min(rMoon)*lStar;
V.earthKm = min(rEarth)*lStar;
if V.moonKm < moonKmMin
    V.reason = sprintf('closest approach to the Moon %.0f km < %g km', V.moonKm, moonKmMin);
    return
end
if V.earthKm < earthKmMin
    V.reason = sprintf('closest approach to the Earth %.0f km < %g km', V.earthKm, earthKmMin);
    return
end
V.ok = true;  V.reason = 'admissible';
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
