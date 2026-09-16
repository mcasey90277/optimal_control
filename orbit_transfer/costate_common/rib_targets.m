function [targets, sDpts] = rib_targets(sD, sD0, dirn)
%% Purpose:
%
%   The unwrapped departure offsets a rib walks from its spine to reach every
%   other phase of a departure list, in walking order.
%
%   A rib holds the arrival phase and steps the departure phase away from the
%   spine sD0 in one sense, so the offsets are ordered from the phase nearest
%   the spine to the farthest, and they are UNWRAPPED (a walk in the -1 sense
%   from sD0 = 0 reaches 0.9583 at offset -0.0417). The list may have any
%   spacing; it need not contain sD0 itself (if it does, that entry is the
%   spine and gets no target).
%
%% References:
%   [1] rib_from_crossing (opts.targets), FINDINGS 39 (why explicit targets
%       beat a derived step count).
%
%% Inputs:
%
%  sD                       [1 x n]                 Departure phases in [0,1),
%                                                   any order or spacing
%
%  sD0                      double                  The spine's departure
%                                                   phase
%
%  dirn                     double                  Walking sense, -1
%                                                   (decreasing phase) or +1
%
%% Outputs:
%
%  targets                  [1 x m]                 Unwrapped offsets from
%                                                   sD0, nearest first, m = n
%                                                   or n-1; negative for dirn
%                                                   = -1, positive for +1
%
%  sDpts                    [1 x m]                 The phases reached, mod 1,
%                                                   in the same order
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

sD = mod(sD(:).', 1);
assert(isscalar(sD0) && isscalar(dirn) && (dirn == -1 || dirn == 1), ...
       'rib_targets: sD0 scalar, dirn -1 or +1');
sDpts = sD(abs(mod(sD - sD0 + 0.5, 1) - 0.5) > 1e-9);     % the spine itself is not a target
if dirn < 0
    off = mod(sD0 - sDpts, 1);                             % distance walked downward
else
    off = mod(sDpts - sD0, 1);                             % distance walked upward
end
[off, order] = sort(off);
sDpts = sDpts(order);
targets = dirn*off;
end
