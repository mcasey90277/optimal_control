function [targets, sDpts] = rib_targets(sD, sD0, dirn)
% RIB_TARGETS  The unwrapped departure offsets a rib walks from its spine
% to reach every other phase of a departure list, in walking order.
%
%   A rib holds the arrival phase and steps the departure phase away from
%   the spine sD0 in one sense, so the offsets are ordered from the phase
%   nearest the spine to the farthest, and they are UNWRAPPED (a walk in the
%   -1 sense from sD0 = 0 reaches 0.9583 at offset -0.0417). The list may
%   have any spacing; it need not contain sD0 itself (if it does, that
%   entry is the spine and gets no target).
%
% INPUTS:
%   sD   - departure phases in [0,1), any order or spacing [1 x n]
%   sD0  - the spine's departure phase [scalar]
%   dirn - walking sense, -1 (decreasing phase) or +1 [scalar]
%
% OUTPUTS:
%   targets - unwrapped offsets from sD0, nearest first [1 x m], m = n or
%             n-1; negative for dirn = -1, positive for +1
%   sDpts   - the phases reached, mod 1, in the same order [1 x m]
%
% REFERENCES:
%   rib_from_crossing (opts.targets), FINDINGS 39 (why explicit targets beat
%   a derived step count).

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
