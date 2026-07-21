function fp = cr3bp_fingerprint(p, extra)
% CR3BP_FINGERPRINT  Build the config fingerprint that determines a solution.
%
% Captures the physics/config a cached artifact depends on, from the
% cr3bp_lt_params struct (so it cannot drift from the actual physics), plus
% caller-specific extras (tf, insertion, epsMin, ...). Consumed by
% check_cr3bp_fp at every cache read (2026-07-21 review triage C5/C6).
%
% INPUTS:
%   p     - cr3bp_lt_params struct (.thrustN .m0kg .ispS .Tmax .c .muStar) [struct]
%   extra - optional struct of run-specific fields to merge in [struct]
% OUTPUTS:
%   fp - fingerprint struct (.thrustN .m0kg .ispS .Tmax .cEx .muStar .pSund
%        + extras) [struct]
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-21-ladder-prep-design.md sec 2.
if nargin < 2 || isempty(extra), extra = struct(); end
fp = struct('thrustN', p.thrustN, 'm0kg', p.m0kg, 'ispS', p.ispS, ...
            'Tmax', p.Tmax, 'cEx', p.c, 'muStar', p.muStar);
if isfield(p, 'pSund'), fp.pSund = p.pSund; end
fn = fieldnames(extra);
for k = 1:numel(fn), fp.(fn{k}) = extra.(fn{k}); end
end
