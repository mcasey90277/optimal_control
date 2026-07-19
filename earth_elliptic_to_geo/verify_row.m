function [pass, info] = verify_row(row, cert, tol)
% VERIFY_ROW  One-sided check of a reproduced Table-3 row against the campaign
% floor: the reproduced solution must be AT LEAST AS GOOD (final mass) as the
% campaign's certified number, within a small slack.
%
% The reproducer is a keep-best-mass OPTIMIZER (minimum-fuel = maximize final
% mass; user decision 2026-07-18, memory tenN-minfuel-razor-basin). So the gate
% is ONE-SIDED on mass: reproduced m_f_kg must be >= cert.m_f_kg - tol.m_f_kg.
% A HIGHER mass (less fuel) always passes and is flagged as an improvement --
% the reproducer is expected to equal or BEAT the campaign. Structure
% (switches, revs) is REPORTED in `info`, NOT gated: the best min-fuel optimum
% can have a different (better) bang-bang structure than the campaign row (10 N:
% the best optimum is 18 sw / 7.56 rev, vs the campaign's 19 sw / 7.326 rev).
%
% Pure function: no file I/O, no solves, no printing (the caller prints from
% `info`). Throws 'verify_row:worseThanCampaign' if the reproduced mass falls
% below the floor -- a regression is a loud failure, never a silent pass.
%
% INPUTS:
%   row  - reproduced row struct with fields .thrustN .m_f_kg .switches .revs [struct]
%   cert - campaign floor row (see table3_certified.m), same fields          [struct]
%   tol  - tolerance struct; only .m_f_kg is used (mass slack below the floor,
%          kg). Any .revsRel/.switchesAbs are ignored (structure not gated)   [struct]
%
% OUTPUTS:
%   pass - true if the reproduced mass clears the floor (throws otherwise)    [logical]
%   info - struct: .massFloor (cert.m_f_kg - tol.m_f_kg) .improvedKg
%          (row.m_f_kg - cert.m_f_kg, positive = beat the campaign) .improved
%          (logical) .switches .campaignSwitches .revs .campaignRevs          [struct]
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md
%       Sec. 5 + Sec. 9b (keep-best-mass, one-sided verify).
%   [2] memory tenN-minfuel-razor-basin (why structure is not gated).

if nargin < 3
    error('verify_row:badInput', 'row, cert, and tol are all required');
end

massFloor = cert.m_f_kg - tol.m_f_kg;
if row.m_f_kg < massFloor
    error('verify_row:worseThanCampaign', ...
        ['T=%g N: reproduced m_f=%.4f kg is BELOW the campaign floor ' ...
         '(%.4f - tol %.4f = %.4f kg) -- a worse min-fuel solution, refusing to pass'], ...
        cert.thrustN, row.m_f_kg, cert.m_f_kg, tol.m_f_kg, massFloor);
end

info = struct('massFloor', massFloor, ...
    'improvedKg', row.m_f_kg - cert.m_f_kg, ...
    'improved', row.m_f_kg > cert.m_f_kg + 1e-3, ...
    'switches', row.switches, 'campaignSwitches', cert.switches, ...
    'revs', row.revs, 'campaignRevs', cert.revs);
pass = true;

end
