function pass = verify_row(row, cert, tol)
% VERIFY_ROW  Assert a (reproduced or candidate) Table-3 row against the
% certified numbers, within tolerance.
%
% Pure function: no file I/O, no solves, no side effects. Compares
% row.m_f_kg vs cert.m_f_kg (absolute tolerance tol.m_f_kg), row.revs vs
% cert.revs (relative tolerance tol.revsRel), and row.switches vs
% cert.switches (absolute tolerance tol.switchesAbs). On any breach,
% throws 'verify_row:mismatch' naming the offending field and its expected
% and actual values; a mismatch is a loud failure, never a silent pass.
%
% INPUTS:
%   row  - candidate row struct with fields .m_f_kg, .revs, .switches     [struct]
%   cert - certified row struct (see table3_certified.m) with the same
%          three fields                                                   [struct]
%   tol  - tolerance struct with fields:
%            .m_f_kg      - absolute tolerance on m_f_kg [kg]              [scalar]
%            .revsRel     - relative tolerance on revs [fraction]          [scalar]
%            .switchesAbs - absolute tolerance on switches (integer count) [scalar]
%
% OUTPUTS:
%   pass - true if row passes all three checks (function throws instead
%          of returning false on any breach)                             [logical]
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md Sec. 5.

if nargin < 3
    error('verify_row:badInput', 'row, cert, and tol are all required');
end

dm = abs(row.m_f_kg - cert.m_f_kg);
if dm > tol.m_f_kg
    error('verify_row:mismatch', ...
        'm_f_kg mismatch: expected %.4f kg (tol %.4f), got %.4f kg (|diff|=%.4f)', ...
        cert.m_f_kg, tol.m_f_kg, row.m_f_kg, dm);
end

drevRel = abs(row.revs - cert.revs) / abs(cert.revs);
if drevRel > tol.revsRel
    error('verify_row:mismatch', ...
        'revs mismatch: expected %.4f (relTol %.4f), got %.4f (relDiff=%.4f)', ...
        cert.revs, tol.revsRel, row.revs, drevRel);
end

dsw = abs(row.switches - cert.switches);
if dsw > tol.switchesAbs
    error('verify_row:mismatch', ...
        'switches mismatch: expected %d (tol %d), got %d (|diff|=%d)', ...
        cert.switches, tol.switchesAbs, row.switches, dsw);
end

pass = true;

end
