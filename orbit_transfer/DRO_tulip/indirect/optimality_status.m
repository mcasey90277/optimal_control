function [code, reason, inferred] = optimality_status(C, opts)
%% Purpose:
%
%   THE status of one certify_root result: 4 sufficient, 3 necessary only,
%   2 conjugate point found, 1 neither, -1 below the floor (not a transfer
%   of the library). One function, called by certify_root (stamping) and by
%   the harvester (legacy records), so the two cannot disagree.
%
%  ASSUMPTIONS / NOTES:
%
% • A result stamped by the current certifier carries .stage (the last
%   stage passed: 1 polish, 2 flight, 3 pointwise, 4 witness, 5 coarse
%   conjugate, 6 gates, 7 H6, 8 dense scan) and .conjFound/.hypAfterConj.
% • A LEGACY result has none of these; its status is inferred from the
%   reason text, conservatively: anything unrecognised takes the lowest
%   recordable tier and inferred = true. A legacy "conjugate test verdict
%   0" is FAIL or UNDETERMINED and is 3 until re-certified.
% • The floor: a finite z8 and a flown miss within opts.gateKm / gateVms.
%
%% Inputs:
%
%  C                        struct                  certify_root result
%  opts                     struct (optional)       .gateKm [100] .gateVms [10]
%
%% Outputs:
%
%  code                     double                  4 | 3 | 2 | 1 | -1
%  reason                   char                    why, in the certifier's words
%  inferred                 logical                 true for a legacy record
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
gateKm = fieldd(opts, 'gateKm', 100);  gateVms = fieldd(opts, 'gateVms', 10);
r = '';  if isfield(C, 'reason') && ischar(C.reason), r = C.reason; end
inferred = ~isfield(C, 'stage');

% ---- the floor --------------------------------------------------------------
z = [];  if isfield(C, 'z'), z = C.z; end
if ~(isnumeric(z) && numel(z) == 8 && all(isfinite(z(:))))
    code = -1;  reason = sprintf('below the floor: no converged root (%s)', r);  return
end
km = fieldd(C, 'flyKm', NaN);  vms = fieldd(C, 'flyVms', NaN);
if ~(isscalar(km) && isscalar(vms) && km <= gateKm && vms <= gateVms)
    code = -1;  reason = sprintf('below the floor: does not fly to the target (%s)', r);  return
end
if isfield(C, 'ok') && isequal(C.ok, true)
    code = 4;  reason = 'full stack passed';  return
end

if ~inferred
    % ---- stamped by the current certifier -----------------------------------
    if C.stage < 4
        code = 1;  reason = ['neither: necessary conditions not established -- ' r];
    elseif isfield(C, 'conjFound') && isequal(C.conjFound, true)
        if strcmp(C.hypAfterConj, 'held')
            code = 2;  reason = ['conjugate point found (gates and H6 held): ' r];
        else
            code = 3;  reason = sprintf(['necessary only: conjugate point found, but its hypotheses were not ' ...
                                         'established (%s) -- %s'], C.hypAfterConj, r);
        end
    else
        code = 3;  reason = ['necessary only: ' r];
    end
    return
end

% ---- legacy: inferred from the reason text -------------------------------------
pointwise = '^(Hamiltonian max|transversality|adjoint equations|minimum-principle|applied throttle|pointwise checks|the tight-tolerance flight|tfMin witness|the witness flight|witness flight inadmissible|witness solution)';
gates = '^(X2:|gates |gate |hypothesis gates|min\|lam_v\||min Q_mt|dim S|accepted lift|independent-field|second gates|lift_margin|H6|H2|H3|DIAGNOSTIC)';
if ~isempty(regexp(r, pointwise, 'once'))
    code = 1;  reason = ['neither (inferred): ' r];
elseif startsWith(r, 'conjugate test verdict')
    code = 3;  reason = ['necessary only (inferred): legacy verdict 0 is FAIL or UNDETERMINED; re-certify to resolve -- ' r];
elseif startsWith(r, 'dense conjugate scan not clear')
    n = str2double(regexp(r, '(\d+) coarse sign change\(s\), (\d+) zero, (\d+) UNRESOLVED, (\d+) multiplicity', 'tokens', 'once'));
    if numel(n) == 4 && (n(1) + n(2) + n(4)) > 0
        code = 2;  reason = ['conjugate point found (inferred; the dense scan runs after gates and H6): ' r];
    else
        code = 3;  reason = ['necessary only (inferred): ' r];
    end
elseif ~isempty(regexp(r, gates, 'once')) || contains(r, 'lower-bound ESTIMATE')
    code = 3;  reason = ['necessary only (inferred): ' r];
else
    code = 1;  reason = ['neither (inferred: reason not recognised): ' r];
end
end

% ---------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
