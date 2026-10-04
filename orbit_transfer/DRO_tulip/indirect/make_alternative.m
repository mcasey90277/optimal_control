function A = make_alternative(C, sD, sA, iD, iA, source, opts)
%% Purpose:
%
%   One row of a catalog's ALTERNATIVES table from a certify_root result: the
%   transfer (z8, junctions, t_f), its optimality status and the numbers
%   behind it, where it sits and where it came from. [] when the result is
%   not a transfer of the library (below the floor) or came through the
%   certifier's test seam.
%
%  ASSUMPTIONS / NOTES:
%
% • A stamp (.status, .status_reason) is trusted only when the result is
%   stamped (optimality_status: a real finite .stage) AND .status is a
%   non-empty real finite scalar; otherwise the result is re-classified. A
%   legacy record whose missing fields were filled with [] is re-classified.
% • DISPLACED PRIMARIES (opts.displacedPrimary): a certified primary pushed
%   out of its cell by a faster root has no flight measured where it is
%   being moved (merge, packaging), so its flyKm/flyVms are NaN unless the
%   caller passes them. Such a row is BUILT DIRECTLY, bypassing the floor,
%   and only for a stamped status 4 with a real finite z8; any other result
%   on this path returns []. The status audit (audit_status_layer) re-flies
%   every alternative and is the safety net.
%
%% Inputs:
%
%  C                        struct                  certify_root result (stamped or legacy)
%  sD, sA                   double                  exact phases
%  iD, iA                   double                  grid cell, 0 = off the grid
%  source                   char                    provenance
%  opts                     struct (optional)       .displacedPrimary [false]
%
%% Outputs:
%
%  A                        struct or []            the row (fields: sD sA iD
%                                                   iA z8 junctions tf_nd
%                                                   status status_reason
%                                                   inferred conj conjVerdict
%                                                   minLamV minQmt dimS
%                                                   h6Margin liftMargin flyKm
%                                                   flyVms source)
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  M. Casey  final review C2/M2: a stamp is trusted only when .status is a
%            real scalar (else re-classified), guarded .status_reason;
%            displaced primaries built directly with a NaN flight 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 7, opts = struct(); end
displaced = isfield(opts, 'displacedPrimary') && isequal(opts.displacedPrimary, true);
A = [];
if isfield(C, 'overridden') && isequal(C.overridden, true), return, end
if displaced
    % a certified primary, displaced: no floor (its flight is not measured
    % here), but only for a stamped status 4 with a real finite z8
    [~, ~, inferred] = optimality_status(C);
    zOk = isfield(C, 'z') && isnumeric(C.z) && numel(C.z) == 8 && isreal(C.z) && all(isfinite(C.z(:)));
    if inferred || ~stampOk(C) || C.status ~= 4 || ~zOk, return, end
    code = 4;  why = reasonOf(C, code);
else
    [code, why, inferred] = optimality_status(C);      % classifies; also applies the floor to a stamped result
    if code < 1, return, end
    if ~inferred && stampOk(C)
        code = double(C.status);  why = reasonOf(C, code);
        if code < 1, A = []; return, end
    end
end
Y = [];  if isfield(C, 'Y'), Y = C.Y; end
g = struct('minLamV', NaN, 'minQmt', NaN, 'dimS', NaN);
if isfield(C, 'g') && isstruct(C.g) && ~isempty(C.g)
    for f = fieldnames(g)', if isfield(C.g, f{1}), g.(f{1}) = C.g.(f{1}); end, end
end
A = struct('sD', mod(sD, 1), 'sA', mod(sA, 1), 'iD', iD, 'iA', iA, 'z8', C.z(:), 'junctions', Y, ...
           'tf_nd', C.z(8), 'status', code, 'status_reason', why, 'inferred', inferred, ...
           'conj', fieldd(C, 'conj', NaN), 'conjVerdict', fieldd(C, 'conjVerdict', ''), ...
           'minLamV', g.minLamV, 'minQmt', g.minQmt, 'dimS', g.dimS, ...
           'h6Margin', fieldd(C, 'h6Margin', NaN), 'liftMargin', fieldd(C, 'liftMargin', NaN), ...
           'flyKm', fieldd(C, 'flyKm', NaN), 'flyVms', fieldd(C, 'flyVms', NaN), 'source', source);
end

% ---------------------------------------------------------------------------
function ok = stampOk(C)
% STAMPOK  .status is a non-empty real finite scalar.  INPUTS: C.  OUTPUTS: ok.
ok = isfield(C, 'status') && isnumeric(C.status) && isscalar(C.status) && isreal(C.status) && isfinite(C.status);
end

function why = reasonOf(C, code)
% REASONOF  The stamped reason, or a stand-in naming the stamped status when
% none was recorded.  INPUTS: C; code.  OUTPUTS: why char.
if isfield(C, 'status_reason') && ischar(C.status_reason) && ~isempty(C.status_reason)
    why = C.status_reason;
else
    why = sprintf('status %d (stamped; no reason recorded)', code);
end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
