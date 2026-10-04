function A = make_alternative(C, sD, sA, iD, iA, source)
%% Purpose:
%
%   One row of a catalog's ALTERNATIVES table from a certify_root result: the
%   transfer (z8, junctions, t_f), its optimality status and the numbers
%   behind it, where it sits and where it came from. [] when the result is
%   not a transfer of the library (below the floor) or came through the
%   certifier's test seam.
%
%% Inputs:
%
%  C                        struct                  certify_root result (stamped or legacy)
%  sD, sA                   double                  exact phases
%  iD, iA                   double                  grid cell, 0 = off the grid
%  source                   char                    provenance
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
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

A = [];
if isfield(C, 'overridden') && isequal(C.overridden, true), return, end
[code, why, inferred] = optimality_status(C);          % classifies; also applies the floor to a stamped result
if code < 1, return, end
if isfield(C, 'status') && isfield(C, 'stage')
    code = C.status;  why = C.status_reason;  inferred = false;
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
           'flyKm', C.flyKm, 'flyVms', C.flyVms, 'source', source);
end

% ---------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
