function [token, why] = backfill_clean_verdict(V)
%% Purpose:
%
%   The ONE decision token of the status-layer backfill (backfill_v2_job,
%   verdict stage): 'CLEAN' only if every condition for adopting the v2
%   catalog holds, 'NOT CLEAN' otherwise, with the failing conditions named.
%   A pure function, so the rule is tested without a run.
%
%  ASSUMPTIONS / NOTES:
%
% • CLEAN requires ALL of: the audit's nBad == 0 (a real finite scalar);
%   coverage complete; primNot4 and movedPrimaries empty; the content key
%   equal to the record's; and the comparison with the record showing no
%   change -- every shared cell agrees (nAgree == nBoth), none missing
%   (nOnlyRef == 0), none slower, none faster.
% • Anything absent or malformed counts as failing (fail closed).
%
%% Inputs:
%
%  V                        struct                  .nBad [scalar] .coverageOk
%                                                   [logical] .primNot4 [1 x n]
%                                                   .movedPrimaries [1 x n]
%                                                   .contentKeyEqual [logical]
%                                                   .cmp (compare_phase_catalogs:
%                                                   .nAgree .nBoth .nOnlyRef
%                                                   .nNewSlower .nNewFaster)
%
%% Outputs:
%
%  token                    char                    'CLEAN' | 'NOT CLEAN'
%  why                      char                    the failing conditions, ''
%                                                   when CLEAN
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

fails = {};
if ~(isfield(V, 'nBad') && isZero(V.nBad)), fails{end+1} = 'audit nBad ~= 0'; end
if ~(isfield(V, 'coverageOk') && isequal(V.coverageOk, true)), fails{end+1} = 'audit coverage incomplete'; end
if ~(isfield(V, 'primNot4') && isempty(V.primNot4)), fails{end+1} = 'primNot4 not empty'; end
if ~(isfield(V, 'movedPrimaries') && isempty(V.movedPrimaries)), fails{end+1} = 'movedPrimaries not empty'; end
if ~(isfield(V, 'contentKeyEqual') && isequal(V.contentKeyEqual, true)), fails{end+1} = 'content key differs from the record'; end
if ~isfield(V, 'cmp') || ~isstruct(V.cmp)
    fails{end+1} = 'no comparison with the record';
else
    C = V.cmp;
    if ~(isfield(C, 'nAgree') && isfield(C, 'nBoth') && isScalarNum(C.nAgree) && isequal(C.nAgree, C.nBoth))
        fails{end+1} = 'not every shared cell agrees with the record';
    end
    for f = {'nOnlyRef', 'nNewSlower', 'nNewFaster'}
        if ~(isfield(C, f{1}) && isZero(C.(f{1}))), fails{end+1} = sprintf('%s ~= 0', f{1}); end
    end
end
if isempty(fails), token = 'CLEAN';  why = '';
else, token = 'NOT CLEAN';  why = strjoin(fails, '; ');
end
end

% ---------------------------------------------------------------------------
function tf = isZero(x)
% ISZERO  A real finite scalar equal to 0.  INPUTS: x.  OUTPUTS: tf.
tf = isScalarNum(x) && x == 0;
end

function tf = isScalarNum(x)
% ISSCALARNUM  A real finite numeric scalar.  INPUTS: x.  OUTPUTS: tf.
tf = isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x);
end
