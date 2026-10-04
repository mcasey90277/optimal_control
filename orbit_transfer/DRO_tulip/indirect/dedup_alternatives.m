function A = dedup_alternatives(A, sheet)
%% Purpose:
%
%   One root, one record: drop every alternative that is the same root
%   (same_root) as the PRIMARY of its cell; within a group of rows that are
%   the same root at the same phases keep ONE -- a stamped row over an
%   inferred one, then the higher established status, then the first.
%   Return the rest sorted by (sD, sA, t_f).
%
%% Inputs:
%
%  A                        struct array            alternatives rows (make_alternative)
%  sheet                    struct                  catalog sheet (.has_solution
%                                                   .entry_index .z8)
%
%% Outputs:
%
%  A                        struct array            deduplicated, sorted
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  M. Casey  final review I2: the kept twin is the stamped, then the
%            higher-status, then the first row                       10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(A), return, end
keep = true(1, numel(A));
near = @(a, b) abs(mod(a - b + 0.5, 1) - 0.5) < 1e-9;
for k = 1:numel(A)
    a = A(k);
    if a.iD > 0 && a.iA > 0 && sheet.has_solution(a.iD, a.iA, 1)
        if same_root(a.z8, sheet.z8(:, sheet.entry_index(a.iD, a.iA, 1))), keep(k) = false;  continue, end
    end
    for m = 1:k-1
        if keep(m) && near(A(m).sD, a.sD) && near(A(m).sA, a.sA) && same_root(A(m).z8, a.z8)
            % one kept row per group: the later row replaces it only if better
            if better(a, A(m)), keep(m) = false; else, keep(k) = false; end
            break
        end
    end
end
A = A(keep);
[~, order] = sortrows([[A.sD].', [A.sA].', [A.tf_nd].']);
A = A(order);
end

% ---------------------------------------------------------------------------
function tf = better(a, b)
% BETTER  Row a is preferred to row b: stamped over inferred, then the higher
% status; a tie is not better (the first is kept).  INPUTS: a; b (rows).
% OUTPUTS: tf logical.
ia = isequal(a.inferred, true);  ib = isequal(b.inferred, true);
if ia ~= ib, tf = ~ia;  return, end
tf = double(a.status) > double(b.status);
end
