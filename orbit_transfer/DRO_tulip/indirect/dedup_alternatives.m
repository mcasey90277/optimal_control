function A = dedup_alternatives(A, sheet)
%% Purpose:
%
%   One root, one record: drop every alternative that is the same root
%   (same_root) as the PRIMARY of its cell, or as an earlier row at the same
%   phases; return the rest sorted by (sD, sA, t_f).
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
            keep(k) = false;  break
        end
    end
end
A = A(keep);
[~, order] = sortrows([[A.sD].', [A.sA].', [A.tf_nd].']);
A = A(order);
end
