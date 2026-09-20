function ok = test_phase_branches()
% TEST_PHASE_BRANCHES  Cells joined by SAFE edges form a branch. On a 4 x 4
% torus: every edge safe is one branch; cutting the two column seams 2|3 and
% 4|1 leaves two; an empty cell belongs to none and does not bridge; a cell
% with no safe edge is a branch of its own; the wrap is an edge like any other.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
T = true(4);
[b, nb] = phase_branches(T, T, T);
ok = chk(ok, nb == 1 && all(b(:) == 1), 'every edge safe: one branch');
sA = T;  sA(:, 2) = false;  sA(:, 4) = false;                       % edges col 2->3 and col 4->1 are jumps
[b, nb] = phase_branches(T, sA, T);
ok = chk(ok, nb == 2 && all(all(b(:, 1:2) == b(1, 1))) && all(all(b(:, 3:4) == b(1, 3))) && b(1, 1) ~= b(1, 3), 'two column seams cut: columns {1,2} and {3,4}');
sA2 = T;  sA2(:, 2) = false;                                        % only ONE seam cut: the wrap still joins them
[~, nb] = phase_branches(T, sA2, T);
ok = chk(ok, nb == 1, 'one seam cut: the torus still connects the two halves the other way round');
has = T;  has(2, 2) = false;
[b, nb] = phase_branches(T, T, has);
ok = chk(ok, b(2, 2) == 0 && nb == 1, 'an empty cell belongs to no branch');
sD = T;  sAi = T;  sD(3, 3) = false;  sD(2, 3) = false;  sAi(3, 3) = false;  sAi(3, 2) = false;   % cell (3,3) cut off on all four sides
[b, nb] = phase_branches(sD, sAi, T);
ok = chk(ok, nb == 2 && nnz(b == b(3, 3)) == 1, 'a cell with no safe edge is a branch of its own');
if ok, fprintf('test_phase_branches: ALL PASS\n'); else, fprintf('test_phase_branches: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
