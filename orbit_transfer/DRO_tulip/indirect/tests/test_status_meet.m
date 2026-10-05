function ok = test_status_meet()
% TEST_STATUS_MEET  status_meet is the meet of the optimality-status lattice
% ordered by what was ESTABLISHED: 1 < 3 < 2 and 1 < 3 < 4, 2 and 4
% incomparable. Full 4 x 4 table, symmetry, idempotence, and refusal of
% anything that is not an integer code 1..4.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
want = [1 1 1 1; ...          % meet(1, b) = 1
        1 2 3 3; ...          % meet(2, 3) = 3, meet(2, 4) = 3
        1 3 3 3; ...          % meet(3, 2) = meet(3, 4) = 3
        1 3 3 4];             % meet(4, 2) = meet(4, 3) = 3
got = zeros(4);
for a = 1:4
    for b = 1:4, got(a, b) = status_meet(a, b); end
end
ok = chk(ok, isequal(got, want), sprintf('the full 4 x 4 table (got %s)', mat2str(got)));
ok = chk(ok, isequal(got, got.'), 'symmetric: meet(a, b) = meet(b, a)');
ok = chk(ok, isequal(diag(got).', 1:4), 'idempotent: meet(x, x) = x');
ok = chk(ok, status_meet(int8(2), 4) == 3 && isa(status_meet(int8(2), 4), 'double'), 'integer classes accepted, double returned');
bad = {0, 5, 2.5, NaN, [], [2 3], '3', 3 + 1i};
ok = chk(ok, all(cellfun(@(v) throws(@() status_meet(v, 3)) && throws(@() status_meet(3, v)), bad)), ...
         'anything not an integer code 1..4 is refused (either argument)');
if ok, fprintf('test_status_meet: ALL PASS\n'); else, fprintf('test_status_meet: FAIL\n'); end
end

function t = throws(f)
% THROWS  True if f() throws.  INPUTS: f (handle).  OUTPUTS: t.
try
    f();  t = false;
catch
    t = true;
end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
