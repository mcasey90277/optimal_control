function ok = test_lib()
%% Purpose:
%
%   The lib/ helpers against hand-computed answers. partial_f_test: a
%   six-point case worked by hand (RSS, F and the p-value from the beta
%   function), an exact-fit case, the :nested refusal when X0 is not in
%   the span of X1 and the :shape refusal.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'lib'));

%% partial_f_test: y = 1 + 2 x + 3 z exactly, X0 drops z -> huge F, p ~ 0
x = [0 1 2 3 4 5].';  z = [1 0 1 0 1 0].';
y = 1 + 2*x + 3*z;
o = partial_f_test(y, [ones(6, 1), x], [ones(6, 1), x, z]);
ok = check(ok, max(abs(o.b1 - [1; 2; 3])) < 1e-10 && o.RSS1 < 1e-18 && o.p < 1e-12, 'exact fit: coefficients recovered, RSS1 ~ 0, p ~ 0');
ok = check(ok, o.df1 == 1 && o.df2 == 3, 'degrees of freedom 1 and 3');

%% hand-worked case: y with noise, compare against the closed forms
y = [1.1; 2.9; 5.2; 6.8; 9.1; 11.0];
X0 = [ones(6, 1), x];  X1 = [ones(6, 1), x, z];
o  = partial_f_test(y, X0, X1);
b0 = X0 \ y;  RSS0 = sum((y - X0*b0).^2);
b1 = X1 \ y;  RSS1 = sum((y - X1*b1).^2);
F  = ((RSS0 - RSS1)/1) / (RSS1/3);
p  = betainc(3/(3 + F), 1.5, 0.5);
ok = check(ok, abs(o.F - F) < 1e-12 && abs(o.p - p) < 1e-12, sprintf('hand-worked F = %.4f, p = %.4f reproduced', F, p));
ok = check(ok, o.p > 0 && o.p < 1 && o.R2_1 >= o.R2_0, 'p in (0, 1) and R2 non-decreasing with the larger model');
ok = check(ok, abs(o.pt1(3) - o.p) < 1e-12 && abs(o.t1(3)^2 - o.F) < 1e-9, 'one added column: its t-test p equals the partial-F p and t^2 = F');
ok = check(ok, isequal(size(o.se1), [3 1]) && all(o.se1 > 0), 'standard errors are a positive [p1 x 1] column');
ok = check(ok, isequal(size(o.ci1), [3 2]) && all(o.ci1(:, 1) < o.b1) && all(o.ci1(:, 2) > o.b1), '95% CI brackets each coefficient');
tq = sqrt(3*(1/betaincinv(0.05, 1.5, 0.5) - 1));
ok = check(ok, abs(tq - 3.1824) < 2e-4 && abs(o.ci1(3, 2) - o.ci1(3, 1) - 2*tq*o.se1(3)) < 1e-12, 'CI half-width is t_{0.975,3} = 3.182 standard errors');
%% exact-fit policies and the rank refusal
o2 = partial_f_test(y, [ones(6, 1), x, z, [1; 0; 0; 1; 0; 1]], [ones(6, 1), x, z, [1; 0; 0; 1; 0; 1], [0; 1; 0; 0; 1; 0]]);
ok = check(ok, o2.F >= 0 && o2.p >= 0 && o2.p <= 1, 'four and five columns on six points: finite F, p in [0, 1]');
yE = 1 + 2*x;
o3 = partial_f_test(yE, [ones(6, 1), x], [ones(6, 1), x, z]);
ok = check(ok, o3.F == 0 && o3.p == 1, 'both models exact (RSS0 = RSS1 = 0): no improvement claimed, p = 1');
try
    partial_f_test(y, [ones(6, 1), x], [ones(6, 1), x, 2*x]);
    ok = check(ok, false, ':rank refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:partial_f_test:rank'), ':rank refusal on a rank-deficient X1');
end

%% the p-value against a tabulated point: F(1, 3) upper tail at F = 10.13 is 0.05
o = struct('F', 10.128, 'df1', 1, 'df2', 3);
ok = check(ok, abs(betainc(o.df2/(o.df2 + o.df1*o.F), o.df2/2, o.df1/2) - 0.05) < 2e-4, 'beta-function tail matches the F(1,3) 5% critical value 10.13');

%% refusals
try
    partial_f_test(y, [ones(6, 1), z], [ones(6, 1), x]);
    ok = check(ok, false, ':nested refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:partial_f_test:nested'), ':nested refusal when X0 is not in the span of X1');
end
try
    partial_f_test(y, X0(1:5, :), X1);
    ok = check(ok, false, ':shape refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:partial_f_test:shape'), ':shape refusal on a row mismatch');
end

%% mirror_state: the hand vector, a matrix, and the refusal
x  = [1; 2; 3; 4; 5; 6];
ok = check(ok, isequal(mirror_state(x), [1; -2; 3; -4; 5; -6]), 'mirror_state flips y, vx, vz');
ok = check(ok, isequal(mirror_state([x, 2*x]), [[1; -2; 3; -4; 5; -6], 2*[1; -2; 3; -4; 5; -6]]), 'mirror_state acts column-wise on [6 x n]');
ok = check(ok, isequal(mirror_state(mirror_state(x)), x), 'mirror_state is an involution');
try
    mirror_state(x(1:5));
    ok = check(ok, false, 'mirror_state :shape refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:mirror_state:shape'), 'mirror_state :shape refusal on a 5-vector');
end

%% symmetry_residual: at an OFF-PLANE point a symmetric orbit maps onto itself; a non-periodic state does not
addpath(fullfile(fileparts(fileparts(here)), 'costate_common'));
[~, ~, rv, info] = evalc('get_family_orbit(''dro'', struct(''tau'', 1))');
x0 = rv(1, :).';  Tp = info.periodND;  mu = 0.012150585609624;
ok = check(ok, max(abs(x0([2 4 6]))) < 1e-12, 'the stored initial state sits on the symmetry plane (so a check there would be vacuous)');
[~, X] = pumpkyn.cr3bp.prop([0 0.37*Tp], x0.', mu);
x1 = X(end, 1:6).';
ok = check(ok, max(abs(x1([2 4 6]))) > 1e-3, 'phi(x0, 0.37 T) is off the plane');
[res, tS] = symmetry_residual(x1, Tp, mu);
ok = check(ok, res < 1e-8, sprintf('DRO tau 1: mirror image of the off-plane point lies on the orbit (residual %.1e at t = %.4f of T = %.4f)', res, tS, Tp));
ok = check(ok, tS >= 0 && tS <= Tp, 'minimizing time inside [0, T]');
xn = x1;  xn(5) = xn(5)*1.02;                              % 2% faster: not on the periodic orbit
resN = symmetry_residual(xn, Tp, mu);
ok = check(ok, resN > 1e-6, sprintf('a non-periodic state is NOT mapped onto an orbit (residual %.1e): the gate can fail', resN));

fprintf('test_lib: %s\n', verdict(ok));
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
function s = verdict(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
