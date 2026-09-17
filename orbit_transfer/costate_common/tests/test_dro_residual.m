function test_dro_residual()
% TEST_DRO_RESIDUAL  Self-test for the accuracy gate's measuring instrument.
%
% dro_residual is now the number that decides whether a direct solution counts
% as an answer, so it needs its own check. A measuring instrument that reports
% large errors is useless if it cannot report a SMALL one on a trajectory known
% to be accurate -- that is the failure mode this test exists to catch.
%
% THREE TESTS
%   1. NULL TEST. Build a trajectory by high-accuracy integration, sample it on
%      a grid, and feed it in. The residual must come back at integrator
%      tolerance (~1e-10), not at scheme tolerance. If it does not, the engine's
%      dynamics disagree with the propagator's and every number it has ever
%      reported is wrong.
%   2. CONTROL RECONSTRUCTION. The same trajectory presented as a
%      Hermite-Simpson solution (with midpoint controls) must not be made WORSE
%      by the quadratic reconstruction than by the linear one. This catches a
%      mis-ordered or mis-scaled Lagrange basis.
%   3. SENSITIVITY. Deliberately corrupt one node by a known amount; the
%      residual on the two intervals touching it must rise by about that amount.
%      A gate that cannot see an injected error cannot see a real one.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per test; errors on failure)
%
% REFERENCES:
%   [1] orbit_transfer/DRO_tulip/direct/certify/dro_residual.m

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));

mu   = 0.012150585609624;
Tmax = 0.1756418;
c    = 8.673746;

% a benign arc, well away from the Moon, so the test measures the ENGINE and
% not the stiffness of a close approach
x0  = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 1.0];
uFix = [0.6; -0.7; 0.39; 1.0];  uFix(1:3) = uFix(1:3)/norm(uFix(1:3));
tf  = 0.30;  N = 60;
tN  = linspace(0, tf, N+1);

odeo = odeset('RelTol',3e-13,'AbsTol',1e-15);
[~, Z] = ode113(@(t,z) rhs(z, uFix, mu, Tmax, c), tN, x0, odeo);

o = struct('X', Z.', 'U', repmat(uFix,1,N+1), 's', tN/tf, 'tf', tf);
nFail = 0;

%% 1. null test -- a truly accurate trajectory must read as accurate
R1 = dro_residual(o, mu, Tmax, c);
tol1 = 1e-9;
ok = R1.RxMax < tol1;
nFail = nFail + report('null test (exact trajectory reads as exact)', R1.RxMax, tol1, ok);

%% 2. quadratic control reconstruction is no worse than linear
% The control here is constant, so the quadratic through three equal values must
% reproduce it exactly and give the same residual as the linear path.
o2 = o;  o2.scheme = 'hermite-simpson';  o2.Um = repmat(uFix,1,N);
R2 = dro_residual(o2, mu, Tmax, c);
ok = abs(R2.RxMax - R1.RxMax) < max(1e-12, 0.05*R1.RxMax);
nFail = nFail + report('quadratic reconstruction matches linear on a constant control', ...
                        abs(R2.RxMax - R1.RxMax), max(1e-12, 0.05*R1.RxMax), ok);

%% 3. sensitivity -- an injected error must be seen, at about its own size
kBad = 30;  bump = 1e-5;
o3 = o;  o3.X(1,kBad) = o3.X(1,kBad) + bump;
R3 = dro_residual(o3, mu, Tmax, c);
seen = max(R3.Rx([kBad-1 kBad]));
ok = seen > 0.3*bump && seen < 30*bump;
nFail = nFail + report(sprintf('sees an injected %.0e position error', bump), seen, bump, ok);

fprintf('\n');
if nFail > 0
    error('test_dro_residual:fail', '%d of 3 tests FAILED', nFail);
end
fprintf('  test_dro_residual: 3/3 PASS\n');
end

% ---------------------------------------------------------------------------
function bad = report(name, val, tol, ok)
% REPORT  One line per test. INPUTS: name; val; tol; ok  OUTPUTS: bad (0 or 1)
if ok, mark = 'PASS'; bad = 0; else, mark = 'FAIL'; bad = 1; end
fprintf('  %-58s %10.3e (ref %8.1e)  %s\n', name, val, tol, mark);
end

% ---------------------------------------------------------------------------
function dz = rhs(z, u, mu, Tmax, c)
% RHS  CR3BP with thrust; must mirror dro_residual/local_rhs exactly -- the null
% test is precisely the check that it does.
% INPUTS: z [7x1]; u [4x1]; mu; Tmax; c   OUTPUTS: dz [7x1]
r = z(1:3);  v = z(4:6);  m = z(7);
al = u(1:3);  al = al/max(norm(al),eps);  th = min(max(u(4),0),1);
dd = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
gr = [r(1) - (1-mu)*(r(1)+mu)/dd^3 - mu*(r(1)-1+mu)/rr^3;
      r(2) - (1-mu)*r(2)/dd^3      - mu*r(2)/rr^3;
           - (1-mu)*r(3)/dd^3      - mu*r(3)/rr^3];
dz = [v; gr + [2*v(2); -2*v(1); 0] + (th*Tmax/m)*al; -(Tmax/c)*th];
end
