function ok = test_direct_ref()
%% Purpose:
%
%   The committed direct-solution fixture is a real solution of the problem
%   the indirect demo solves: right shape, boundary conditions met, its own
%   trapezoidal defects small, and a control that stayed inside the relaxed
%   bound (so it is the UNCONSTRAINED optimum the PMP solve is compared to).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));

ok = chk(ok, isequal(size(R.X), [4 201]) && isequal(size(R.U), [1 201]) ...
             && isequal(size(R.muDefect), [4 200]) && numel(R.tN) == 201, ...
         'shapes: X [4 x 201], U [1 x 201], muDefect [4 x 200]');
ok = chk(ok, R.tN(1) == 0 && abs(R.tN(end) - 5) < 1e-12 && abs(R.tf - 5) < 1e-12, ...
         'the grid spans the fixed 5 s horizon');
ok = chk(ok, max(abs(R.X(:,1) - [0;0;0;0])) < 1e-8, 'starts at rest, pendulum down');
ok = chk(ok, max(abs(R.X(:,end) - [0;pi;0;0])) < 1e-6, ...
         sprintf('ends at rest, pendulum up (miss %.1e)', max(abs(R.X(:,end) - [0;pi;0;0]))));
ok = chk(ok, max(abs(R.U)) < 0.9*R.uMax, ...
         sprintf('the control stayed off the relaxed bound: max |u| = %.1f N of %.0f N', ...
                 max(abs(R.U)), R.uMax));

% its own trapezoidal defects, recomputed here from the stored nodes
h = diff(R.tN);
F = zeros(4, 201);
for k = 1:201
    F(:,k) = ref_field(R.X(:,k), R.U(k), R.p);
end
d = R.X(:,2:end) - R.X(:,1:end-1) - (h/2).*(F(:,2:end) + F(:,1:end-1));
ok = chk(ok, max(abs(d(:))) < 1e-6, sprintf('trapezoidal defects small: %.1e', max(abs(d(:)))));
ok = chk(ok, R.J > 0 && isfinite(R.J), sprintf('a finite positive cost: J = %.6f', R.J));

if ok, fprintf('TEST_DIRECT_REF: ALL PASS\n');
else,  fprintf('TEST_DIRECT_REF: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function dx = ref_field(x, u, p)
%% Purpose:
%
%   The example's own dynamics, written out here so the fixture is checked
%   against the problem statement rather than against the code under test.
%
q2 = x(2);  q1d = x(3);  q2d = x(4);
s = sin(q2);  c = cos(q2);
D1 = p.m1 + p.m2*(1 - c^2);
D2 = p.L*(p.m1 + p.m2)*(1 - (p.m2/(p.m1 + p.m2))*c^2);
dx = [q1d; q2d;
      (p.L*p.m2*s*q2d^2 + u + p.m2*p.g*c*s)/D1;
      (p.L*p.m2*c*s*q2d^2 + u*c + (p.m1 + p.m2)*p.g*s)/D2];
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
