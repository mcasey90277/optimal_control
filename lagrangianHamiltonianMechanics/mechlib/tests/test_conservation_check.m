function ok = test_conservation_check()
%% Purpose:
%
%   mech.conservation_check on exact circular motion of a free particle
%   held on a circle (passes), an injected energy ramp (fails with the
%   energy reason), an injected constraint violation (fails with the
%   constraint reason), and the transposed-input refusal.
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
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

syms x y xd yd real
sys = mech.lagrange_dalembert([x; y], [xd; yd], (xd^2 + yd^2)/2, [x, y]);
t = linspace(0, 4*pi, 401);
Y = [cos(t); sin(t); -sin(t); cos(t)];               % exact circle, |v| = 1

chk = mech.conservation_check(t, Y, sys);
ok = check(ok, chk.ok && chk.driftRel < 1e-12 && chk.resid < 1e-12, ...
           sprintf('exact circle: drift %.1e, residual %.1e, ok = %d', chk.driftRel, chk.resid, chk.ok));
ok = check(ok, isequal(size(chk.E), [1 401]) && max(abs(chk.E - 0.5)) < 1e-12, 'E is [1 x N] and equals 1/2');
ok = check(ok, isempty(chk.why), 'why is empty on a pass');

Yr = Y;  Yr(3:4, :) = Yr(3:4, :).*(1 + 1e-6*t/t(end));   % energy ramp 2e-6 relative
chkR = mech.conservation_check(t, Yr, sys);
ok = check(ok, ~chkR.ok && ~chkR.okEnergy && chkR.okConstraint, 'injected energy ramp FAILS the energy bar only');
ok = check(ok, contains(chkR.why, 'energy'), sprintf('why names energy: "%s"', chkR.why));

Yc = Y;  Yc(3:4, :) = Yc(3:4, :) + 1e-6*Y(1:2, :);         % radial velocity 1e-6
chkC = mech.conservation_check(t, Yc, sys);
ok = check(ok, ~chkC.ok && ~chkC.okConstraint && chkC.okEnergy, 'injected radial velocity FAILS the constraint bar only');
ok = check(ok, contains(chkC.why, 'constraint'), sprintf('why names the constraint: "%s"', chkC.why));

chkB = mech.conservation_check(t, Yr, sys, struct('energy', 1e-3));
ok = check(ok, chkB.ok, 'a looser energy bar passes the ramp');

ok = check(ok, throws(@() mech.conservation_check(t, Y.', sys), 'mech:conservation_check:shape'), ...
           'transposed Y refused with :shape');
ok = check(ok, throws(@() mech.conservation_check(t(1:10), Y, sys), 'mech:conservation_check:shape'), ...
           't and Y of different lengths refused with :shape');

if ok, fprintf('TEST_CONSERVATION_CHECK: ALL PASS\n'); else, fprintf('TEST_CONSERVATION_CHECK: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
