function ok = test_drone_physics()
%% Purpose:
%
%   THE INDEPENDENT PHYSICAL ORACLE for the drone models. Every other test
%   in this folder compares the field against itself (a complex step of it,
%   a generated Jacobian of it); this one derives its expectations from
%   mechanics and geometry and never from the implemented rows. The house
%   rule behind it: on 2026-09-17 a cart-pole sign error survived nine
%   self-comparing reviews, and a power-balance identity found it.
%
%   Checks:
%     1. kinematic: F = 0 and B = I exactly (the model is its own oracle).
%     2. pointmass: free fall vdot = -g e3; hover a = g e3 holds v; power
%        balance d/dt(1/2 |v|^2 + g pz) = a . v pointwise.
%     3. rigid6dof hover: f_k = m g / 4, Theta = 0, omega = 0, v = 0 is an
%        equilibrium.
%     4. rigid6dof free fall: f = 0 gives vdot = -g e3 at every attitude.
%     5. rigid6dof POWER BALANCE, pointwise at random states and controls:
%            d/dt (1/2 m |v|^2 + m g pz + 1/2 omega' I omega)
%              = T (R e3) . v + tau . omega
%        The gyroscopic term omega x I omega does no work, so it cancels;
%        a sign error in any row breaks the identity.
%     6. rigid6dof KINEMATIC CONSISTENCY, pointwise: Rdot = R [omega]x with
%        Rdot = sum_k dR/dTheta_k Thetadot_k -- pins W independently of the
%        rates row.
%     7. rigid6dof torque-free flight (equal rotor thrusts) conserves the
%        INERTIAL angular momentum R I omega; a hover-trim flight stays put.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(31);  restore = onCleanup(@() rng(rngWas));
e3 = [0; 0; 1];

%% 1. kinematic
pk = drone_params('kinematic');
[F, B] = drone_field(randn(3,1), pk);
ok = chk(ok, isequal(F, zeros(3,1)) && isequal(B, eye(3)), 'kinematic: F = 0, B = I exactly');

%% 2. pointmass
pq = drone_params('pointmass');
worstFall = 0;  worstHover = 0;  worstPow = 0;
for trial = 1:20
    x = randn(6,1);  a = randn(3,1);
    [F, B] = drone_field(x, pq);
    xd = F + B*a;
    worstFall  = max(worstFall, max(abs(F(4:6) + pq.g*e3)));
    xdH = F + B*(pq.g*e3);
    worstHover = max(worstHover, max(abs(xdH(4:6))));
    v = x(4:6);
    lhs = v.'*xd(4:6) + pq.g*xd(3);                 % d/dt (1/2|v|^2 + g pz) along the field
    rhs = a.'*v;
    worstPow = max(worstPow, abs(lhs - rhs)/max(1, abs(rhs)));
end
ok = chk(ok, worstFall < 1e-15, sprintf('pointmass free fall: vdot = -g e3 (%.1e)', worstFall));
ok = chk(ok, worstHover < 1e-15, sprintf('pointmass hover a = g e3 holds v (%.1e)', worstHover));
ok = chk(ok, worstPow < 1e-13, sprintf('pointmass power balance d/dt(1/2|v|^2 + g pz) = a.v (%.1e)', worstPow));

%% 3. rigid6dof hover equilibrium
pr = drone_params('rigid6dof');
xH = [randn(3,1); zeros(9,1)];
fH = (pr.m*pr.g/4)*ones(4,1);
[F, B] = drone_field(xH, pr);
ok = chk(ok, max(abs(F + B*fH)) < 1e-15, 'rigid6dof hover (f_k = m g/4, level, at rest) is an equilibrium');

%% 4-6. rigid6dof pointwise identities at random states and controls
worstFall = 0;  worstPow = 0;  worstKin = 0;
Imat = diag(pr.I);
for trial = 1:50
    x = randn(12,1);  x(7:8) = x(7:8)*0.5;           % roll/pitch within +/- ~1.2 rad
    f = pr.fMax*rand(4,1);
    [F, B] = drone_field(x, pr);
    xd = F + B*f;
    worstFall = max(worstFall, max(abs(F(4:6) + pr.g*e3)));
    v = x(4:6);  Th = x(7:9);  om = x(10:12);
    [R, ~, dR] = euler_zyx(Th);
    T   = sum(f);
    tau = [pr.l*(f(2) - f(4)); pr.l*(f(3) - f(1)); pr.kappa*(f(1) - f(2) + f(3) - f(4))];
    % 5. power balance: E = 1/2 m|v|^2 + m g pz + 1/2 om' I om, dE/dt along the field
    lhs = pr.m*(v.'*xd(4:6)) + pr.m*pr.g*xd(3) + om.'*(Imat*xd(10:12));
    rhs = T*((R*e3).'*v) + tau.'*om;
    worstPow = max(worstPow, abs(lhs - rhs)/max(1e-3, abs(rhs)));
    % 6. kinematic consistency: Rdot from Thetadot vs R [omega]x
    Thd = xd(7:9);
    Rdot = dR{1}*Thd(1) + dR{2}*Thd(2) + dR{3}*Thd(3);
    omx = [0 -om(3) om(2); om(3) 0 -om(1); -om(2) om(1) 0];
    worstKin = max(worstKin, max(max(abs(Rdot - R*omx))));
end
ok = chk(ok, worstFall < 1e-15, sprintf('rigid6dof free fall: vdot = -g e3 at every attitude (%.1e)', worstFall));
ok = chk(ok, worstPow < 1e-10, sprintf('rigid6dof POWER BALANCE dE/dt = T (R e3).v + tau.omega (%.1e rel)', worstPow));
ok = chk(ok, worstKin < 1e-12, sprintf('rigid6dof KINEMATICS Rdot = R [omega]x (%.1e)', worstKin));

%% 7. torque-free flight conserves R I omega; hover trim stays put
x0 = [zeros(6,1); 0.3; -0.2; 0.5; 1.0; -0.7; 0.4];
fEq = 0.05*ones(4,1);                                % equal thrusts: tau = 0, T = 0.2 N
oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, X] = ode113(@(t, x) fieldAt(x, fEq, pr), [0 2], x0, oo);
L0 = euler_zyx(x0(7:9))*(Imat*x0(10:12));
Lend = euler_zyx(X(end,7:9).')*(Imat*X(end,10:12).');
ok = chk(ok, max(abs(Lend - L0))/max(abs(L0)) < 1e-8, ...
         sprintf('torque-free flight conserves inertial angular momentum R I omega (%.1e rel over 2 s)', ...
                 max(abs(Lend - L0))/max(abs(L0))));
[~, XH] = ode113(@(t, x) fieldAt(x, fH, pr), [0 2], xH, oo);
ok = chk(ok, max(abs(XH(end,:).' - xH)) < 1e-12, 'hover trim holds for 2 s');

if ok, fprintf('TEST_DRONE_PHYSICS: ALL PASS\n'); else, fprintf('TEST_DRONE_PHYSICS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function xd = fieldAt(x, u, p)
%% Purpose:
%
%   xdot = F(x) + B(x) u for the integrator.
%
[F, B] = drone_field(x, p);
xd = F + B*u;
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
