function ok = test_huber_saltation()
%% Purpose:
%
%   The STM returned by cr3bp_minfuel_prop must equal the true sensitivity
%   of the state map, INCLUDING across a throttle DISCONTINUITY. The 'huber'
%   family jumps at Q = 1 (s: p*Q -> 1); integrating Phi_dot = A*Phi with
%   the branchwise AD Jacobian misses the switching-time sensitivity, and
%   the correct STM carries the saltation update
%       Phi+ = [ I + (F+ - F-) n' / (n' F-) ] Phi-,   n = grad Q,
%   at each crossing (review 2026-09-05, P0.2, sol + Astra). The 'eps'
%   family is continuous at its clip knees, so its STM needs no update --
%   it serves as the control case for the finite-difference harness.
%
%   Fixture: the field's own demo state, dt = 0.6, huber p = 0.3 -- the
%   arc crosses Q = 1 once, at t ~ 0.153 (measured). Reference: central
%   finite differences of the 14-state map with the same integrator
%   tolerances (RelTol 1e-10), which resolve the switch by step control.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All cases passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
if isempty(which('cr3bp_minfuel_prop'))
    addpath(fileparts(fileparts(mfilename('fullpath'))));
end
mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
y0  = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
       15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
dt  = 0.6;
h   = 1e-6;

fams = {struct('family','eps','p',0.3),    'eps    (continuous, control case)', 1e-5; ...
        struct('family','huber','p',0.3),  'huber  (Q=1 jump at t~0.153)',       1e-5; ...
        struct('family','huberc','p',0.3), 'huberc (hybrid: continuous ramp, no saltation needed)', 1e-5};
for k = 1:size(fams, 1)
    sm = fams{k,1};
    [~, PHI] = cr3bp_minfuel_prop(dt, y0, true, T_, c_, mu_, sm);
    PHIfd = zeros(14);
    for j = 1:14
        e = zeros(14,1);  e(j) = h;
        yp = cr3bp_minfuel_prop(dt, y0 + e, false, T_, c_, mu_, sm);
        ym = cr3bp_minfuel_prop(dt, y0 - e, false, T_, c_, mu_, sm);
        PHIfd(:, j) = (yp - ym) / (2*h);
    end
    err = max(abs(PHI(:) - PHIfd(:))) / max(abs(PHIfd(:)));
    ok = chk(ok, err < fams{k,3}, ...
             sprintf('%s: max|PHI - PHI_fd| / max|PHI_fd| = %.2e (tol %.0e)', ...
                     fams{k,2}, err, fams{k,3}));
end

%% Astra review #2 (2026-09-06): switch-surface kinematics and edge cases
% (a) Qdot closed form: Qdot = -Tmax lam_v'lam_r / (m |lam_v|), independent of
%     the throttle -- so n'F+ == n'F- at every switch.
%     Two independent checks: the chain rule n'F with the eps field, and a
%     SHORT forward difference (Q is strongly curved at this state:
%     dt = 1e-3 is 10x off, dt = 1e-7 agrees to 1e-3 -- measured).
Qd = cr3bp_minfuel_qdot(y0, T_, c_);
Fe = cr3bp_minfuel_pmp(y0, T_, c_, mu_, struct('family','eps','p',0.3));
rho0 = norm(y0(11:13));  n0 = zeros(14,1);
n0(7) = -T_*rho0/y0(7)^2;  n0(11:13) = T_*y0(11:13)/(y0(7)*rho0);  n0(14) = T_/c_;
Qof = @(y) T_*(norm(y(11:13))/y(7) + y(14)/c_);
yh = cr3bp_minfuel_prop(1e-7, y0, false, T_, c_, mu_, struct('family','eps','p',0.3));
Qfd = (Qof(yh) - Qof(y0)) / 1e-7;
ok = chk(ok, abs(Qd - n0'*Fe) < 1e-10 && abs(Qd - Qfd) < 2e-3*abs(Qd), ...
         sprintf('closed-form Qdot %.6f == n''F %.6f; short FD %.6f', Qd, n0'*Fe, Qfd));

% (b) starting EXACTLY on the switch surface with Qdot > 0 must take the HI
%     branch immediately (s = 1 just after t = 0), not default to 'lo'.
y1 = y0;  y1(14) = y1(14) + (1 - T_*(norm(y0(11:13))/y0(7) + y0(14)/c_))*c_/T_;   % set Q(y1) = 1 exactly
assert(abs(T_*(norm(y1(11:13))/y1(7) + y1(14)/c_) - 1) < 1e-12, 'fixture: Q != 1');
assert(cr3bp_minfuel_qdot(y1, T_, c_) > 0, 'fixture: need Qdot > 0 at the start');
[~, ~, Ts, Ys] = cr3bp_minfuel_prop(0.05, y1, false, T_, c_, mu_, struct('family','huber','p',0.3));
[~, ~, axs] = cr3bp_minfuel_pmp(Ys(2,:)', T_, c_, mu_, struct('family','huber','p',0.3));
ok = chk(ok, abs(axs.s - 1) < 1e-12 && Ts(end) > 0.05 - 1e-12, ...
         sprintf('start at Q = 1 with Qdot > 0: hi branch immediately (s = %.4f), reached dt', axs.s));

% (c) backward propagation is not supported by the event-split path:
try
    cr3bp_minfuel_prop(-0.1, y0, false, T_, c_, mu_, struct('family','huber','p',0.3));
    ok = chk(ok, false, 'huber dt < 0 rejected');
catch E
    ok = chk(ok, contains(E.identifier, 'cr3bp_minfuel_prop'), sprintf('huber dt < 0 rejected (%s)', E.identifier));
end

if ok, fprintf('TEST_HUBER_SALTATION: ALL PASS\n');
else,  fprintf('TEST_HUBER_SALTATION: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
