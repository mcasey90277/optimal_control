function ok = test_cartpole_pmp()
%% Purpose:
%
%   The whole indirect solve, and the claims the spec makes for it:
%     1. it converges, and the terminal state is hit to 1e-9;
%     2. the Hamiltonian H = u^2 + lam'(F + G u) is CONSTANT along the arc
%        -- an autonomous, fixed-t_f PMP extremal must conserve it; this
%        catches a wrong costate rate, a wrong seed basin, or a mis-scaled
%        control that the pointwise stationarity identity (below) cannot;
%     3. dH/du = 0 at the reported control, checked by a CENTRED FINITE
%        DIFFERENCE that perturbs u directly -- independent of the
%        u* = -lam'G/2 formula that built U, unlike an algebraic
%        recomputation of the same expression;
%     4. flying the recovered control through the true dynamics
%        (oc.fly_control) arrives where the solve says -- the G1b idea;
%     5. the cost agrees with the direct fixture to 1% (the gap is the
%        direct method's discretization error; the indirect solve is the
%        more accurate of the two);
%     6. the control profiles agree in shape to 2% RMS;
%     7. regression: lam0 reproduces the stored reference to 1e-8;
%     8. the same answer comes back with a different segment count -- the
%        solution is the problem's, not the discretization's.
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
root = fileparts(fileparts(here));                 % optimal_control
addpath(here, fullfile(root, 'oclib'), fullfile(root, 'orbit_transfer', 'costate_common'));

out = run_cartpole_pmp(struct('K', 8, 'plot', false));
R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));
p = R.p;

ok = chk(ok, out.info.converged, 'ms_bvp reports converged');
ok = chk(ok, out.missTerminal < 1e-9, sprintf('terminal miss %.2e < 1e-9', out.missTerminal));

%% Hamiltonian constancy: H = u^2 + lam'(F + G u) is conserved along an
%% autonomous, fixed-t_f extremal. Built from out.X/out.Lam/out.U directly
%% through cartpole_field -- an independent consequence of the solve, not
%% an identity of how U was computed.
M = numel(out.t);
H = zeros(1, M);
for k = 1:M
    [Fk, Gk] = cartpole_field(out.X(:,k), p);
    H(k) = out.U(k)^2 + out.Lam(:,k).'*(Fk + Gk*out.U(k));
end
relRangeH = (max(H) - min(H)) / max(abs(H));
ok = chk(ok, relRangeH < 1e-9, ...
         sprintf('Hamiltonian constant along the arc: relative range %.2e', relRangeH));

%% dH/du = 0 by CENTRED FINITE DIFFERENCE, perturbing u directly at a
%% handful of sample times -- does not reuse the u* = -lam'G/2 formula.
du = 1e-4;
idxFD = round(linspace(2, M-1, 7));
worstFD = 0;
for k = idxFD
    [Fk, Gk] = cartpole_field(out.X(:,k), p);
    u0 = out.U(k);
    Hp = (u0 + du)^2 + out.Lam(:,k).'*(Fk + Gk*(u0 + du));
    Hm = (u0 - du)^2 + out.Lam(:,k).'*(Fk + Gk*(u0 - du));
    worstFD = max(worstFD, abs((Hp - Hm)/(2*du)));
end
ok = chk(ok, worstFD < 1e-6, ...
         sprintf('dH/du = 0 by finite difference: worst |dH/du| = %.2e (curvature d2H/du2 = 2)', worstFD));

ok = chk(ok, out.missFlown < 1e-6, sprintf('flown control arrives: %.2e', out.missFlown));

relJ = abs(out.J - R.J)/R.J;
ok = chk(ok, relJ < 0.01, sprintf('cost agrees with the direct solution: J = %.6f vs %.6f (%.3f%%)', ...
                                  out.J, R.J, 100*relJ));

Uref = interp1(R.tN, R.U, out.t, 'linear');
rmsU = sqrt(mean((out.U - Uref).^2))/sqrt(mean(Uref.^2));
ok = chk(ok, rmsU < 0.02, sprintf('control profiles agree: %.3f%% RMS', 100*rmsU));

ref = load(fullfile(here, 'data', 'cartpole_pmp_ref.mat'));
ok = chk(ok, max(abs(out.lam0 - ref.lam0)) < 1e-8, ...
         sprintf('regression: lam0 reproduces (%.1e)', max(abs(out.lam0 - ref.lam0))));

out16 = run_cartpole_pmp(struct('K', 16, 'plot', false));
ok = chk(ok, max(abs(out16.lam0 - out.lam0)) < 1e-7, ...
         sprintf('K = 16 finds the same extremal (%.1e)', max(abs(out16.lam0 - out.lam0))));

if ok, fprintf('TEST_CARTPOLE_PMP: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
