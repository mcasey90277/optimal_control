function ok = test_cartpole_pmp()
%% Purpose:
%
%   The whole indirect solve, and the claims the spec makes for it:
%     1. it converges, and the terminal state is hit to 1e-9;
%     2. the Hamiltonian H = u^2 + lam'(F + G u) is CONSTANT along the arc
%        -- but this flow is HAMILTONIAN, so H is conserved along EVERY
%        trajectory of it, whatever lam(0) is; constancy therefore cannot
%        tell a converged extremal from a wrong seed basin. What it does
%        catch: a wrong costate rate or a mis-scaled control, either of
%        which breaks the conservation law that the pointwise stationarity
%        identity (below) cannot see;
%     3. an ARITHMETIC-CONSISTENCY check on the reported (U, Lam) pair, NOT
%        an optimality test: out.U was built as u = -lam'G/2 from the same
%        Lam, X and cartpole_field this check reuses, and H is exactly
%        quadratic in u, so a centred finite difference of H in u equals
%        2u + lam'G identically -- it cannot fail by construction. Kept as
%        a cross-check that the stored (X, Lam, U) triple is internally
%        self-consistent with the field, not as evidence of optimality;
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
addpath(here, fullfile(root, 'oclib'));

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

%% Arithmetic-consistency check, NOT an optimality test (see item 3 above):
%% a CENTRED FINITE DIFFERENCE of H in u, evaluated at the stored u0 =
%% out.U(k) = -lam'G/2, is 2*u0 + lam'G identically -- zero by construction
%% for this exactly-quadratic H, regardless of whether (X, Lam) solve the
%% true PMP. What a nonzero result WOULD catch is a mismatch between out.U
%% and out.Lam/out.X (e.g. a stale or mis-indexed pairing), not a wrong
%% extremal.
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
         sprintf(['arithmetic consistency of (U, Lam), not optimality: centred FD of ', ...
                  '2u + lam''G = %.2e (identically zero by construction)'], worstFD));

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
