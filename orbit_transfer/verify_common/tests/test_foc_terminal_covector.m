% TEST_FOC_TERMINAL_COVECTOR  Mapped terminal covector (Hager 2000) vs the raw
% one-sided interval dual and its I5 linear extrapolation, on a toy double
% integrator with an analytically known zero terminal costate.
%
% min integral 0.5*u^2 dt, xdot=v, vdot=u, x(0)=v(0)=0, x(T)=-1, v(T) FREE,
% T=2. Extends the test_foc_check_toy pattern (same hand-built opti+creg, no
% campaign deps) with a genuinely free terminal state that carries NO
% terminal cost -- v never appears in the objective, only in xdot=v -- which
% is exactly the man.massFreeAtTf semantics foc_check gates on.
%
% PMP: H = 0.5*u^2 + lam_x*v + lam_v*u => u = -lam_v (stationarity).
% lam_x_dot = 0 => lam_x = c1 (const). lam_v_dot = -lam_x = -c1
%   => lam_v(t) = c1*(T-t).
% Transversality lam_v(T) = 0 is therefore AUTOMATIC in this parametrization
% (not an extra condition to solve for) -- x(T)=-1 alone fixes c1 via
% x(T) = -c1*T^3/3 = -1  =>  c1 = 3/T^3 = 3/8 (T=2). With this sign, u(t) =
% c1*(t-T) <= 0 and vdot = u <= 0 on [0,T], so v(t) is MONOTONICALLY
% NONINCREASING -- required by foc_check's massRow-nonincreasing manifest
% guard, which is why v (not x) stands in for "mass" here.
%
% WHY THE MESH IS DELIBERATELY COARSE (N=20, vs N=60 in test_foc_check_toy).
% The raw one-sided interval dual is not itself a nodal costate sample: it is
% Lambda_N, the multiplier of the LAST DEFECT INTERVAL, and it behaves like a
% sample of lam_v near sigma_end - h_N/2 rather than AT sigma_end (see
% foc_dual_to_costate.m header). Since lam_v is linear with slope -c1 here,
% that interval-vs-endpoint offset alone predicts a raw reading of order
% c1*h_N/2 -- an O(h) artifact that shrinks (and gets harder to demonstrate
% cleanly against solver/AD floor noise) as the mesh refines. A coarse mesh
% makes the offset large and unambiguous relative to the mapped covector,
% which this test asserts sits at KKT-stationarity floor regardless of h.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_verify_common;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));  import casadi.*

N = 20;  T = 2;  sg = linspace(0, 1, N+1).';  ds = diff(sg).';
opti = casadi.Opti();
X = opti.variable(2, N+1);  U = opti.variable(1, N+1);
creg = struct('label',{},'rows',{});
f = [X(2,:); U(1,:)];                              % [v; u]
r0 = size(opti.g,1)+1;
for k = 1:N
    opti.subject_to(X(:,k+1) - X(:,k) - (ds(k)*T/2)*(f(:,k)+f(:,k+1)) == 0);
end
creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1));
opti.subject_to(X(:,1) == [0;0]);
opti.subject_to(X(1,end) == -1);                   % v(T) left FREE, no terminal cost
opti.minimize(sum((ds*T/2).*(0.5*U(1,1:N).^2 + 0.5*U(1,2:N+1).^2)));
opti.set_initial(X, [linspace(0,-1,N+1); -0.375*ones(1,N+1)]);
opti.set_initial(U, -0.1*ones(1,N+1));
opti.solver('ipopt', struct('print_time',false), struct('print_level',0));
sol = opti.solve();
out = struct('X', full(sol.value(X)), 'U', full(sol.value(U)), ...
             'model', struct('opti', opti, 'creg', creg));

% Ad hoc manifest: v (state row 2) plays the "mass" role -- free at t_f, no
% terminal cost, nonincreasing (see PMP derivation above). No thrust/direction
% rows apply to this smooth (non-bang-bang) toy, so those checks are skipped
% ([] fields), exactly per foc_manifest's own "skip, don't fake" convention.
man = struct('name','toy_freemass', 'nx',2, 'nu',1, 'dirRows',[], 'thrRow',[], ...
    'massRow',2, 'timeRow',[], 'autonomous',true, 'horizonKind','none', ...
    'massFreeAtTf',true, 'throttleCostKind','affine');

rep = foc_check(out, sg, man, struct());

fprintf('test_foc_terminal_covector: N=%d (coarse)\n', N);
fprintf('  mapped        |lam_v(tf)| rel : %.3e\n', rep.lamMassEndMapped);
fprintf('  extrapolated  |lam_v(tf)| rel : %.3e   [I5]\n', rep.lamMassEndRel);
fprintf('  raw one-sided |lam_v(tf)| rel : %.3e   [legacy]\n', rep.lamMassEndRelOneSided);

% (a) the mapped covector's mass component is the exact discrete object and
% must sit at KKT-stationarity floor, independent of mesh coarseness.
assert(rep.lamMassEndMapped < 1e-8, ...
    'test_foc_terminal_covector: mapped terminal covector should be ~0 by KKT stationarity, got %.3e', ...
    rep.lamMassEndMapped);

% (b) the raw one-sided interval dual must be demonstrably LARGER on this
% coarse mesh -- proving the two genuinely differ and that "mapped" is not
% just a relabeling of the same number.
assert(rep.lamMassEndMapped * 10 < rep.lamMassEndRelOneSided, ...
    ['test_foc_terminal_covector: mapped (%.3e) should be more than 10x smaller ' ...
     'than the raw one-sided dual (%.3e) on this deliberately coarse (N=%d) mesh -- ' ...
     'if this fails the O(h) endpoint offset was not reproduced, try a coarser N'], ...
    rep.lamMassEndMapped, rep.lamMassEndRelOneSided, N);

fprintf('test_foc_terminal_covector: ALL PASS\n');
