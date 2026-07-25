% TEST_FOC_CHECK_TOY  foc_check on a self-built bang-bang toy (no campaign deps).
% min int u dt, xdot=v, vdot=u, u in [0,1], x(0)=v(0)=0, x(T)=1, v(T) FREE, T=2.
% Optimal: u=1 then u=0 (single switch) -- lam_v linear hits 0 slope, S=1+lam_v.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_verify_common;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));  import casadi.*
N = 60;  T = 2;  sg = linspace(0, 1, N+1).';  ds = diff(sg).';
opti = casadi.Opti();
X = opti.variable(2, N+1);  U = opti.variable(1, N+1);
creg = struct('label',{},'rows',{});
f = [X(2,:); U(1,:)];                              % [v; u]
r0 = size(opti.g,1)+1;
for k = 1:N
    opti.subject_to(X(:,k+1) - X(:,k) - (ds(k)*T/2)*(f(:,k)+f(:,k+1)) == 0);
end
creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1));
r0 = size(opti.g,1)+1;  opti.subject_to(U(1,:).' >= 0);
creg(end+1) = struct('label','thrLo','rows',r0:size(opti.g,1));
r0 = size(opti.g,1)+1;  opti.subject_to(U(1,:).' <= 1);
creg(end+1) = struct('label','thrHi','rows',r0:size(opti.g,1));
opti.subject_to(X(:,1) == [0;0]);  opti.subject_to(X(1,end) == 1);
opti.minimize(sum((ds*T/2).*(U(1,1:N)+U(1,2:N+1))));
opti.set_initial(X, [linspace(0,1,N+1); 0.5*ones(1,N+1)]);
opti.set_initial(U, 0.5*ones(1,N+1));
opti.solver('ipopt', struct('print_time',false), struct('print_level',0));
sol = opti.solve();
out = struct('X', full(sol.value(X)), 'U', full(sol.value(U)), ...
             'model', struct('opti', opti, 'creg', creg));
rep = foc_check(out, sg, foc_manifest('toy'), struct());
assert(rep.kktStatInf < 1e-6, 'KKT stationarity: %.2e', rep.kktStatInf);
assert(rep.signPct >= 99, 'sign law: %.1f%%', rep.signPct);
assert(rep.nSwitches == 1, 'expected 1 switch, got %d', rep.nSwitches);
assert(rep.singularArcNodes == 0, 'no singular arc expected');
assert(rep.sdotMinRel > 1e-3, 'Sdot at the switch must be regular');
assert(rep.pass, 'advisory pass expected on the toy');
fprintf('test_foc_check_toy: ALL PASS\n');
