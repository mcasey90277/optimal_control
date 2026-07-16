% TEST_P2_HOMOTOPY  Paper toy P2: energy->fuel homotopy vs analytic bang-off-bang.
cp = getenv('CASADI_PATH'); if isempty(cp), cp = fullfile(getenv('HOME'),'casadi-3.7.0'); end
addpath(cp);
N = 200;  dt = 2/N;
opti = casadi.Opti();
X  = opti.variable(2, N+1);                 % [x; v]
Up = opti.variable(1, N+1);  Um = opti.variable(1, N+1);
u  = Up - Um;
opti.subject_to(Up >= 0);  opti.subject_to(Up <= 1);
opti.subject_to(Um >= 0);  opti.subject_to(Um <= 1);
for k = 1:N     % trapezoid defects for [xdot; vdot] = [v; u]
    fk  = [X(2,k);   u(k)];
    fk1 = [X(2,k+1); u(k+1)];
    opti.subject_to(X(:,k+1) - X(:,k) - (dt/2)*(fk+fk1) == 0);
end
opti.subject_to(X(:,1)   == [0; 0]);
opti.subject_to(X(:,end) == [0.5; 0]);
opti.set_initial(X, [linspace(0,0.5,N+1); zeros(1,N+1)]);
runc = @(w,epsv) w - epsv*(w.*(1-w));       % per-component homotopy integrand
for epsv = [1 0.6 0.3 0.12 0.04 0.01 0]
    g = runc(Up,epsv) + runc(Um,epsv);
    opti.minimize( sum((dt/2)*(g(1:N)+g(2:N+1))) );
    opti.solver('ipopt', struct('print_time',0), struct('print_level',0,'max_iter',800));
    sol = opti.solve();
    opti.set_initial(X,  sol.value(X));
    opti.set_initial(Up, sol.value(Up));
    opti.set_initial(Um, sol.value(Um));
end
uv = sol.value(Up) - sol.value(Um);
tg = linspace(0,2,N+1);
cost = trapz(tg, abs(uv));
assert(abs(cost - (2-sqrt(2))) < 3e-3, 'P2 cost mismatch: %.5f', cost);
assert(uv(1) > 0.95 && uv(end) < -0.95, 'not bang at ends');
assert(abs(uv(round(N/2))) < 0.05, 'not coasting at midpoint');
kSw = find(uv < 0.5, 1);                    % first departure from the +1 arc
assert(abs(tg(kSw) - (1-1/sqrt(2))) < 0.03, 'switch time off');
fprintf('test_p2_homotopy: ALL PASS (cost %.6f vs %.6f)\n', cost, 2-sqrt(2));
