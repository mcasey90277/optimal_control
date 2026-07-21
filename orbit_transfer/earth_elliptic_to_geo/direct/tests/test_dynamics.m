% TEST_DYNAMICS  Ballistic invariants + thrust mass-rate exactness.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
p = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;  a0 = P0/(1-0.75^2);
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0.0612, 0, pi, p.mu);
x0 = [r0; v0; 1; 0];
% (a) one ballistic period: energy/|h| conserved, state returns
T0 = 2*pi*a0^1.5;
oo = odeset('RelTol',1e-12,'AbsTol',1e-13);
[~, xx] = ode113(@(t,x) lt2b_rhs_time(x, [1;0;0;0], p), [0 T0], x0, oo); % s=0: alpha moot
xe = xx(end,:).';
E  = @(x) 0.5*dot(x(4:6),x(4:6)) - p.mu/norm(x(1:3));
assert(abs(E(xe) - E(x0)) < 1e-9, 'energy drift');
assert(norm(cross(xe(1:3),xe(4:6)) - cross(r0,v0)) < 1e-9, 'h drift');
assert(norm(xe(1:6) - x0(1:6)) < 1e-6, 'period return failed');
assert(abs(xe(7) - 1) < 1e-14 && abs(xe(8) - T0) < 1e-9, 'm/t states wrong');
% (b) full thrust for 1 TU: exact linear mass, energy increases (tangential)
odef = @(t,x) lt2b_rhs_time(x, [x(4:6)/norm(x(4:6)); 1], p);
[~, xt] = ode113(odef, [0 1], x0, oo);
assert(abs(xt(end,7) - (1 - p.Tmax/p.c)) < 1e-10, 'mass rate wrong');
assert(E(xt(end,:).') > E(x0), 'tangential thrust must raise energy');
fprintf('test_dynamics: ALL PASS\n');
