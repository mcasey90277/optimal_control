% TEST_TERMINAL  Fixed-GEO state properties + manifold residuals.
p = kepler_lt_params(10, 1500, 2000);
tf1 = geo_terminal('fixed', p, 2.3);
rv = tf1.rvf;  r = rv(1:3);  v = rv(4:6);
assert(abs(norm(r)-1) < 1e-12 && abs(norm(v)-1) < 1e-12);
assert(abs(r(3)) < 1e-12 && abs(v(3)) < 1e-12 && abs(dot(r,v)) < 1e-12);
tm = geo_terminal('manifold', p, []);
assert(max(abs(tm.resid(rv))) < 1e-12, 'GEO state must satisfy manifold');
rvBad = rv;  rvBad(1) = rvBad(1)*1.01;       % radius + radial-rate violated
res = tm.resid(rvBad);
assert(abs(res(3)) > 1e-3, 'radius constraint insensitive');
fprintf('test_terminal: ALL PASS\n');
