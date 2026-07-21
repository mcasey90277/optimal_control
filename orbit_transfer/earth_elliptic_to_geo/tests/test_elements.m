% TEST_ELEMENTS  Forward conversion invariants + roundtrip.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
p  = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;  e0 = 0.75;  a0 = P0/(1-e0^2);
% (a) paper initial state: apogee, i=7 deg
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0.0612, 0, pi, p.mu);
assert(abs(norm(r0) - a0*(1+e0)) < 1e-10, 'not at apogee radius');
assert(abs(dot(r0,v0)) < 1e-10,           'radial rate nonzero at apsis');
vis = sqrt(p.mu*(2/norm(r0) - 1/a0));
assert(abs(norm(v0) - vis) < 1e-10,       'vis-viva violated');
hv = cross(r0, v0);
assert(abs(acosd(hv(3)/norm(hv)) - 7.0052) < 1e-3, 'inclination wrong');
% (b) GEO check: equatorial circular prograde at any L
[rg, vg] = elements_to_cart(1, 0, 0, 0, 0, 0.7, p.mu);
assert(abs(norm(rg)-1) < 1e-12 && abs(norm(vg)-1) < 1e-12);
assert(abs(rg(3)) < 1e-12 && abs(vg(3)) < 1e-12 && abs(dot(rg,vg)) < 1e-12);
hg = cross(rg, vg);  assert(hg(3) > 0, 'retrograde GEO');
% (c) roundtrip on a grid (incl. wrap-aware L compare)
rng(7);
for kk = 1:50
    el = [0.2+rand, 0.7*(rand-0.5), 0.7*(rand-0.5), 0.2*(rand-0.5), ...
          0.2*(rand-0.5), 2*pi*rand-pi];
    [rr, vv] = elements_to_cart(el(1), el(2), el(3), el(4), el(5), el(6), p.mu);
    eb = cart_to_elements(rr, vv, p.mu);
    assert(max(abs([eb.P eb.ex eb.ey eb.hx eb.hy] - el(1:5))) < 1e-10, 'roundtrip elems');
    dL = mod(eb.L - el(6) + pi, 2*pi) - pi;
    assert(abs(dL) < 1e-10, 'roundtrip L');
end
fprintf('test_elements: ALL PASS\n');
