r1 = [1;0;0]; r2 = [0;1.2;0]; mu = 1;
[V1, V2, Ns, out] = lambert_uv_multirev(r1, r2, 40, mu, +1, 20);
out.Nmax
size(V1,2)
disp('tmins:'); disp(out.tmins)

disp('D2')
[~,~,~,o7] = lambert_uv_multirev(r1, r2, 7.0, mu, +1, 20);
[o7.Nmax, o7.tmins(1)]
