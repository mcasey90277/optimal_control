r1 = [1;0;0]; r2 = [0;1.2;0]; mu = 1;
[V1, V2, Ns, out] = lambert_uv_multirev(r1, r2, 40, mu, +1, 20);
fprintf('%.2f\n', out.zs(1))
for i=1:out.Nmax
  idx = 2 + 2*(i-1);
  fprintf('%.2f/%.2f\n', out.zs(idx), out.zs(idx+1))
end
