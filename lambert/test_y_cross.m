r1n=1; r2n=1.2; A=sqrt(1.2); mu=1; dt=2.0;
for z=-5:0.5:-2
  fprintf('z=%f, t=%f\n', z, lambert_tof(z, r1n, r2n, A, mu));
end
