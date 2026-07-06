r1=[1;0;0]; r2=[0;1.2;0]; mu=1;
[v1,v2,info] = lambert_uv(r1, r2, 0.01, mu, +1);
info.z
