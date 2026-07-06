r1 = [1;0;0]; r2 = [0;1.2;0]; mu = 1;

[v1, v2, info] = lambert_uv(r1, r2, 2.0, mu, +1);
v1', v2', info.z

pk_v1 = [ 1.7292152422437140e-01;  9.9659460660332588e-01; 0];
pk_v2 = [-8.3049550550277162e-01; -6.8224231238170905e-03; 0];
disp([norm(v1-pk_v1), norm(v2-pk_v2)]);

f2b = @(t,x) [x(4:6); -mu*x(1:3)/norm(x(1:3))^3];
sol = ode89(f2b, [0 2], [r1; v1], odeset('RelTol',3e-14,'AbsTol',1e-14));
xf = deval(sol, 2);  
disp([norm(xf(1:3)-r2), norm(xf(4:6)-v2)]);

disp('C2-C3')
[v1h,~,ih] = lambert_uv(r1, r2, 0.5, mu, +1);   % faster than parabolic
disp([ih.z, norm(v1h)^2/2 - mu/1])                       % energy > 0?
[v1r,~,ir] = lambert_uv(r1, r2, 2.0, mu, -1);
ir.dtheta                                        % expect 3*pi/2

disp('C4')
[v1v, v2v] = lambert_uv([15945.34;0;0], ...
    [12214.83899;10249.46731;0], 4560, 398600.4418, +1);
v1v', v2v'

