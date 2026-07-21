function xdot = lt2b_rhs_time(x, u, par)
% LT2B_RHS_TIME  Time-domain EOM: inertial 2-body gravity + low thrust (8-state).
%
% x = [r(3); v(3); m; t], u = [alpha(3); s] with ||alpha||=1, s in [0,1].
% Thrust accel = (Tmax/m)*s*alpha; mdot = -(Tmax/c)*s; tdot = 1 (time carried
% as a state so the Sundman solver can pin t(tau_f)). Written without
% norm/abs/max so it evaluates on BOTH numeric doubles and CasADi MX.
%
% INPUTS:
%   x   - State vector [r(3); v(3); m; t] [8x1]
%   u   - Control [alpha(3); s] with ||alpha||=1, s in [0,1] [4x1]
%   par - Structure from kepler_lt_params (contains mu, Tmax, c)
%
% OUTPUTS:
%   xdot - Time derivative of state [8x1] = d/dt [r; v; m; t]
%
% REFERENCES:
%   [1] process/DESIGN.md sec 2 (problem statement).

r = x(1:3);  v = x(4:6);  m = x(7);
rn2 = r(1)^2 + r(2)^2 + r(3)^2 + 1e-12;      % softened, AD/CS-safe
acc = -par.mu * r * rn2^(-1.5) + (par.Tmax/m) * u(4) * u(1:3);
xdot = [v; acc; -(par.Tmax/par.c)*u(4); 1];
end
