function zdot = ocp_dynamics(~, z, mu)
% OCP_DYNAMICS  Coupled state-costate dynamics for the min-energy transfer.
%
%   State z = [x; lam] in R^8 with x = [r; v] (position, velocity) and
%   lam = [lam_r; lam_v] the costates. The optimal control has been
%   substituted in closed form: u* = -lam_v (the primer vector).
%
%     r'     = v
%     v'     = g(r) + u*,     u* = -lam_v
%     lam_r' = -G(r)' * lam_v
%     lam_v' = -lam_r
%
% INPUTS:
%   ~   - Time (unused; the system is autonomous)
%   z   - Stacked state-costate vector [8x1]: [x1 x2 x3 x4 l1 l2 l3 l4]'
%   mu  - Gravitational parameter [scalar]
%
% OUTPUTS:
%   zdot - Time derivative [8x1]
%
% REFERENCES:
%   [1] myLatex/notes/cov_pmp_orbit_transfer.tex, Sec. 6 (worked problem).

r_vec = z(1:2);
v_vec = z(3:4);
lam_r = z(5:6);
lam_v = z(7:8);

u     = -lam_v;                          % optimal control = primer vector
g     = two_body_accel(r_vec, mu);
G     = gravity_gradient(r_vec, mu);

zdot = [ v_vec;
         g + u;
         -G' * lam_v;
         -lam_r ];
end
