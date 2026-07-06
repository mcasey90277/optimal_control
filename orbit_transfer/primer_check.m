function [err_primer, err_uv] = primer_check(sol, tf, mu)
% PRIMER_CHECK  Verify Lawden's primer equation along the indirect solution.
%
%   Independent verification of the costate block structure: extract the
%   primer initial conditions from the converged BVP solution,
%
%     p(0)    = -lam_v(0),      pdot(0) = +lam_r(0)
%
%   (from the definition p = -lam_v and the costate equation
%   lam_v-dot = -lam_r: pdot = -lam_v-dot = +lam_r),
%
%   propagate the second-order primer ODE  p'' = G(r(t)) p  as its own
%   4-state system alongside nothing but the STATE trajectory r(t), and
%   compare the result with -lam_v(t) from the full state-costate flow.
%   Agreement to integration tolerance proves p'' = Gp holds along the
%   trajectory -- the primer really does propagate under the gravity
%   gradient of the nominal path.
%
% INPUTS:
%   sol - ode45 solution structure from SOLVE_INDIRECT (8-state flow)
%   tf  - Final time [scalar]
%   mu  - Gravitational parameter [scalar]
%
% OUTPUTS:
%   err_primer - max over the grid of ||p_ode(t) - (-lam_v(t))|| [scalar]
%   err_uv     - max over the grid of ||u(t) - p(t)|| with u = -lam_v and
%                p from the primer ODE; identical quantity by construction
%                of u*, so this equals err_primer -- reported separately to
%                make the "thrust rides the primer" statement explicit.
%
% REFERENCES:
%   [1] Lawden, "Optimal Trajectories for Space Navigation", 1963.
%   [2] myLatex/notes/cov_pmp_orbit_transfer.tex, Secs. 6-7.

z0    = deval(sol, 0);
p0    = -z0(7:8);                    % p(0)    = -lam_v(0)
pdot0 = z0(5:6);                     % pdot(0) = lam_r(0)

opts  = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
psol  = ode45(@(t, w) primer_ode(t, w, sol, mu), [0 tf], [p0; pdot0], opts);

tq = linspace(0, tf, 1001);
zq = deval(sol,  tq);
wq = deval(psol, tq);

p_from_costate = -zq(7:8, :);        % -lam_v(t) from the coupled flow
p_from_ode     = wq(1:2, :);         % p(t) from the standalone primer ODE
u_traj         = -zq(7:8, :);        % u(t) = -lam_v(t)

err_primer = max(vecnorm(p_from_ode - p_from_costate));
err_uv     = max(vecnorm(u_traj     - p_from_ode));
end

% ----------------------------------------------------------------------
function wdot = primer_ode(t, w, sol, mu)
% First-order form of p'' = G(r(t)) p; r(t) read from the nominal solution.
z     = deval(sol, t);
G     = gravity_gradient(z(1:2), mu);
wdot  = [ w(3:4);
          G * w(1:2) ];
end
