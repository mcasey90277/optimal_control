function out = solve_drones_direct(x0, xT, geom, opts)
%% Purpose:
%
%   The multi-agent kinematic drone problem (Zhang, Liu, Darbon &
%   Karniadakis, arXiv:2610.30388, eq. 19) as a trapezoidal direct
%   transcription solved by IPOPT through CasADi:
%
%       min  integral_0^T sum_i 1/2 |v_i|^2 dt
%       s.t. xdot_i = v_i,  x_i(0) = x0_i,  x_i(T) = xT_i,
%            |x_i - x_j| >= 2 r  (mode 'hard': at every node and every
%                                 Hermite midpoint, squared form)
%            or the paper's C1 log barrier added to the running cost
%                                (mode 'barrier': the smooth twin Route 2
%                                 shoots; no inequalities)
%            |x_i - c| >= R_ball  when geom.ball is set (3D)
%
%   2D layouts (geom.d == 2) pin v_z = 0 by equality so no agent can dodge
%   in z. The initial guess is the straight line plus a tangential bump of
%   amplitude opts.bump with a common rotation sense: every circle layout
%   is degenerate at the straight line (all agents meet at the centre at
%   T/2), and the mirror image of any solution is a solution.
%
%  ASSUMPTIONS / NOTES:
%
% - The defect constraints are added FIRST so their multipliers are
%   opti.lam_g(1:3MN), read by row range (never opti.dual: canonicalised
%   orientation, see memory casadi-optidual-sign-bug).
% - The barrier expression is a second copy of barrier_shape's formulas
%   (CasADi needs if_else); out.Unodes lets a test pin it to barrier_terms.
% - A sampled path constraint holds between samples only to O(h^2); the
%   clearance gate in verify_drones_direct reports the fine-grid minimum.
%
%% Inputs:
%
%  x0, xT                   [3 x M]                 From drone_layouts
%
%  geom                     struct                  From drone_layouts
%
%  opts                     struct (optional)       .N [100] .T [10]
%                                                   .mode ['hard'|'barrier']
%                                                   .sense [+1] .bump
%                                                   [0.05] .pairs
%                                                   [pair_list(M)]
%                                                   .midpoints [true]
%                                                   .maxIter [3000]
%                                                   .print [0] .X0 .V0
%                                                   (warm start)
%                                                   .bumpScale [ones(1,M)]
%
%   opts.bumpScale          [1 x M]                 Per-agent multipliers on
%                                                   the bump amplitude of the
%                                                   default guess (agent a's
%                                                   bump is bump*bumpScale(a));
%                                                   every entry must be > 0,
%                                                   else the error
%                                                   solve_drones_direct:bump.
%                                                   Alternating values seed a
%                                                   symmetry-broken basin
%                                                   (drone_table1's multistart)
%
%   opts.X0, opts.V0        [3M x N+1] each         Warm start. BOTH or
%                                                   neither: one alone, or
%                                                   either of the wrong shape,
%                                                   is refused by the error
%                                                   solve_drones_direct:warm.
%                                                   With a warm start IPOPT
%                                                   runs warm_start_init_point
%                                                   = 'yes' and mu_init = 1e-4
%                                                   (a cold restart from
%                                                   mu = 0.1 leaves a saddle
%                                                   guess by luck, not by
%                                                   descent); cold solves keep
%                                                   IPOPT's defaults
%
%% Outputs:
%
%  out                      struct                  .X .V [3M x N+1] .tN
%                                                   .J .JLB .success
%                                                   .ipoptStatus
%                                                   .iterations .maxDefect
%                                                   .lamDefect [3M x N]
%                                                   .minHnodes .Unodes
%                                                   .mode .N .T .geom
%                                                   .pairs .sense .bump
%                                                   .bumpScale .warm .wall
%
%   out.lamDefect           [3M x N]                opti.lam_g on the defect
%                                                   rows (row range, never
%                                                   opti.dual). Column k
%                                                   samples +V at the MIDPOINT
%                                                   of interval k: the
%                                                   adjacent-column average
%                                                   (mu_{k-1} + mu_k)/2
%                                                   equals V at node k
%                                                   (V-stationarity of the
%                                                   Lagrangian): 2.8e-17 on
%                                                   the barrier fixtures; in
%                                                   hard mode ~1e-7 except at
%                                                   nodes next to an ACTIVE
%                                                   midpoint clearance (the
%                                                   Hermite midpoint depends
%                                                   on V, so its multiplier
%                                                   enters: 1.8e-4 on the
%                                                   M = 8 hard fixture).
%                                                   For Route 2 the PMP
%                                                   costate is lambda =
%                                                   -lamDefect stationed at
%                                                   the interval midpoints;
%                                                   the matching oclib map is
%                                                   oc.duals_to_costates with
%                                                   scheme 'trapezoid-nodal'
%                                                   (adjacent-interval average
%                                                   at the N+1 nodes), NOT
%                                                   'trapezoid' (left-node
%                                                   stationing: a half-step
%                                                   phase error)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the paper's M = 4 circle, hard clearance:
     [a, b, g] = drone_layouts('paper', 4);
     o = solve_drones_direct(a, b, g);
     fprintf('M = 4: J = %.6f (J_LB %.4f, OCFormer 0.2003), %s\n', o.J, o.JLB, o.ipoptStatus);
     return
end
if nargin < 4, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
casadi_path();
import casadi.*

d = @(f, v) fieldd(opts, f, v);
N      = d('N', 100);
T      = d('T', 10);
mode   = lower(d('mode', 'hard'));
sense  = d('sense', 1);
bump   = d('bump', 0.05);
M      = geom.M;
pairs  = d('pairs', pair_list(M));
midpts = d('midpoints', true);
maxIt  = d('maxIter', 3000);
prnt   = d('print', 0);
bScale = d('bumpScale', ones(1, M));
if ~(bump > 0)
    error('solve_drones_direct:bump', 'the tangential bump must be positive: the straight line is degenerate');
end
if numel(bScale) ~= M || ~all(bScale(:) > 0)
    error('solve_drones_direct:bump', 'opts.bumpScale must hold M = %d entries, every one > 0', M);
end
bScale = reshape(bScale, 1, M);
if ~any(strcmp(mode, {'hard', 'barrier'}))
    error('solve_drones_direct:mode', 'mode must be ''hard'' or ''barrier''');
end

nX = 3*M;  nN = N + 1;  h = T/N;
tN = linspace(0, T, nN);
P  = size(pairs, 2);
t0 = tic;

%% Variables and the defects (FIRST, for the multiplier rows)
opti = Opti();
X = opti.variable(nX, nN);
V = opti.variable(nX, nN);
opti.subject_to(X(:, 2:end) - X(:, 1:end-1) - (h/2)*(V(:, 1:end-1) + V(:, 2:end)) == 0);
opti.subject_to(X(:, 1) == x0(:));
opti.subject_to(X(:, end) == xT(:));
if geom.d == 2
    opti.subject_to(V(3:3:end, :) == 0);
end

%% Cost: trapezoid of 1/2 |v|^2
L = 0.5*sum1(V.^2);                           % [1 x nN]
J = h*sum2(L) - (h/2)*(L(1) + L(end));

%% Pair geometry at the nodes and (optionally) the Hermite midpoints
ia = reshape(3*(pairs(1,:) - 1) + [1; 2; 3], [], 1);
ib = reshape(3*(pairs(2,:) - 1) + [1; 2; 3], [], 1);
Xmid = (X(:, 1:end-1) + X(:, 2:end))/2 + (h/8)*(V(:, 1:end-1) - V(:, 2:end));
if midpts, Xs = [X, Xmid]; else, Xs = X; end
Unode = [];
if P > 0
    Dp = Xs(ia, :) - Xs(ib, :);                                        % [3P x ns]
    sq = Dp(1:3:end, :).^2 + Dp(2:3:end, :).^2 + Dp(3:3:end, :).^2;   % [P x ns]
    switch mode
        case 'hard'
            opti.subject_to(sq(:) >= (2*geom.r)^2);
        case 'barrier'
            hh = sqrt(sq(:, 1:nN)) - 2*geom.r;                        % nodes only
            Ub = barrier_casadi(hh, geom.eps, geom.l);
            Unode = sum1(Ub);                                          % [1 x nN]
    end
end
if ~isempty(geom.ball)
    Db = Xs - repmat(geom.ball.c, M, size(Xs, 2));
    sqb = Db(1:3:end, :).^2 + Db(2:3:end, :).^2 + Db(3:3:end, :).^2;  % [M x ns]
    switch mode
        case 'hard'
            opti.subject_to(sqb(:) >= geom.ball.R^2);
        case 'barrier'
            hb = sqrt(sqb(:, 1:nN)) - geom.ball.R;
            Ubb = sum1(barrier_casadi(hb, geom.eps, geom.l));
            if isempty(Unode), Unode = Ubb; else, Unode = Unode + Ubb; end
    end
end
if strcmp(mode, 'barrier') && ~isempty(Unode)
    J = J + h*sum2(Unode) - (h/2)*(Unode(1) + Unode(end));
end
opti.minimize(J);

%% Initial guess: a warm start (both X0 and V0), or straight lines plus the tangential bump
X0 = d('X0', []);  V0 = d('V0', []);
warm = ~isempty(X0) || ~isempty(V0);
if warm && (isempty(X0) || isempty(V0) || ~isequal(size(X0), [nX nN]) || ~isequal(size(V0), [nX nN]))
    error('solve_drones_direct:warm', ['a warm start needs BOTH opts.X0 and opts.V0, each [3M x N+1] = [%d x %d] ' ...
          '(got X0 %s, V0 %s)'], nX, nN, mat2str(size(X0)), mat2str(size(V0)));
end
if ~warm
    [X0, V0] = bumped_lines(x0, xT, tN, bump, sense, bScale);
    if strcmp(mode, 'barrier') && P > 0
        [~, rhoG] = clearance_terms(X0, pairs, geom);
        if min(rhoG(:)) < 1e-9
            error('solve_drones_direct:guess', ['barrier mode: the default initial guess has coincident agents ' ...
                  '(ring layouts); pass a warm start opts.X0/opts.V0, e.g. the hard-mode solution']);
        end
    end
end
opti.set_initial(X, X0);
opti.set_initial(V, V0);

%% Solve
ip = struct('print_level', prnt, 'max_iter', maxIt, 'tol', 1e-10, ...
            'acceptable_tol', 1e-8, 'linear_solver', 'mumps', ...
            'bound_relax_factor', 0);   % constraints hold to tol, not to a relaxed bound; the nodal clearance gate is -1e-9
if warm
    % a warm start must stay warm: start the barrier parameter small instead
    % of IPOPT's 0.1 (which re-centres the iterate and, from a guess on a
    % symmetric saddle, picks a basin by luck)
    ip.warm_start_init_point = 'yes';
    ip.mu_init = 1e-4;
end
opti.solver('ipopt', struct('print_time', false), ip);
try
    sol = opti.solve();  solved = true;
catch
    sol = opti.debug;    solved = false;
end

%% Report, with the defects re-evaluated from the returned numbers
out = struct();
out.X  = full(sol.value(X));
out.V  = full(sol.value(V));
out.tN = tN;
out.J  = full(sol.value(J));
out.JLB = lower_bound(x0, xT, T);
out.ipoptStatus = opti.stats().return_status;
out.iterations  = opti.stats().iter_count;
out.success = solved && any(strcmp(out.ipoptStatus, {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));
def = out.X(:, 2:end) - out.X(:, 1:end-1) - (h/2)*(out.V(:, 1:end-1) + out.V(:, 2:end));
out.maxDefect = max(abs(def(:)));
lamAll = full(sol.value(opti.lam_g));
out.lamDefect = reshape(lamAll(1:nX*N), nX, N);
if P > 0
    hn = clearance_terms(out.X, pairs, geom);
    out.minHnodes = min(hn(:));
else
    out.minHnodes = Inf;
end
if isempty(Unode), out.Unodes = []; else, out.Unodes = full(sol.value(Unode)); end
out.mode = mode;  out.N = N;  out.T = T;  out.geom = geom;  out.pairs = pairs;
out.sense = sense;  out.bump = bump;  out.bumpScale = bScale;  out.midpoints = midpts;
out.warm = warm;
out.wall = toc(t0);
end

% ------------------------------------------------------------------------
function U = barrier_casadi(hh, epsB, l)
%% Purpose:
%
%   The paper's C1 log barrier as a CasADi expression (a copy of
%   barrier_shape's two formulas; barrier_shape itself branches on
%   real(h), which has no symbolic meaning). Pinned to barrier_terms by
%   test_direct_small through out.Unodes.
%
import casadi.*
U = if_else(hh > l, -epsB*log(hh), -epsB*log(l) + (epsB/2)*(((hh - 2*l)/l).^2 - 1));
end

function [X0, V0] = bumped_lines(x0, xT, tN, bump, sense, bScale)
%% Purpose:
%
%   Straight lines at constant speed plus a tangential bump
%   A_a sin(pi t/T) n_perp with a common rotation sense about e3, so the
%   agents start swirling the same way instead of meeting at the centre;
%   A_a = bump*bScale(a) (bScale all ones = the uniform guess).
%
M  = size(x0, 2);  T = tN(end);  s = tN/T;
X0 = zeros(3*M, numel(tN));  V0 = X0;
for a = 1:M
    dv = xT(:, a) - x0(:, a);
    L = sqrt(sum(dv.^2));
    if L > 0, dir = dv/L; else, dir = [1; 0; 0]; end
    perp = [-dir(2); dir(1); 0];                       % e3 x dir
    if sqrt(sum(perp.^2)) < 1e-6
        perp = [0; -dir(3); dir(2)];                   % e1 x dir, for a z-aligned line
    end
    perp = sense*perp/sqrt(sum(perp.^2));
    ia = 3*(a-1) + (1:3);
    A  = bump*bScale(a);
    X0(ia, :) = x0(:, a) + dv*s + A*perp*sin(pi*s);
    V0(ia, :) = dv/T + A*perp*(pi/T)*cos(pi*s);
end
end

function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
