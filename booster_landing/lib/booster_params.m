function P = booster_params()
% BOOSTER_PARAMS  Falcon-9-class 3-DOF landing-burn parameter set.
%
% Single source of truth for physical constants, boundary conditions,
% solver grid and Monte-Carlo settings. Public F9 estimates; every number
% adjustable here and only here.
%
% INPUTS:  none
% OUTPUTS: P - parameter struct (fields documented inline below)
%
% REFERENCES:
%   [1] Blackmore, Acikmese, Scharf, "Minimum-Landing-Error Powered-Descent
%       Guidance for Mars Landing Using Convex Optimization," JGCD 2010.
%   [2] docs/superpowers/specs/2026-08-08-booster-landing-design.md

%% Vehicle (public Falcon 9 block-5 estimates):
P.mdry   = 25600;              % dry mass [kg]
P.m0     = 30000;              % mass at landing-burn start [kg]
P.Tmax   = 845e3;              % one Merlin 1D, sea level [N]
P.Tmin   = 0.40 * P.Tmax;      % ~40 percent min throttle [N]
P.etaT   = 0.93;               % GUIDANCE thrust de-rate [-] (ADJUDICATED
                                % 2026-08-08, task-7b). The guidance solves
                                % against etaT*Tmax = 785.85 kN; the tracker
                                % and the truth sim keep the FULL [Tmin,Tmax]
                                % annulus, so the 7% is reserved headroom for
                                % feedback rather than lost capability.
                                % Rationale: the min-fuel solution rode Tmax
                                % for 54% of its nodes (certify G5 bound
                                % fraction 0.9917) with only ~2% of net
                                % deceleration to spare, so there was NOTHING
                                % left to add during braking -- a -5% thrust
                                % dispersion was unrecoverable by ANY tracker
                                % around that nominal (measured 67.4 m/s
                                % touchdown, identical under two different
                                % controllers, i.e. arithmetic not tracking).
                                % De-rating the guidance buys back real
                                % authority: the braking arc's net decel goes
                                % from ~21.5 to ~19.3 m/s^2 nominal, leaving
                                % the tracker ~2.2 m/s^2 to command against
                                % dispersions. Costs propellant (see the
                                % task-7 report's 7b section for the measured
                                % fuel delta) -- that is the trade the
                                % adjudication accepted.
P.Isp    = 282;                % sea-level Isp [s]
P.g0     = 9.80665;            % standard gravity [m/s^2]
P.gvec   = [0; 0; -P.g0];      % flat-Earth gravity, z up [m/s^2]

%% Boundary conditions (pad at origin, z up):
P.r0     = [500; 100; 2000];   % post-entry-burn position [m]
P.v0     = [-30; 0; -180];     % descending ~180 m/s [m/s]
P.vf     = [0; 0; -1.5];       % terminal velocity target [m/s] (ADJUDICATED
                                % 2026-08-08, task-7 fix report round 3: was
                                % v(tf)=0. With P.Tmin exceeding vehicle
                                % weight at every mass on this trajectory,
                                % the v=0-at-z=0 endpoint is singular -- the
                                % closed loop cannot complete the last
                                % stretch of descent under pure vertical
                                % thrust and instead arrests just above the
                                % pad (measured: 0.53 m). Touching down
                                % descending at 1.5 m/s (legs absorb it;
                                % P.vtd_max=2.0 mission gate unchanged)
                                % substantially improves closed-loop
                                % touchdown behavior (a genuine landing near
                                % this target IS achievable, which the old
                                % v(tf)=0 endpoint never allowed at any
                                % gain). CORRECTION (round-3 report,
                                % retracted there, corrected here per round-
                                % 4 review): this does NOT "remove the
                                % arrest state entirely" as originally
                                % claimed -- P.Tmin still exceeds weight at
                                % every mass regardless of the target
                                % velocity, so closed-loop tracking error
                                % that falls behind the guidance's own
                                % descent-rate schedule can still arrest
                                % near the ground under some gain/grid
                                % settings (see sim_closed_loop.m's
                                % TERMINAL-PHASE note and P.zTermBand below
                                % for the structural fix that round 4
                                % applied on top of this BC change).

%% Path constraints:
P.gs_deg        = 30;          % glideslope min elevation angle [deg]
P.theta_max_deg = Inf;         % thrust-pointing cone half-angle, Inf = off

%% Discretization / solver:
P.N      = 60;                 % Hermite-Simpson segments (collocation)
P.Nconv  = 120;                % trapezoid nodes for the convex solver
P.tf_lo  = 10;                 % free-final-time bracket [s]
P.tf_hi  = 22;                 % ADAPTATION (task-7 fix report round 3,
                                % 2026-08-08): was 50. A dedicated sweep
                                % after the P.vf terminal-BC change (above)
                                % found solve_pdg_convex's fixed-tf
                                % subproblem hits a genuine, reproducible
                                % Infeasible_Problem_Detected wall around
                                % tf~27 s (loitering that long burns past
                                % P.mdry -- Tmin exceeds weight everywhere,
                                % so there is no free coast), with a
                                % numerically marginal/nondeterministic
                                % band just below that (tf~25-26 s:
                                % lossless_gap measured both ~6e-5 tight
                                % and ~0.45 untight for the IDENTICAL call
                                % in separate runs -- IPOPT/BLAS-threading
                                % sensitivity, not a real feasibility
                                % change). The old tf_hi=50 golden-section
                                % bracket's fixed initial probes (at the
                                % golden-ratio points of [tf_lo,tf_hi])
                                % landed exactly in that marginal-to-
                                % infeasible zone and intermittently failed
                                % outright. Both solvers agree the true
                                % mass-optimal tf is ~15.6-16 s (unimodal,
                                % confirmed by a tf sweep); 22 keeps
                                % generous headroom around that peak while
                                % moving the golden search's own probe
                                % points (~14.6, ~17.4) into the
                                % consistently well-behaved region. Also
                                % narrows solve_pdg_colloc's free-tf upper
                                % bound and its Tc=(tf_lo+tf_hi)/2 ND time
                                % scale (30->16 s) closer to the true
                                % answer -- a tighter, more physically
                                % grounded characteristic scale, not a
                                % looser one. NOTE for Phase 2 (P.drag.on):
                                % this bracket was validated in vacuum only
                                % -- drag will change the fuel-vs-tf curve
                                % (and hence where the tf~27 s infeasibility
                                % wall sits), so tf_hi may need widening
                                % when drag is switched on; re-run the tf
                                % sweep in this comment's block rather than
                                % assume 22 still has headroom.

P.zTermBand = 150;             % terminal-phase altitude gate [m] (task-7
                                % fix report round 4, 2026-08-08): below
                                % this altitude, sim_closed_loop.m's
                                % control_law stops feeding z-position error
                                % into the thrust command and instead tracks
                                % the guidance's own v(z) profile -- the
                                % structural fix for the closed-loop arrest
                                % mechanism (see that file's TERMINAL-PHASE
                                % note). 150 m is the reviewer-suggested
                                % starting point (this trajectory's descent
                                % rate near the end is ~15-20 m/s, so 150 m
                                % is roughly the last ~8-10 s of NOMINAL
                                % flight time before the final high-gain
                                % braking arc -- not tightly tuned; a sweep
                                % of this value was not run this round).

%% Atmosphere (Phase 2, OFF by default -- vacuum keeps convexification exact):
P.drag.on   = false;
P.drag.rho0 = 1.225;           % sea-level density [kg/m^3]
P.drag.H    = 8500;            % scale height [m]
P.drag.Cd   = 1.0;             % landing-leg config drag coefficient [-]
P.drag.A    = 10.75;           % reference area, 3.7 m diameter [m^2]

%% Success criteria / Monte Carlo:
P.pad_radius = 15;             % landing accuracy requirement [m]
P.vtd_max    = 2.0;            % max touchdown speed [m/s]
P.seed       = 42;             % rng seed, all random draws
end
