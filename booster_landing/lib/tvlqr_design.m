function ctrl = tvlqr_design(sol, P, opts)
% TVLQR_DESIGN  Time-varying LQR about a PDG guidance trajectory.
%
% Linearizes pdg_dynamics along (x*(t), T*(t)) and integrates the
% differential Riccati equation backward:
%   -Pdot = A'P + PA - PBR^{-1}B'P + Q,   P(tf) = Qf
% Gains K(t) = R^{-1} B(t)' P(t), stored on a dense grid; controller is
%   T_cmd(t) = T*(t) - K(t) (x - x*(t))  (saturation applied by the sim).
%
% The state weight is PHASE-SCHEDULED (see ADAPTATION 5): Q = Q(t), high
% position weight while the guidance rides Tmin (thrust margin available
% for lateral correction), position weight collapsed across the switch
% into the Tmax braking arc (where every degree of tilt is stolen straight
% out of vertical braking).
%
% INPUTS:
%   sol  - collocation solution (Task 3 interface)
%   P    - booster_params
%   opts - (optional) overrides:
%          .omegaA .zetaA  phase-A lateral pole placement [def 0.70, 1.00]
%                    -- QA is DERIVED from these plus R and the phase-A
%                       mass (see below); zetaA >= 1/sqrt(2) required
%          .QA (7x7) margin-phase state weight   [def: derived, see below]
%          .QB (7x7) braking-phase state weight  [def: see below]
%          .Q  (7x7) legacy CONSTANT weight -- sets QA=QB=Q, i.e. disables
%                    the phase schedule (kept so callers/tests written
%                    against the pre-schedule interface still work)
%          .R  (3x3) control weight [def anisotropic, see ADAPTATION 5]
%          .Qf (7x7) terminal weight
%          .tSwitch [s] override for the auto-detected Tmin->Tmax switch
%          .tBlend  [s] tanh half-width of the QA->QB blend
% OUTPUTS:
%   ctrl - .tgrid(1xM) .K(3x7xM) .Pt(7x7xM) .xnom(@t->7x1)
%          .Tnom(@t_scalar->3x1) -- SCALAR-t ONLY (see note below); every
%          caller in this repo (ricrhs below, sim_closed_loop's
%          control_law) already calls it with a scalar t.
%          .tSwitch .tBlend .Qfun(@t->7x7) -- the schedule actually flown
%          (diagnostics; not required by any downstream task interface).
%
% ADAPTATION (task-7 fix report, 2026-08-08): ctrl.Tnom was originally a
% global pchip spline threaded through all nodes+midpoints. certify_pdg.m's
% G2 gate already disproved that exact representation for reconstructing a
% Hermite-Simpson control (task-5 fix report): a global spline is not the
% per-segment quadratic Simpson's rule is actually built against, and
% measured a non-shrinking ~0.6-0.8 kg mass "floor" that the per-segment
% reconstruction collapsed to ~0.0001-0.0005 kg. ctrl.Tnom now calls the
% shared lib/hs_quad_ctrl.m (promoted out of certify_pdg.m, byte-identical
% logic, regression-checked by test_certify_nominal.m) instead, so the
% TVLQR feedforward flies the same control G2 already certified accurate
% to ~1.6e-4 m, not a spline known to smear the bang-bang switch and
% front-load braking near touchdown. hs_quad_ctrl is scalar-t only (see
% its own header) -- ctrl.Tnom inherits that restriction; ctrl.xnom stays
% a pchip state interpolant (states are smooth between HS nodes, no
% representation mismatch there).
%
% ADAPTATION 2 (task-7 fix report round 2, 2026-08-08): default R raised
% from 1e-10 to 1.5e-9 at the time, under a v(tf)=0 terminal BC. Root
% cause: under the deep terminal-arc saturation this trajectory forces
% (P.Tmin exceeds vehicle weight at every mass -- see sim_closed_loop.m),
% a saturated control law only has DIRECTION left to work with, not
% magnitude; at R=1e-10 the direction is chosen almost entirely to null a
% large dispersed horizontal error, leaving little of the saturated
% thrust vector's direction "budget" for arresting vertical speed.
%
% ADAPTATION 3 (task-7 fix report round 3, 2026-08-08): default R lowered
% again, to 7e-10, after P.vf changed the guidance terminal velocity from
% 0 to [0;0;-1.5] m/s (booster_params.m; the v=0 endpoint was singular
% under P.Tmin>weight). A fresh mini-sweep under the new BC (R from 1e-12
% to 2e-9, plus a Qf(4:6) probe at R=7e-10 -- task-7 fix report round 3
% for the full table) found the dispersed-50m-safe window narrowed to
% R in [6.5e-10, 7.0e-10] (R<=6e-10 fails vtd>2; R>=7.5e-10 reverts to an
% arrest, td.landed=false, even though its vtd number looks fine).
% R=7e-10 sits at the good end of that window with markedly more margin
% than round 2's R=1.5e-9 gave under the old BC (dispersed miss=9.64 m
% vs pad_radius=15, vtd=0.88 vs vtd_max=2.0 -- round 2's setting was
% razor-thin at miss=14.77/vtd=1.377).
%
% The zero-dispersion case STILL does not reach a genuine touchdown at
% ANY R tested in the dispersed-safe window (arrest ~0.65 m up, vtd~0.64
% m/s, far from the P.vf=-1.5 m/s target) -- removing the v(tf)=0
% singularity substantially improved nominal tracking (a LOW R, ~1e-11 to
% 7e-11, now DOES achieve a genuine landing near the -1.5 m/s target) but
% did not eliminate the arrest possibility under closed-loop tracking
% error, and that low-R band is disjoint from the dispersed-safe band by
% over an order of magnitude. See the task-7 fix report round 3 for the
% full sweep. -- ROUND-4 UPDATE: this round-3 "disjoint band" finding
% turned out to be partly an N=30-grid ARTIFACT, not the physics itself
% (see below); the ROUND-3 default (R=7e-10) is superseded here.
%
% ADAPTATION 4 (task-7 fix report round 4, 2026-08-08): default R lowered
% to 1e-9. Root cause of round 3's "nominal never lands" finding, found
% by a reviewer and confirmed here: at the N=30 grid the CLOSED-LOOP TEST
% used (not the N=P.N=60 production grid the shipped defaults are meant
% for), hs_quad_ctrl's PRE-round-4 reconstruction dipped up to 18% below
% Tmin over 11.7% of the flight between nodes/midpoints -- invisible to
% every gate that existed before round 4's new G2ff check (certify_pdg.m).
% The closed-loop sim's saturating clamp turned that into a real,
% grid-dependent +0.6-1 m open-loop altitude bias; two rounds of R tuning
% were unknowingly buying feedback gain to cancel that feedforward defect,
% not tuning against real closed-loop physics. hs_quad_ctrl is now
% annulus-feasible BY CONSTRUCTION at any grid (see that function), and
% the closed-loop test now uses N=P.N=60 (tests/test_closed_loop_nominal.m)
% -- at that combination, the NOMINAL case lands cleanly (alt~0, vtd in
% [1.4,2.1]) at EVERY R from 7e-10 to ~1e-7, confirmed by a fresh sweep.
%
% The DISPERSED-case R antagonism, however, is REAL and SURVIVES both the
% hs_quad_ctrl fix and sim_closed_loop.m's then-current altitude-scheduled
% terminal-phase remedy (zeroed z-position feedback + v(z) tracking below
% P.zTermBand). HISTORICAL: that remedy no longer exists -- task-7b
% replaced it with full altitude-indexed guidance, which makes the
% zeroed-altitude-channel and the "past-tf freeze" both structurally
% impossible rather than patched, and retires P.zTermBand entirely. See
% sim_closed_loop.m's ALTITUDE-INDEXED GUIDANCE note for what actually
% runs. At N=60, post-both-fixes: dispersed miss<15
% needs R<~1.5e-9 (R=1.5e-9: miss=13.99, right at the edge; R=3e-9:
% miss=20.4, over) while dispersed vtd<2 needs R>~5e-8 (R=3e-8: vtd=1.68;
% R=1e-8: vtd=2.80, over) -- a ~20-30x gap with NO overlap in an
% extensive sweep (task-7 fix report round 4 for the full table). Per
% that round's explicit instruction ("if the antagonism survives the
% structural remedy, STOP with the sweep table"), this is reported, not
% forced: R=1e-9 is chosen as the shipped default because it gets the
% NOMINAL case fully right (landed=true, vtd=1.51 m/s, comfortable margin
% under vtd_max=2.0) and the DISPERSED miss comfortably inside pad_radius
% (10.85 m vs 15 m) -- the geometric "did it reach the pad" requirement --
% while being honest that dispersed vtd (6.82 m/s) badly fails the
% vtd_max=2.0 safety gate at this setting. No R in the tested range
% satisfies both dispersed gates at once; see the report for candidate
% structural remedies beyond a weight/schedule retune (anisotropic R,
% integral action, or revisiting whether a single 3-DOF thrust vector can
% ever simultaneously correct a 50 m lateral offset AND arrest vertical
% speed to 2 m/s under P.Tmin>weight -- a genuine actuator-authority
% question, not a tuning one).
%
% ADAPTATION 5 (task-7 fix report round 5, 2026-08-08): the round-4
% "genuine actuator-authority limit" verdict on the dispersed 58 m case is
% WRONG, and the sweeps that produced it were searching the wrong axis.
% Measured root cause (diag, N=60, dsp.dr0=[50;-30;0], old defaults):
%
%   * the closed-loop LATERAL position-error pole is BOUNDED ABOVE by
%     sqrt(qPos/qVel), for EVERY R. Exactly (see the ARE inversion in
%     part (a) below),
%         pole = a/b = sqrt( qPos / (2 m sqrt(qPos r) + qVel) ),
%     which is strictly DECREASING in r and rises to its supremum
%     sqrt(qPos/qVel) only as r -> 0. At the old
%     Q=diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0]) that supremum is
%     sqrt(1e-4/1e-2) = 0.1 rad/s -- a 10 s time constant on a 15.6 s
%     flight, and NO choice of R can beat it. (At the shipped r=1e-9 the
%     2 m sqrt(qPos r) term is 0.0187, i.e. ~1.9x LARGER than qVel=0.01,
%     so the pole is actually 0.059 rad/s, further below the bound;
%     measured 0.084 at t=0, where the Riccati is still off its
%     steady-state value.) THAT is why every R sweep in rounds 2-4 moved
%     `miss` monotonically and never fixed it: lowering R buys gain
%     MAGNITUDE on both channels but cannot push the pole past a ceiling
%     fixed by the qPos:qVel ratio alone. The search axis was orthogonal
%     to the defect.
%   * consequently the tracker asks for K(1,1)*58 m ~ 18 kN of lateral
%     thrust at t=0 while the guidance sits AT Tmin=338 kN with 507 kN of
%     unused magnitude margin -- it declines free authority. The 58 m
%     offset was still 32 m at the Tmin->Tmax switch (t~6.8 s) and 10.9 m
%     at touchdown; the residual lateral demand then had to be paid out of
%     the saturated Tmax braking arc, where a 16 deg tilt (measured)
%     costs cos(16 deg) of vertical braking -- which is exactly the
%     dispersed vtd=6.82 m/s failure. Not an authority limit: an
%     authority-SCHEDULING failure.
%
% Fix, two parts:
%  (a) PHASE-SCHEDULED Q. Phase A (t < tSwitch, the Tmin margin arc) gets
%      lateral weights designed by pole placement rather than guessed. For
%      the lateral double integrator (A=[0 1;0 0], B=[0;1/m]) the LQR
%      closed loop is s^2 + (Kv/m)s + (Kr/m) with
%        Kr/m = sqrt(qPos/r)/m,   Kv/m = sqrt((2m sqrt(qPos r)+qVel)/r)/m,
%      so asking for natural frequency w and damping zeta inverts to
%        qPos = m^2 r w^4,        qVel = m^2 r w^2 (4 zeta^2 - 2)
%      (note qVel>=0 forces zeta>=1/sqrt(2) -- LQR's own Butterworth floor).
%      Shipped: w=0.70 rad/s, zeta=1.00, with m taken as the mean mass
%      over the margin arc and r=R(1,1) read from the ACTUAL weight, both
%      at call time (opts.omegaA/.zetaA expose the design point). At the
%      default R that evaluates to qPos=0.2099, qVel=0.8569. The design
%      formula is validated, not assumed:
%      the Riccati solution measures K(1,1)/m = 0.4808 at t=0 against the
%      predicted w^2 = 0.49 (and K(1,4)/m = 1.377 against 2*zeta*w = 1.40).
%      Feasibility check that sets the ceiling on w: the initial lateral
%      command is w^2*|dr| = 0.49*58 = 28 m/s^2 ~ 840 kN, which together
%      with the ~338 kN vertical feedforward sits just outside the
%      Tmax=845 kN rim -- so w=0.70 spends the margin arc's spare authority
%      and just barely more, with sim_closed_loop's vertical-priority
%      allocate_thrust absorbing the small overshoot without touching the
%      vertical channel. The correction is essentially complete by t~5 s,
%      well before tSwitch=7.16 s. w much beyond 0.75 buys nothing (the
%      command is then annulus-limited, not gain-limited) and w below ~0.6
%      leaves residual offset for the braking arc to pay for.
%      Phase B (t > tSwitch, the Tmax braking arc) keeps the PRE-round-5
%      weights verbatim, so the braking arc still tracks the velocity
%      profile with the small position gain it always had -- once the
%      offset is dead, a few metres of residual error costs <1 deg of
%      tilt, and nothing needed changing there. The blend is a tanh over
%      +-tBlend so the Riccati RHS stays smooth; tSwitch is auto-detected
%      from sol (|T*| crossing mid-annulus, one switch by G5).
%  (b) ANISOTROPIC position weight: only the LATERAL (x,y) channels are
%      re-scheduled -- BOTH their position weights (qP) and their velocity
%      weights (qV). (CORRECTED, external code review 2026-08-09: this
%      paragraph used to claim "only the LATERAL (x,y) POSITION weights"
%      are re-scheduled and that "ALL velocity weights keep their
%      pre-round-5 values", which the code below contradicts --
%      opts.QA = diag([qP qP 1e-4 qV qV qVelZ 0]) sets QA(4,4)=QA(5,5)=qV
%      from the same pole-placement derivation. The CODE is right and the
%      old comment was wrong: qP alone does not place a pole, it only sets
%      the undamped frequency; qV is what buys the specified damping
%      ratio zetaA, and pricing position without repricing velocity would
%      give phase A an underdamped lateral loop, not the critically
%      damped one this design specifies.) The z-POSITION weight (1e-4)
%      and the vertical VELOCITY weight (qVelZ, set by the separate
%      vertical-loop design below) are common to both phases and are NOT
%      re-scheduled, and R stays the isotropic
%      1e-9*eye(3). This is deliberately more surgical than the
%      anisotropic-R form of the same idea: the vertical channel already
%      demonstrably works (nominal vz tracked the guidance to <0.3 m/s
%      through the whole braking arc at the old weights), so it is left
%      untouched and only the axis that was measurably starved is fed.
%      Weighting position anisotropically in Q is equivalent to pricing
%      the axes differently in R for a diagonal problem, without also
%      re-scaling the vertical loop.
%
% Measured at N=60, dsp.dr0=[50;-30;0] (task-7 report round 5 for the full
% (w,zeta) sweep and the dispersion battery):
%   before: dispersed miss 10.851 m, vtd 6.819 m/s, sat_frac 0.653
%   after : dispersed miss  0.944 m, vtd 1.363 m/s, sat_frac ~0.13
% i.e. both gates (miss<15, vtd<2) pass together with ~14 m and ~0.6 m/s
% of margin, on the same R the round-4 report proved could not satisfy
% both at once -- because R was never the axis that controlled it.
% Nominal is unchanged-clean (landed, miss 0.052 m, vtd 1.480 m/s).
%
% KNOWN REMAINING DEFECT, NOT THIS ROUND'S (measured, orthogonal to the
% above): a pure ALTITUDE dispersion is still tracked far too weakly --
% dr0=[0;0;-50] touches down at 20.4 m/s, dr0=[0;0;+50] arrests 9.2 m up,
% and both numbers are IDENTICAL before and after this round's change, so
% they are a pre-existing defect on a different axis, neither caused nor
% cured here. Same root cause shape (the margin arc's authority is not
% asked for), and a probe (task-7 report round 5) shows scheduling the
% z-position weight up in phase A does move it a long way -- e.g. a
% wz=0.45 rad/s z-loop takes dr0=[0;0;-50] from 20.4 to 0.85 m/s -- but it
% trades into vertical ARRESTS (arriving below the guidance's own descent
% rate is unrecoverable when Tmin exceeds weight), so it needs its own
% design round with the terminal v(z) schedule, not a weight bump. This
% matters for task 8, whose 1-sigma r0 draw includes 50 m of altitude.
%
% REFERENCES:
%   [1] Anderson & Moore, "Optimal Control: Linear Quadratic Methods."
if nargin < 3, opts = struct(); end
% Phase-scheduled state weight. qPos/qVel ratio sets the lateral
% position-error pole (see ADAPTATION 5); the overall qPos level sets the
% commanded correction accel, a = sqrt(qPos/r)/m.
if ~isfield(opts,'R'), opts.R = 1e-9*eye(3); end

%% Nominal interpolants (state: pchip over nodes; thrust: shared
%% per-segment HS quadratic reconstruction, scalar t only -- see the
%% ADAPTATION note above):
Nn = size(sol.X,2) - 1;  h = sol.tf/Nn;
ctrl.xnom = @(t) interp1(sol.t.', sol.X.', clampt(t, sol.tf), 'pchip').';
ctrl.Tnom = @(t) hs_quad_ctrl(clampt(t, sol.tf), sol.U, sol.Um, h, Nn, P.Tmin, P.etaT*P.Tmax);

%% Phase schedule: Tmin-arc -> Tmax-arc switch, auto-detected from sol.
%% tBlend is resolved FIRST because annulus_switch's no-switch fallback is
%% expressed in units of it (external code review, 2026-08-09 -- see that
%% local function's NO-SWITCH FALLBACK note).
if isfield(opts,'tBlend'),  ctrl.tBlend  = opts.tBlend;
else,                       ctrl.tBlend  = 0.4;                     end
if isfield(opts,'tSwitch'), ctrl.tSwitch = opts.tSwitch;
else,   ctrl.tSwitch = annulus_switch(sol, P, ctrl.tBlend);         end

%% Phase-A lateral weights BY POLE PLACEMENT (see ADAPTATION 5 part (a)).
%% Derived from (omegaA, zetaA) and the ACTUAL R and phase-A mass rather
%% than hardcoded, so that an opts.R override re-derives the same closed
%% loop instead of silently detuning phase A (omega scales as r^-1/4, so
%% a 100x R change would otherwise move the pole by 3.2x).
if ~isfield(opts,'omegaA'), opts.omegaA = 0.70; end   % lateral wn [rad/s]
if ~isfield(opts,'zetaA'),  opts.zetaA  = 1.00; end   % lateral damping
inA = sol.t <= ctrl.tSwitch;
if ~any(inA), inA(1) = true; end                      % no-switch fallback
                                                      % guard (ts<0 then,
                                                      % so no node is inA)
mA  = mean(sol.X(7, inA));                            % mass over the margin arc
mB  = mean(sol.X(7, ~inA));                           % mass over the braking arc
if isnan(mB), mB = sol.X(7,end); end

%% VERTICAL VELOCITY LOOP (task-7b co-design). Under altitude indexing the
%% altitude-error channel is identically zero (see sim_closed_loop), so the
%% vertical axis is a PURE FIRST-ORDER velocity loop, vzdot = Tz/m - g. Its
%% scalar ARE gives gain/m = sqrt(qVelZ/r)/m, so a target bandwidth inverts
%% to qVelZ = (m*omegaVz)^2 * r. omegaVz=12 rad/s is ~100x the pre-7b value
%% (measured 0.12 rad/s), which is what the terminal arc actually needs: the
%% descent-rate reference sweeps at ~19.8 m/s^2, so a loop of bandwidth w
%% carries a standing lag of 19.8/w m/s, and 0.12 rad/s means the vehicle is
%% simply never tracking the last 150 m at all.
if ~isfield(opts,'omegaVz'), opts.omegaVz = 12; end   % vertical vel loop [rad/s]
qVelZ = (mB * opts.omegaVz)^2 * opts.R(3,3);

if ~isfield(opts,'QB'), opts.QB = diag([1e-4 1e-4 1e-4 1e-2 1e-2 qVelZ 0]); end
if ~isfield(opts,'QA')
    rLat = opts.R(1,1);                               % lateral control price
    wA   = opts.omegaA;  zA = opts.zetaA;
    % qVel>=0 needs zetaA>=1/sqrt(2) (LQR's own Butterworth floor, see the
    % ADAPTATION 5 part (a) derivation above) -- assert it here rather than
    % let a below-floor zetaA silently produce a negative qV a few lines
    % down (bonus, final-review, 2026-08-09).
    assert(zA >= 1/sqrt(2) - eps, ...
        'tvlqr_design:zetaAFloor', ...
        'opts.zetaA=%.6g is below the LQR Butterworth floor 1/sqrt(2)=%.6g -- the derived qVel would go negative.', ...
        zA, 1/sqrt(2));
    qP   = mA^2 * rLat * wA^4;
    qV   = mA^2 * rLat * wA^2 * (4*zA^2 - 2);         % >=0 needs zetaA>=1/sqrt(2)
    opts.QA = diag([qP qP 1e-4 qV qV qVelZ 0]);
end
%% Qf(6,6) must carry the STEADY-STATE RICCATI VALUE of the vertical
%% velocity channel, or the boundary layer silently undoes the design
%% exactly where it matters. With the original Qf(6,6)=1 the terminal
%% vertical gain collapsed to sqrt(1/r)/m = 1.19 1/s -- and since altitude
%% indexing evaluates K at t*(z) -> tf as z -> 0, the vehicle flew the LAST
%% METRE on a 1.19 1/s loop no matter how the rest was tuned.
%%
%% The right value is P_ss, NOT qVelZ. Q is a running cost (cost per unit
%% state^2 per unit TIME); Qf is a terminal cost (cost per unit state^2).
%% Setting Qf(6,6)=qVelZ mixes the two units. The scalar ARE for this
%% channel, -p^2/(r m^2) + q = 0, gives
%%       P_ss = m*sqrt(qVelZ*r),
%% which reproduces the design bandwidth at tf: K/m = P_ss/(r m^2) =
%% sqrt(qVelZ/r)/m = omegaVz. Using qVelZ instead put P(tf) 12x high
%% (111.29 vs 9.27) and spiked K(3,6)/m from 12.6 to 159 1/s over the last
%% 0.14 s -- a saturation-dominated overshoot layer, not the documented
%% first-order loop. Corrected here (task-7b polish, reviewer catch).
if ~isfield(opts,'Qf')
    PssZ = mB * sqrt(qVelZ * opts.R(3,3));        % ARE steady state, not qVelZ
    opts.Qf = diag([1e-2 1e-2 1e-2 1 1 PssZ 0]);
end
if isfield(opts,'Q'), opts.QA = opts.Q;  opts.QB = opts.Q; end   % legacy
ctrl.Qf = opts.Qf;   % exposed so tests assert P(tf)=Qf without a magic number
QA = opts.QA;  QB = opts.QB;  ts = ctrl.tSwitch;  tb = ctrl.tBlend;
ctrl.Qfun = @(t) (1 - 0.5*(1 + tanh((t - ts)/tb)))*QA ...
                 +    0.5*(1 + tanh((t - ts)/tb)) *QB;

%% Backward Riccati on vec(P), dense output grid:
M     = 4*(Nn+1);
ctrl.tgrid = linspace(0, sol.tf, M);
Rinv  = inv(opts.R);
ric   = @(t, pv) ricrhs(t, pv, ctrl, P, ctrl.Qfun, Rinv);
oo    = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, PV] = ode45(ric, fliplr(ctrl.tgrid), opts.Qf(:), oo);   % integrates tf->0
PV    = flipud(PV);                                          % re-order 0->tf
ctrl.Pt = zeros(7,7,M);  ctrl.K = zeros(3,7,M);
for k = 1:M
    Pk = reshape(PV(k,:), 7, 7);  Pk = (Pk + Pk.')/2;
    [~, ~, Bk] = pdg_dynamics(ctrl.xnom(ctrl.tgrid(k)), ...
                              ctrl.Tnom(ctrl.tgrid(k)), P);
    ctrl.Pt(:,:,k) = Pk;
    ctrl.K(:,:,k)  = Rinv * (Bk.' * Pk);
end
end

function t = clampt(t, tf), t = min(max(t, 0), tf); end

function ts = annulus_switch(sol, P, tBlend)
% ANNULUS_SWITCH  Time the nominal thrust magnitude crosses mid-annulus.
%
% The min-fuel solution is bang-bang on |T*| (certify_pdg's G5 gate counts
% exactly one switch on this trajectory), so a single mid-annulus crossing
% cleanly separates the Tmin coast-down arc from the Tmax braking arc.
%
% TODO (Phase 2): this detects the FIRST UPWARD crossing only. On a
% max-min-max profile (which drag, P.drag.on, could produce) that is the
% SECOND switch, so phase A would wrongly cover the leading Tmax arc.
% Generalize to "last upward crossing" or to a per-arc schedule if G5 ever
% certifies more than one interior switch.
%
% INPUTS:  sol    - collocation solution;  P - booster_params (Tmin,Tmax)
%          tBlend - the caller's tanh blend half-width [s], needed only to
%                   express the no-switch fallback (see below)
% OUTPUTS: ts     - switch time [s]
%
% NO-SWITCH FALLBACK: phase B (the conservative braking-arc weights) over
% the whole flight. An all-Tmax profile has NO margin arc to schedule
% against, so there is nothing to be aggressive with; returning sol.tf
% instead -- as this function did when first written -- would apply the
% aggressive phase-A LATERAL weights across the entire flight INCLUDING
% the braking arc, which is precisely the failure mode the schedule exists
% to prevent.
%
% BUG FIX (external code review, 2026-08-09): the fallback used to return
% ts=0, which does NOT deliver that intent. The caller's schedule is
%   Qfun(t) = (1 - w(t))*QA + w(t)*QB,   w(t) = 0.5*(1 + tanh((t-ts)/tb)),
% so ts=0 gives tanh(0)=0 and w(0)=0.5 at t=0 -- a 50/50 QA/QB blend,
% i.e. HALF the aggressive phase-A lateral gains injected into the first
% moments of flight, exactly where the fallback was supposed to be
% conservative. (The error decays over ~tBlend, so it corrupts the opening
% transient rather than the whole flight, but the opening transient is
% where a large lateral dispersion is corrected.) The fix places the
% switch far enough in the PAST that the tanh saturates: at ts=-10*tBlend,
% w(0) = 0.5*(1+tanh(10)) = 1 - 2.1e-9, i.e. pure QB to nine digits.
% Expressed in units of tBlend rather than as a bare constant so it stays
% saturated for any caller-supplied blend width.
Tmag = sqrt(sum(sol.U.^2, 1));
mid  = 0.5*(P.Tmin + P.etaT*P.Tmax);   % GUIDANCE ceiling, not the engine's
kx   = find(Tmag(1:end-1) < mid & Tmag(2:end) >= mid, 1, 'first');
if isempty(kx)
    ts = -10 * tBlend;  return
end
w  = (mid - Tmag(kx)) / (Tmag(kx+1) - Tmag(kx));
ts = sol.t(kx) + w*(sol.t(kx+1) - sol.t(kx));
end

function pdot = ricrhs(t, pv, ctrl, P, Qfun, Rinv)
Pk = reshape(pv, 7, 7);
[~, A, B] = pdg_dynamics(ctrl.xnom(t), ctrl.Tnom(t), P);
Pd   = -(A.'*Pk + Pk*A - Pk*B*Rinv*B.'*Pk + Qfun(t));
pdot = Pd(:);
end
