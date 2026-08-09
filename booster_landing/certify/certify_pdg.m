function rep = certify_pdg(solC, solV, P, tolScale)
% CERTIFY_PDG  Gates G1-G5 for a PDG solution pair (colloc + convex).
%
% G1 (discrete):     max Hermite-Simpson defect, re-evaluated in plain
%                    MATLAB pdg_dynamics (independent of the CasADi build
%                    the NLP was actually solved with).
% G2 (continuous):   ode45 re-integration of the RECONSTRUCTED colloc
%                    control against the colloc state trajectory --
%                    "defect is not accuracy" (house lesson): G1 only
%                    proves the discrete constraints were satisfied, G2 is
%                    the first real measurement of continuous-time error.
%                    Control reconstruction is the exact per-segment
%                    quadratic Lagrange interpolant through (U_k, Um_k,
%                    U_{k+1}) -- see the "G2 control reconstruction" note
%                    below for why this replaces a global pchip spline.
%                    G2ff (below) is the SECOND HALF of G2, not a separate
%                    gate: same reconstructed control, but asked a
%                    different question -- is it ever outside [Tmin,Tmax]
%                    BETWEEN samples (feedforward feasibility) rather than
%                    how far the flown trajectory ends up from the target
%                    (continuous residual). See G2ff_pass in OUTPUTS below.
% G3 (cross-method): final mass / final time / position-trajectory
%                    agreement between the independently-formulated
%                    collocation (nonconvex annulus) and convexified
%                    (lossless relaxation) solvers. ALL THREE are scored
%                    (the trajectory L-inf row was info-only before the
%                    2026-08-09 external-review fix wave -- see the "G3
%                    trajectory row" note below). The |dmf| threshold is
%                    a MEASUREMENT of Taylor-bound model error (adjudicated
%                    2026-08-08, see the "G3 |dmf| gate" note below), not
%                    an arbitrary agreement tolerance.
% G4 (losslessness): the convex relaxation gap (||u||-sigma) is tight.
% G5 (PMP structure): throttle bang-bang (>=95% of node+midpoint samples
%                    on a bound, <=2 interior switches, max-thrust at
%                    touchdown) and primer-vector alignment (MIDPOINT
%                    thrust parallel to the segment's velocity-costate
%                    direction -- see the "G5 primer TIME BASE" note
%                    below for why midpoint, not node), the two textbook
%                    necessary conditions for a mass-optimal 3-DOF
%                    landing burn (Lawden/PMP).
%
% G2 control reconstruction (CHANGED from the brief's literal global-pchip
% spec, documented in the task-5 fix report): the brief's original G2 flew
% a single pchip spline threaded through ALL nodes+midpoints. That measured
% a mass residual that PLATEAUED around 0.6-0.8 kg for N>=60 and did not
% shrink with further mesh refinement (checked out to N=240) -- exactly
% the profile of a control-REPRESENTATION mismatch, not a real
% continuous-time defect: a global spline does not reconstruct the same
% control HS's own defect equations (Simpson's rule) are built against,
% which is a SEPARATE quadratic per segment, not one spline across segment
% boundaries. Decisive experiment (task-5 fix report): re-flying the SAME
% saved N=60 and a fresh N=240 solution with the per-segment quadratic
% below collapses the residual from ~0.79/0.73 kg to ~0.0001-0.0005 kg at
% BOTH grids (a >1000x drop) -- the floor was an artifact of the pchip
% choice, and is REMOVED (not merely reduced) by using the representation
% HS actually assumes. G2 no longer needs any tolScale accommodation as a
% result (see tolScale doc below).
%
% G5 primer TIME BASE -- the midpoint fix, and the re-tightening it
% earned (external code review, 2026-08-09; SUPERSEDES the "G5 primer
% LOOSENING under drag" note this file carried from task-11, whose
% explanation is now known to have been WRONG).
%
% What the earlier note claimed: the primer angle measured 0.57 deg in
% vacuum and 2.61 deg under drag at the production grid (N=60), the drag
% number did not shrink with mesh refinement the way the vacuum number
% did, and the residual was therefore "most likely" a genuine drag-model
% effect in how the discrete dual approximates the continuous costate --
% left as an open question, with the drag primer gate loosened from 1 deg
% to 10 deg to accommodate it.
%
% ROOT CAUSE (found by the external review, then MEASURED here): the
% comparison was made at the wrong instant. solC.lam_defect(:,k) is the
% multiplier of the Hermite-Simpson defect of SEGMENT k, and the discrete
% stationarity condition it satisfies is the one for that segment's
% MIDPOINT control Um(:,k), not for either end node. Um(:,k) enters the
% NLP through exactly one constraint block -- segment k's defect, via the
% Simpson weight (4h/6)*f(xm,Um) -- plus its own annulus/cone rows, whose
% gradients are parallel to Um itself. So
%     0 = dL/dUm_k = -(4h/6) * B(xm)' * lam_k + (annulus mult) * Um_k
% forces lam_v,k EXACTLY parallel to Um_k (B's velocity rows are I/m and
% its mass row is -T'/(|T| Isp g0), both radial in T). A NODE control
% U(:,k), by contrast, appears in TWO adjacent defect blocks AND inside
% both neighbours' xm = 0.5(x_k+x_k+1) + (h/8)(f_k - f_k+1) terms, so its
% stationarity condition mixes lam_{k-1} and lam_k -- comparing lam_k to
% U(:,k) is an O(h) time-shifted comparison, and O(h) is exactly what
% the old numbers were.
%
% MEASURED (2026-08-09, this file's own solutions, node vs midpoint
% comparison against the same duals):
%     N     case   node_deg    mid_deg
%    20     vac     1.7238     8.5e-07
%    30     vac     1.1424     8.5e-07
%    40     vac     0.8534     8.5e-07
%    60     vac     0.5741     8.5e-07
%    20    drag     8.6275     1.2e-06
%    30    drag     5.1122     1.2e-06
%    40    drag     3.6264     1.2e-06
%    60    drag     2.6051     1.2e-06
% The node column halves as N doubles (O(h), the signature of a time-base
% error); the midpoint column is CONSTANT at 1.2e-06 deg, which is not a
% physical number at all -- it is the double-precision floor of acos near
% 1: acosd(1-eps) = 1.20742e-06 deg. Alignment is exact to floating point,
% in vacuum AND under drag, at every grid. There is no drag effect, no
% discretization artifact, and nothing left open: both were one bug.
%
% RE-TIGHTENED THRESHOLD: 0.01 deg, ONE gate for vacuum and drag alike
% (the drag branch is deleted, not merely reduced -- there is no measured
% drag/vacuum difference left to accommodate). Why 0.01 and not "3x the
% measurement": the measurement IS machine epsilon, so 3x it would be 3x
% numerical noise. What the angle now measures is the KKT stationarity
% residual of the solve, so the honest calibration is against solver
% CONVERGENCE DEPTH, measured here at N=60 by varying opts.tol alone:
%     tol      vac_deg      drag_deg
%    1e-04     0.03704      92.113
%    1e-06     5.09e-05      0.00758
%    1e-08     1.21e-06      4.62e-04
%    1e-09     8.54e-07      2.40e-05
%    1e-11     1.21e-06      1.21e-06
% Both shipped defaults (1e-9 vacuum, 1e-11 drag -- see
% solve_pdg_colloc.m) sit on the acos floor. 0.01 deg therefore carries
% ~4 orders of margin over every shipped configuration and tolerates a
% full 3 orders of degradation in convergence depth, while still FAILING
% on a genuinely unconverged solve (the tol=1e-4 row, which the old 10 deg
% drag gate would have passed at 0.037 deg and failed only in the drag
% column) and on any sign-flip or gross-misalignment regression (90-180
% deg, nowhere near the boundary). Net effect: the gate is 1000x tighter
% than the vacuum semantics it replaces and 1000x tighter than the drag
% loosening it deletes, and it is now a real first-order-optimality
% check rather than a discretization-noise tolerance.
%
% Structure (bang-bang, <=2 switches, max-last) is unaffected and still
% gates every run, drag or vacuum.
%
% G3 trajectory row (external code review, 2026-08-09): G3_traj_Linf --
% the worst position disagreement between the two solvers' trajectories,
% sampled at 200 common times -- was COMPUTED and printed but never
% entered G3_pass, so two solutions could agree on final mass and final
% time (the only scored rows) while taking materially different paths to
% the pad. It is now scored at 5.0 m, alongside |dmf| and |dtf|. The
% number is dimensioned from the mission, not from the measurement: 5.0 m
% is a third of the P.pad_radius = 15 m landing gate, i.e. the two
% independently-formulated solvers are required to agree on WHERE the
% vehicle is to well within the accuracy that decides mission success.
% Measured at the production grid: 0.373 m, a 13x margin.
%
% G5 primer INVALID under a finite pointing cone (final-review fix,
% 2026-08-09): the primer-vector sub-check above (thrust parallel to the
% velocity costate) is a PMP necessary condition for the UNCONSTRAINED
% control problem. The moment P.theta_max_deg is finite, the maximizing T
% is no longer +lambda_v/|lambda_v| whenever that direction falls outside
% the pointing cone -- by KKT, the optimal control there is instead the
% PROJECTION of the unconstrained direction onto the cone boundary, i.e.
% the solution sits ON the cone (T at angle theta_max from vertical) with
% a nonzero multiplier on the cone constraint itself. T is then generally
% NOT parallel to lambda_v by design, not by bug: the naive primer angle
% (this file's Tdir.pdir check) measures the GAP between the constrained
% and unconstrained directions, which is a property of the active
% constraint, not a defect. MEASURED (2026-08-09, P.theta_max_deg=15,
% nominal grid): both solvers independently ride the cone boundary at
% ~15.000 deg from vertical (collocation and convexification agree to
% ~4e-3 deg -- itself a strong cross-method structural check, just not
% the primer's), while the naive primer angle reads ~6.27 deg at the
% same solution -- comfortably over the vacuum <1 deg gate on a
% perfectly-converged optimum, which is exactly the false-fail this guard
% exists to prevent. RE-MEASURED after the midpoint time-base fix
% (2026-08-09, P.theta_max_deg=15, N=60): the corrected MIDPOINT primer
% angle reads 7.322 deg -- essentially identical to the node value
% (7.322 deg) at the same solution, i.e. the cone gap is a REAL
% constrained-vs-unconstrained direction offset that does not go away
% with the time-base fix, unlike the vacuum/drag artifact. The guard is
% therefore still required, and now against a 0.01 deg gate it matters
% far MORE than it did against 1/10 deg. (Midpoint tilt from vertical on
% that solution: max 15.000 deg, min 2.262 deg -- the solution really
% does ride the cone boundary over part of the arc.) FIX: when
% P.theta_max_deg is finite, the primer
% sub-check is not evaluated as a pass/fail gate at all --
% rep.G5_primer_deg is still computed and reported (info only, for
% inspection), rep.G5_primer_mode is set to 'skipped-cone' (the same
% pattern G3/G4 use when solV is []; see print_certify_report.m), and
% rep.G5_pass rests on rep.G5_structOk alone. The 0.01-deg primer bound
% above only ever applies in the P.theta_max_deg=Inf branch -- a cone run
% cannot also be scored against it. Structure (bang-bang,
% <=2 switches, max-last) is unaffected and still gates every run. This
% experiment is otherwise out of scope for this fix wave (P.theta_max_deg
% still defaults to Inf everywhere in this campaign's certified results);
% this guard only prevents a future pointing-cone experiment from either
% false-failing G5 or silently mis-scoring it.
%
% G3 |dmf| gate (ADJUDICATED 2026-08-08, threshold CHANGED from the
% brief's literal 0.1 kg to 1.0 kg, documented in the task-5 fix report):
% G3's |dmf| is reclassified from an "agreement tolerance" (something the
% two solvers should be tuned/refined to close) to a MEASUREMENT of a
% real, understood, bounded physical/modeling effect -- the convex-side
% solver's Taylor-linearized mass bound (see solve_pdg_convex.m's
% mu1/mu2/dz construction) is an approximation, not exact, and its cost is
% the ~0.70 kg offset measured here (about 0.02% of the ~3473 kg of fuel
% burned). This is NOT a discretization artifact: a dedicated refinement
% sweep (fixed tf, Nconv=120/180/240/300, plus a tighter golden-search
% tolTf=0.01) held the gap in a 0.70-0.94 kg band with no downward trend
% -- see the task-5 fix report's "Important 7" section for the full sweep
% and for why Nconv=360's apparent outlier there is IPOPT non-convergence,
% not evidence of a competing optimum (a convex program has none). The
% user adjudicated the number after reviewing that evidence: 1.0 kg is
% the new nominal gate, chosen with headroom over the measured ~0.70-0.94
% kg floor. rep.G3_dmf still reports the raw, unrounded number either way
% -- this only changes what counts as PASS, not what is measured.
%
% INPUTS:
%   solC     - sol struct from solve_pdg_colloc [see that function]
%   solV     - sol struct from solve_pdg_convex, or [] to skip G3/G4
%              (Phase 2 drag runs have no convex twin, since the
%              convexification is only exact in vacuum)
%   P        - booster_params struct
%   tolScale - (optional, def 1) scales ONLY the G3 cross-method
%              agreement tolerances; for COARSE-grid test calls only (see
%              tests/test_certify_nominal.m) -- never used to loosen the
%              nominal-grid (tolScale=1) tolerance in production reports.
%              G2 does NOT use tolScale (see the reconstruction note
%              above: with the corrected control representation its
%              residual is far under gate even at coarse grids, so
%              scaling it is unnecessary and would only hide a real
%              regression). G3's |dmf| gate is now 1.0 kg at tolScale=1
%              (adjudicated 2026-08-08, see the "G3 |dmf| gate" note
%              above) -- a measurement threshold with headroom over the
%              genuine ~0.70-0.94 kg Taylor-bound floor, not merely an
%              agreement tolerance -- so most callers no longer need
%              tolScale>1 for G3 either; it remains available for
%              COARSE-grid calls where |dmf| runs a bit higher still.
%
% OUTPUTS:
%   rep - gate report struct:
%     .tolScale  - the tolScale this report was computed with (so the
%                  struct is self-describing; print_certify_report reads
%                  this to show EFFECTIVE thresholds, not nominal ones)
%     .G1_defect, .G1_pass
%     .G2_pos, .G2_vel, .G2_dm, .G2_pass
%     .G2ff_below_tmin, .G2ff_above_tmax, .G2ff_pass (NEW, task-7 fix
%       report round 4: feedforward-feasibility -- is the flyable
%       reconstructed control ever outside [Tmin,Tmax] BETWEEN
%       nodes+midpoints, dense-sampled; see the G2b note above G2ff_pass)
%     .G3_dmf, .G3_dtf, .G3_traj_Linf, .G3_pass ('skipped' if solV==[])
%       -- all three rows are SCORED (traj_Linf was info-only before the
%       2026-08-09 review fix; see the "G3 trajectory row" note above)
%     .G4_gap, .G4_pass ('skipped' if solV==[])
%     .G5_bound_frac, .G5_switches, .G5_structOk, .G5_primer_deg, .G5_pass
%     .G5_primer_mode - 'scored' (default, P.theta_max_deg=Inf) or
%       'skipped-cone' (P.theta_max_deg finite -- primer_deg is info-only,
%       G5_pass rests on G5_structOk alone; see the "G5 primer INVALID
%       under a finite pointing cone" header note above)
%     .G0_tf_err, .G0_dt_err, .G0_tf_tol, .G0_dt_tol, .G0_pass, .G0_msg
%       - TIME-BASE consistency
%       (G1 uses h=tf/N, G2/G5 read solC.t; see the G0 block below). These
%       were two bare `assert` calls until the 2026-08-09 external review
%       pointed out they contradicted this file's own never-throws
%       contract; they are now scored like any other gate, with the
%       assertion text preserved verbatim in .G0_msg.
%     .all_pass - logical AND of all gates (skipped gates excluded)
%
% This function is report-only and never throws; callers (tests, the
% front door) decide what a FAIL means. This is now literally true: the
% two time-base assertions that used to violate it are gate G0 (see
% above). A malformed or failed-solver solC produces a report with failed
% gates, not an exception with no report at all.
%
% REFERENCES:
%   [1] Acikmese & Ploen, "Convex Programming Approach to Powered Descent
%       Guidance for Mars Landing," JGCD 2007. (lossless convexification)
%   [2] Lawden, "Optimal Trajectories for Space Navigation," 1963.
%       (primer vector theory / PMP bang-bang structure)
if nargin < 4, tolScale = 1; end
rep.tolScale = tolScale;
% isfield guard (task-11 close-out review, parity with solve_pdg_colloc.m's
% own "if isfield(P,'drag') && P.drag.on" tol-default check): a bare
% P.drag.on would throw for any caller whose P predates the drag field
% instead of just treating it as vacuum.
rep.drag_on  = isfield(P,'drag') && P.drag.on;
                            % Reported so the struct is self-describing.
                            % It no longer selects a primer THRESHOLD: the
                            % drag loosening was deleted once the primer
                            % time-base bug was fixed (see the "G5 primer
                            % TIME BASE" header note) -- vacuum and drag
                            % are now scored against the same 0.01 deg.

%% GUIDANCE thrust ceiling (task-7b, P.etaT). Every gate below certifies a
%% GUIDANCE solution, so its upper annulus bound is the de-rated etaT*Tmax,
%% not the engine's full Tmax: the control reconstruction clamp (G2/G2ff),
%% the G2ff over-ceiling check, the G4 tightness scale and G5's "on the
%% upper bound" detector all key off TmaxG. The tracker/sim keep the full
%% annulus -- see sim_closed_loop.m's allocate_thrust.
TmaxG = P.etaT * P.Tmax;
rep.TmaxG = TmaxG;

%% G1: HS defects re-evaluated in plain MATLAB (independent of CasADi):
N = size(solC.X,2) - 1;  h = solC.tf / N;  dmax = 0;
for k = 1:N
    fk  = pdg_dynamics(solC.X(:,k),   solC.U(:,k),   P);
    fk1 = pdg_dynamics(solC.X(:,k+1), solC.U(:,k+1), P);
    xm  = 0.5*(solC.X(:,k)+solC.X(:,k+1)) + (h/8)*(fk - fk1);
    fm  = pdg_dynamics(xm, solC.Um(:,k), P);
    d   = solC.X(:,k+1) - solC.X(:,k) - (h/6)*(fk + 4*fm + fk1);
    dmax = max(dmax, max(abs(d)));
end
rep.G1_defect = dmax;   rep.G1_pass = dmax < 1e-6;

%% G0: time-base consistency. G1 above uses h=tf/N; G2/G5 below use solC.t
%% directly. If a future solver's solC.t ever drifted from a uniform
%% h=tf/N grid, G1 and G2/G5 would silently disagree by an h-INDEPENDENT
%% offset invisible to G1 (which never reads solC.t at all) -- exactly
%% the failure mode this campaign burned time on before ("defect is not
%% accuracy"). Check it explicitly rather than assume it.
%%
%% SCORED, NOT ASSERTED (external code review, 2026-08-09): these were two
%% bare `assert` calls, which contradicted this function's own documented
%% "report-only and never throws" contract -- and did so in the worst
%% direction: a malformed or failed-solver solC (exactly the input a
%% caller most needs a report for) threw instead of producing one, so the
%% other four gates' diagnostics were lost with it. They are now gate G0.
%% The failure MESSAGES are preserved verbatim in rep.G0_msg, so nothing
%% a caller could previously read in the exception text is lost; it is
%% simply delivered as data instead of as a throw. Downstream gates still
%% run on a G0 failure (their numbers are then suspect, which is the
%% point of reporting rather than hiding them), and rep.all_pass includes
%% G0 so no caller can pass certification on an inconsistent time base.
% Thresholds are RELATIVE (1e-9 * max(1,|scale|)), so they are stored in
% the report rather than left for print_certify_report to guess -- the
% same discipline the tolScale/EFFECTIVE-threshold machinery exists for.
rep.G0_tf_tol = 1e-9 * max(1, abs(solC.tf));
rep.G0_dt_tol = 1e-9 * max(1, abs(h));
rep.G0_tf_err = abs(solC.t(end) - solC.tf);
rep.G0_dt_err = max(abs(diff(solC.t) - h));
g0msg = {};
if ~(rep.G0_tf_err < rep.G0_tf_tol)
    g0msg{end+1} = sprintf(['solC.t(end)=%.10g does not match solC.tf=%.10g -- ' ...
        'G1 (h=tf/N) and G2/G5 (built from solC.t) would silently disagree.'], ...
        solC.t(end), solC.tf);
end
if ~(rep.G0_dt_err < rep.G0_dt_tol)
    g0msg{end+1} = sprintf(['solC.t is not uniformly spaced at h=tf/N=%.10g ' ...
        '(max spacing error %.3g) -- G1/G2/G5 time bases would silently disagree.'], ...
        h, rep.G0_dt_err);
end
rep.G0_pass = isempty(g0msg);
rep.G0_msg  = strjoin(g0msg, ' | ');

%% Node+midpoint control samples in chronological order (node_1, mid_1,
%% node_2, mid_2, ..., node_{N+1}), used by G5's bang-bang statistics
%% below. mid_k sits strictly between node_k and node_{k+1} in time, so
%% this interleaving is already time-ordered without needing a sort.
TU = zeros(3, 2*N + 1);
TU(:,1:2:end) = solC.U;  TU(:,2:2:end) = solC.Um;

%% G2: continuous residual -- fly the RECONSTRUCTED control with ode45
%% ("defect is not accuracy"). See the reconstruction note in the header:
%% exact per-segment quadratic Lagrange through (U_k, Um_k, U_{k+1}), the
%% control representation Simpson's rule (hence the HS defects above) is
%% actually built against, not a global spline across segment boundaries.
odef = @(tt, xx) pdg_dynamics(xx, hs_quad_ctrl(tt, solC.U, solC.Um, h, N, P.Tmin, TmaxG), P);
% Integration structure from the shared engine (oclib move 2, span mode =
% this gate's single-call ode45 convention); dynamics + the
% annulus-feasible control reconstruction stay THIS campaign's closures.
zG2 = oc.fly_control(solC.X(:,1), [0 solC.tf], odef, ...
    struct('mode','span', 'solver',@ode45, 'RelTol',1e-10, 'AbsTol',1e-10));
ef = zG2 - solC.X(:,end);
rep.G2_pos = sqrt(sum(ef(1:3).^2));  rep.G2_vel = sqrt(sum(ef(4:6).^2));
rep.G2_dm  = abs(ef(7));
rep.G2_pass = rep.G2_pos < 1 && rep.G2_vel < 0.1 && rep.G2_dm < 0.5;

%% G2b (task-7 fix report round 4, 2026-08-08, NEW): feedforward
%% FEASIBILITY -- is the FLYABLE reconstructed control (the same
%% hs_quad_ctrl call G2 just integrated, and the exact function Task 6's
%% TVLQR feedforward and Task 7's closed-loop truth sim also call) ever
%% outside [Tmin,Tmax] BETWEEN nodes+midpoints, not just at the sampled
%% points G5's bound-fraction check below already covers? Root-cause gate
%% for a real bug: before hs_quad_ctrl's round-4 direction/magnitude split
%% (see that function's own ADAPTATION note), the single-vector
%% reconstruction dipped up to 18% below Tmin over 11.7% of the flight at
%% the N=30 test grid -- invisible to G1 (discrete defects only), G2's OWN
%% pass/fail (a feedforward that under-thrusts for a while and
%% over-thrusts elsewhere can still land close in aggregate), and G5
%% (node+midpoint samples only). The closed-loop sim's saturating clamp
%% turned that invisible feedforward defect into a real +0.6-1m open-loop
%% altitude bias that two rounds of TVLQR weight tuning were unknowingly
%% fighting. Densely sampled (20 points/segment, not just nodes+
%% midpoints) so a dip strictly BETWEEN samples cannot hide. Gated at
%% every grid (not just nominal): the direction/magnitude split makes a
%% violation impossible by construction regardless of resolution, so a
%% real regression here should never require a coarse-grid exemption.
tdense = linspace(0, solC.tf, 20*N + 1);
TmagD  = zeros(size(tdense));
for kk = 1:numel(tdense)
    Tvk = hs_quad_ctrl(tdense(kk), solC.U, solC.Um, h, N, P.Tmin, TmaxG);
    TmagD(kk) = sqrt(sum(Tvk.^2));
end
rep.G2ff_below_tmin = max(max(0, P.Tmin - TmagD) / P.Tmin);
rep.G2ff_above_tmax = max(max(0, TmagD - TmaxG) / TmaxG);
rep.G2ff_pass = rep.G2ff_below_tmin < 1e-6 && rep.G2ff_above_tmax < 1e-6;

%% G3/G4: cross-method agreement + losslessness (skip if no convex twin):
if isempty(solV)
    rep.G3_pass = 'skipped';  rep.G4_pass = 'skipped';
else
    rep.G3_dmf = abs(solC.mf - solV.mf);
    rep.G3_dtf = abs(solC.tf - solV.tf);
    tq  = linspace(0, min(solC.tf, solV.tf), 200);
    rC  = interp1(solC.t.', solC.X(1:3,:).', tq.', 'pchip');
    rV  = interp1(solV.t.', solV.X(1:3,:).', tq.', 'pchip');
    rep.G3_traj_Linf = max(sqrt(sum((rC - rV).^2, 2)));
    % |dmf| threshold 1.0 kg (was 0.1 kg, adjudicated 2026-08-08 -- see
    % the "G3 |dmf| gate" note in this function's header): a measurement
    % gate over the genuine, non-shrinking ~0.70-0.94 kg Taylor mass-bound
    % model error, not an agreement tolerance the solvers should close.
    % traj_Linf is SCORED at 5.0 m (external review, 2026-08-09 -- it was
    % computed and printed but never gated, so two solvers could agree on
    % mf and tf while flying materially different paths). 5.0 m is a third
    % of P.pad_radius=15 m: cross-method agreement on WHERE the vehicle is
    % must be well inside the accuracy that decides mission success. See
    % the "G3 trajectory row" header note.
    rep.G3_pass = rep.G3_dmf < 1.0*tolScale && rep.G3_dtf < 0.2*tolScale ...
                  && rep.G3_traj_Linf < 5.0*tolScale;
    rep.G4_gap  = solV.lossless_gap;
    rep.G4_pass = rep.G4_gap < 1e-4 * TmaxG / P.m0;
end

%% G5: PMP structure. Throttle bang-bang + primer alignment.
%% Bound/switch statistics run over NODES *AND* MIDPOINTS (TU, built
%% above) -- solC.Um are independent NLP decision variables and G2 flies
%% them too, so a mid-annulus MIDPOINT is a real, undetected chatter
%% signal if only nodes are checked (a node-only check could show
%% bound_frac=1 while a midpoint sits strictly inside the annulus).
Tmag = sqrt(sum(TU.^2, 1));
onLo = Tmag < P.Tmin * 1.001;   onHi = Tmag > TmaxG * 0.999;
rep.G5_bound_frac = mean(onLo | onHi);
segs = diff([onHi(1), onHi]);                  % max<->min transitions
rep.G5_switches = sum(segs ~= 0);
% Structure check (ADAPTATION from the brief's literal "onHi(1) && onHi(end)"
% [max-first AND max-last], documented in the task-5 report): the measured
% nominal solution is min-first / max-last (onHi(1)=false, single switch
% Tmin->Tmax at t~7.06 s of tf~15.69 s, held to touchdown), not max-first.
% This is genuine physics, not a bug: task-3's own report established that
% P.Tmin (338 kN) already exceeds the hover-equivalent thrust at any mass
% in [mdry,m0] (m0*g0=294.2 kN, mdry*g0=251.1 kN), so the vehicle cannot
% loiter even at minimum throttle -- coasting at Tmin early, then a single
% hard Tmax brake just before the pad, is the fuel-optimal bang-bang
% structure for THIS entry state, and PMP does not mandate a specific
% opening throttle (only: bang-bang, no singular arc, and -- for a landing
% burn that must arrest velocity before touchdown -- max thrust at the
% very end). We therefore require max-LAST (a physical necessity: some
% peak-thrust braking is needed right before touchdown) but not max-FIRST.
% onHi(end) here is TU's last column, i.e. solC.U(:,end) (TU interleaves
% node,mid,...,node so its last of 2N+1 columns is the final node) -- the
% same physical quantity the brief's original node-only check used.
rep.G5_structOk = rep.G5_bound_frac >= 0.95 && rep.G5_switches <= 2 && onHi(end);
% Primer: velocity costate from defect duals; thrust must align with the
% costate direction. SIGN (ADAPTATION, documented per the brief's own
% instruction: "if G5_primer_deg comes out near 180 deg ... flip pdir sign
% ONCE"): the brief's literal pdir = -lamv/|lamv| gave primer_deg ~179.7-
% 179.9 deg (thrust anti-parallel to -lam_v). Flipping to pdir=+lamv/|lamv|
% gives ~0.15-1.2 deg (thrust parallel to +lam_v) -- so for the sign
% convention solve_pdg_colloc's opti.lam_g extraction actually returns
% (house lesson: opti.lam_g, never opti.dual -- see MEMORY
% casadi-optidual-sign-bug), the correctly-aligned primer direction is
% +lam_v, not -lam_v. Applied ONCE here, no abs() anywhere in this file.
%
% TIME BASE (external code review, 2026-08-09): the comparison is against
% the segment MIDPOINT control solC.Um(:,k), NOT the node control
% solC.U(:,k). solC.lam_defect(:,k) is segment k's defect multiplier, and
% the only control whose discrete stationarity condition involves that
% multiplier ALONE is Um(:,k) -- a node control appears in two adjacent
% defect blocks and in both neighbours' xm terms, so scoring lam_k
% against U(:,k) is an O(h) time-shifted comparison. That single line was
% the entire source of this campaign's "vacuum discretization artifact"
% (0.57 deg) and "drag primer degradation" (2.61 deg); both collapse to
% the acos machine-precision floor (1.2e-6 deg) under the corrected
% comparison, at every grid. Full measurement tables and the resulting
% threshold re-tightening are in the "G5 primer TIME BASE" header note.
% Um is 3xN, exactly matching lam_defect's N columns -- no truncation
% needed (the old node form had to drop the last node).
% STATION ASSOCIATION + SIGN come from the shared library home
% (oclib/+oc/duals_to_costates) since 2026-08-09: this gate independently
% re-invented, and independently re-fixed, the Hermite-Simpson
% midpoint-vs-node rule that orbit_transfer's covector mapping already
% owned -- one home per subtle rule ends that failure mode. The library
% resolves the global dual sign by its primer-vs-control vote (this
% campaign's duals carry +lam_v parallel to thrust, the opposite of the
% orbit convention; the vote absorbs the difference), and returns
% costates whose PRIMER -lam_v/|lam_v| aligns with thrust. The angle
% computed below is therefore convention-free.
if isempty(which('oc.duals_to_costates'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'oclib'));
end
[lamMap, ~, dgMap] = oc.duals_to_costates(struct( ...
    'scheme', 'hermite-simpson', 'mu', solC.lam_defect, ...
    'tNodes', solC.t, 'uDir', solC.Um));
Tdir = solC.Um ./ sqrt(sum(solC.Um.^2,1));     % 3xN, MIDPOINT directions
pdir = -lamMap(4:6,:) ./ max(sqrt(sum(lamMap(4:6,:).^2,1)), 1e-30);
cosang = sum(Tdir .* pdir, 1);
rep.G5_primer_deg = max(acosd(min(1, max(-1, cosang))));
rep.G5_voteMargin = dgMap.voteMargin;          % info: 1.00 = unanimous
% Primer threshold RE-TIGHTENED to 0.01 deg for vacuum AND drag alike
% (external code review, 2026-08-09) -- the task-11 drag loosening to
% 10 deg is deleted, because the drag/vacuum difference it accommodated
% was the time-base bug fixed just above, not a drag-model effect. See
% the "G5 primer TIME BASE" header note for the two measurement tables
% (grid sweep and solver-tolerance sweep) behind the number. Structure
% (bang-bang + max-last) still gates every run regardless.
% CONE GUARD (final-review fix, 2026-08-09): under a finite pointing cone
% the primer sub-check is invalid (see the "G5 primer INVALID under a
% finite pointing cone" header note above) -- report primer_deg as
% info-only and gate on structure alone, the same 'skipped'-style pattern
% G3/G4 use above when there is no convex twin. The 0.01 deg bound only
% ever applies in the uncapped (P.theta_max_deg=Inf) branch, since a cone
% run cannot also be scored against it -- and it matters MORE now that the
% bound is 0.01 rather than 1/10 deg, because the cone gap (7.32 deg
% measured) is far above it.
if isfinite(P.theta_max_deg)
    rep.G5_primer_mode = 'skipped-cone';
    rep.G5_pass = rep.G5_structOk;
else
    rep.G5_primer_mode = 'scored';
    rep.G5_pass = rep.G5_structOk && rep.G5_primer_deg < 0.01;
end

%% Verdict (G0 included -- see the G0 block for why it is a gate and not
%% an assertion):
gates = [rep.G0_pass, rep.G1_pass, rep.G2_pass, rep.G2ff_pass, ...
         isequal(rep.G3_pass,true) || isequal(rep.G3_pass,'skipped'), ...
         isequal(rep.G4_pass,true) || isequal(rep.G4_pass,'skipped'), ...
         rep.G5_pass];
rep.all_pass = all(gates);
end

% NOTE (task-7 fix report, 2026-08-08): hs_quad_ctrl was promoted from a
% local function here to lib/hs_quad_ctrl.m (byte-identical logic) so
% Task 6's tvlqr_design.m could build ctrl.Tnom from the same per-segment
% quadratic reconstruction G2 above certifies, instead of a global pchip
% spline. test_certify_nominal.m re-run clean after the move (both blocks
% still all_pass) -- see that test's output in the task-7 fix report for
% the regression check.

% NOTE (deviation from the brief's "nested function" suggestion, documented
% in the task-5 report): print_certify_report is called directly by
% tests/test_certify_nominal.m and (Task 10) the front door -- both outside
% this file. MATLAB local functions are only visible within their own
% function file, so a copy nested here would NOT be callable from either
% caller; the brief's own fallback ("expose it by making
% certify/print_certify_report.m a 3-line wrapper if needed elsewhere")
% applies, and it IS needed elsewhere. Rather than maintain two copies (one
% nested/dead, one real) the single real implementation lives in its own
% file, certify/print_certify_report.m, right beside this one.
