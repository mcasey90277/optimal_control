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
%                    (lossless relaxation) solvers. The |dmf| threshold is
%                    a MEASUREMENT of Taylor-bound model error (adjudicated
%                    2026-08-08, see the "G3 |dmf| gate" note below), not
%                    an arbitrary agreement tolerance.
% G4 (losslessness): the convex relaxation gap (||u||-sigma) is tight.
% G5 (PMP structure): throttle bang-bang (>=95% of node+midpoint samples
%                    on a bound, <=2 interior switches, max-thrust at
%                    touchdown) and primer-vector alignment (thrust
%                    parallel to the velocity-costate direction), the two
%                    textbook necessary conditions for a mass-optimal
%                    3-DOF landing burn (Lawden/PMP).
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
% G5 primer LOOSENING under drag (task-11 fix, Phase 2, 2026-08-09;
% CORRECTED in the task-11 close-out review round, 2026-08-09): documented
% per the task-11 brief step 2's own fallback instruction, "if primer
% alignment fails ONLY in the drag case, relax G5 for drag runs to the
% bang-bang structure check alone." NOTE this is a DIFFERENT brief
% instruction from the primer SIGN-flip rule applied below at the pdir
% computation ("if G5_primer_deg comes out near 180 deg ... flip pdir sign
% ONCE") -- an earlier version of this note conflated the two quotes as
% if one were the "companion" of the other; that framing was wrong and is
% fixed here. What "relax" means was ALSO corrected: the first pass at
% this fix set rep.G5_pass = rep.G5_structOk alone under drag, i.e.
% REMOVED the primer check entirely rather than loosening it. Review
% caught this as a real regression risk: with the primer check gone, a
% future bug that flips pdir's sign (or otherwise sends primer_deg toward
% ~90-180 deg) would silently PASS G5 under drag -- exactly the class of
% bug the sign-flip fix earlier in this file exists to catch, and exactly
% the kind of "loosen a gate until it's not a gate" failure this campaign
% has otherwise been careful to avoid (see the G3 |dmf| gate note below for
% the house standard: a relaxed threshold must still be a MEASUREMENT with
% headroom, not an escape hatch). FIX: rep.G5_pass now requires
% rep.G5_structOk AND rep.G5_primer_deg < 10 (not < 1) under drag -- the
% detector stays, only the threshold loosens. MEASURED, not assumed: the
% same coarse grid (N=30) that gives a vacuum primer_deg of 1.14 deg (a
% known, PRE-EXISTING discretization artifact -- test_certify_nominal.m's
% own note: N>=40 is needed to clear <1 deg in vacuum, N=20 measures ~1.7
% deg) gives 5.11 deg under drag; at the PRODUCTION grid (N=60, where
% vacuum's own primer clears comfortably under 1 deg) drag still measures
% 2.61 deg -- i.e. finer meshing does NOT cure it the way it cures the
% vacuum artifact, so this is a genuine drag-model effect, not (only) a
% discretization one. Analytically the direction condition should be
% UNCHANGED by drag: the Hamiltonian's control-dependent term,
% lambda_v.T/m - lambda_m|T|/(Isp g0), does not involve aD (drag enters
% only through the velocity/altitude/mass BLOCK of the dynamics Jacobian
% A, never through B = d(xdot)/dT), so the maximizing T direction is still
% +lambda_v/|lambda_v| regardless of P.drag.on. The measured degradation
% is therefore most likely in how well the discrete lam_defect dual
% approximates the CONTINUOUS costate once the defect residual's own
% curvature grows under the added |v|v/m nonlinearity (drag couples
% velocity components through vmag = sqrt(sum(v.^2)) in a way vacuum's
% purely-linear-in-v dynamics do not) -- a genuine open question, not
% resolved here; left as this file's flag for the min-fuel paper's future
% work list rather than force-closed by a grid search that was not run to
% convergence in this task. The 10 deg threshold was chosen for margin,
% not tightness: it sits >3x above the measured 2.61 deg production-grid
% value (comfortable headroom against normal run-to-run noise) while still
% catching a gross misalignment or sign-flip regression, which would show
% up in the 90-180 deg range, nowhere near the pass/fail boundary.
% Structure (bang-bang, <=2 switches, max-last) is unaffected and still
% gates every run, drag or vacuum.
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
% the primer's), while the OLD naive primer angle reads ~6.27 deg at the
% same solution -- comfortably over the vacuum <1 deg gate on a
% perfectly-converged optimum, which is exactly the false-fail this guard
% exists to prevent. FIX: when P.theta_max_deg is finite, the primer
% sub-check is not evaluated as a pass/fail gate at all --
% rep.G5_primer_deg is still computed and reported (info only, for
% inspection), rep.G5_primer_mode is set to 'skipped-cone' (the same
% pattern G3/G4 use when solV is []; see print_certify_report.m), and
% rep.G5_pass rests on rep.G5_structOk alone. The drag-loosened 10-deg
% primer bound above only ever applies in the P.theta_max_deg=Inf branch
% -- a cone run cannot also be scored against it. Structure (bang-bang,
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
%     .G4_gap, .G4_pass ('skipped' if solV==[])
%     .G5_bound_frac, .G5_switches, .G5_structOk, .G5_primer_deg, .G5_pass
%     .G5_primer_mode - 'scored' (default, P.theta_max_deg=Inf) or
%       'skipped-cone' (P.theta_max_deg finite -- primer_deg is info-only,
%       G5_pass rests on G5_structOk alone; see the "G5 primer INVALID
%       under a finite pointing cone" header note above)
%     .all_pass - logical AND of all gates (skipped gates excluded)
%
% This function is report-only and never throws; callers (tests, the
% front door) decide what a FAIL means.
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
                            % print_certify_report reads this to print the
                            % G5 primer row against the LOOSENED (10 deg,
                            % not 1 deg) threshold under drag -- see the
                            % "G5 primer LOOSENING under drag" header note
                            % above.

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

%% Time-base consistency check: G1 uses h=tf/N; G2/G5 below use solC.t
%% directly. If a future solver's solC.t ever drifted from a uniform
%% h=tf/N grid, G1 and G2/G5 would silently disagree by an h-INDEPENDENT
%% offset invisible to G1 (which never reads solC.t at all) -- exactly
%% the failure mode this campaign burned time on before ("defect is not
%% accuracy"). Assert it explicitly rather than assume it.
assert(abs(solC.t(end) - solC.tf) < 1e-9 * max(1, abs(solC.tf)), ...
    'certify_pdg:timebase', ...
    'solC.t(end)=%.10g does not match solC.tf=%.10g -- G1 (h=tf/N) and G2/G5 (built from solC.t) would silently disagree.', ...
    solC.t(end), solC.tf);
assert(max(abs(diff(solC.t) - h)) < 1e-9 * max(1, abs(h)), ...
    'certify_pdg:timebase', ...
    'solC.t is not uniformly spaced at h=tf/N=%.10g (max spacing error %.3g) -- G1/G2/G5 time bases would silently disagree.', ...
    h, max(abs(diff(solC.t) - h)));

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
oo   = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~, XX] = ode45(odef, [0 solC.tf], solC.X(:,1), oo);
ef = XX(end,:).' - solC.X(:,end);
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
    rep.G3_pass = rep.G3_dmf < 1.0*tolScale && rep.G3_dtf < 0.2*tolScale;
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
% (Primer is evaluated at NODES only -- solC.lam_defect has no midpoint
% dual, so there is nothing to add here the way G5's bound/switch check
% added Um above.)
lamv = solC.lam_defect(4:6, :);                % 3xN, node-adjacent duals
Tdir = solC.U(:,1:end-1) ./ sqrt(sum(solC.U(:,1:end-1).^2,1));
pdir = lamv ./ max(sqrt(sum(lamv.^2,1)), 1e-30);
cosang = sum(Tdir .* pdir, 1);
rep.G5_primer_deg = max(acosd(min(1, max(-1, cosang))));
% Primer threshold LOOSENED (not removed) under drag (task-11, Phase 2;
% corrected in the task-11 close-out review round) -- see the "G5 primer
% LOOSENING under drag" header note above for the measured evidence (fails
% ONLY in the drag case at the production grid, does not shrink with mesh
% refinement the way the vacuum discretization artifact does) and for why
% dropping the check entirely (an earlier version of this fix) was a real
% regression risk. Structure (bang-bang + max-last) still gates every run
% regardless.
% CONE GUARD (final-review fix, 2026-08-09): under a finite pointing cone
% the primer sub-check is invalid (see the "G5 primer INVALID under a
% finite pointing cone" header note above) -- report primer_deg as
% info-only and gate on structure alone, the same 'skipped'-style pattern
% G3/G4 use above when there is no convex twin. The drag-loosened bound
% only ever applies in the uncapped (P.theta_max_deg=Inf) branch, since a
% cone run cannot also be scored against it.
if isfinite(P.theta_max_deg)
    rep.G5_primer_mode = 'skipped-cone';
    rep.G5_pass = rep.G5_structOk;
elseif isfield(P,'drag') && P.drag.on
    rep.G5_primer_mode = 'scored';
    rep.G5_pass = rep.G5_structOk && rep.G5_primer_deg < 10;
else
    rep.G5_primer_mode = 'scored';
    rep.G5_pass = rep.G5_structOk && rep.G5_primer_deg < 1;
end

%% Verdict:
gates = [rep.G1_pass, rep.G2_pass, rep.G2ff_pass, isequal(rep.G3_pass,true) || ...
         isequal(rep.G3_pass,'skipped'), isequal(rep.G4_pass,true) || ...
         isequal(rep.G4_pass,'skipped'), rep.G5_pass];
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
