# Consult prompt — IFS min-fuel indirect solve (for GPT-5.6 sol)

> Paste everything below the line into GPT-5.6 sol. Self-contained.

---

I'm an engineer working on a low-thrust minimum-fuel trajectory optimization
problem and I've hit a well-characterized numerical wall on the **indirect
(costate) solve**. I want you to (1) tell me if my diagnosis is wrong, and (2)
recommend the single best next technique — with specifics, not a survey. Be
skeptical; if I've mis-diagnosed the failure, say so and why.

## Problem

CR3BP (Earth–Moon, mu* = 0.012150585609624), continuous low thrust.
GTO (350 × 35786 km, argument of perigee −25°) → a south-pole "tulip" orbit
point. m0 = 15 kg, T_max = 25 mN, Isp = 2100 s. The transfer is a **~40-rev
spiral**; minimum time-of-flight is t_f,min = 6.290694 nondimensional
(27.88 days). Minimum-**fuel** (mass-optimal) solution is bang-bang with roughly
10–47 thrust switches depending on the allotted t_f.

I use **Sundman regularization**: independent variable is τ with
dt/dτ = κ = r1^1.5 (r1 = distance from the primary), and physical time t is
carried as an extra state; τ_f is fixed. This tames the 40-rev perigee
sensitivity.

I already have a **working direct solver** (direct collocation + a PMP-steered
mesh-refinement stage) that produces good bang-bang solutions across a band of
t_f. What does NOT yet work is the indirect finishing solve below.

## The indirect method ("IFS")

Take the direct bang-bang solution's **switch structure as fixed** and solve the
first-order PMP conditions by **multiple shooting in Sundman τ with explicit
switch nodes** — each switch is its own shooting node, so there is **no
saltation matrix** (a switch's sensitivity is just the ordinary endpoint
sensitivity of its two neighboring arcs). Each arc runs at a *known constant*
throttle (u=1 burn / u=0 coast), so there is **no smoothing parameter ε** and no
1/ε layer — this deliberately removes the crawl that killed a prior
ε-smoothed indirect attempt.

- Augmented state (16-dim): [r(3); v(3); m; t; λr(3); λv(3); λm; λt].
- Switching function S = 1 − ‖λv‖·c/m − λm  (burn where S<0, coast where S>0;
  c = exhaust-velocity constant). Primer vector α = −λv/‖λv‖.
- Terminal BCs (square system): r,v rendezvous + λm(τ_f)=0 transversality,
  t_f fixed.
- Unknowns: initial costates λ0 plus the interior node states and switch times.
- Residual solved by a scaled Gauss–Newton: two-sided Jacobian equilibration +
  truncated/rank-revealing SVD step + Levenberg fallback + α-floor line search +
  adaptive truncation continuation. Complex-step Jacobians (h = 1e-20).
- Costate seed comes from the direct solution via a **dual→costate map** (recover
  λ from the NLP's KKT multipliers on the dynamics defects).

## What I've found (the wall)

**1. Cold dual-map seed never converges, at ANY t_f in the good band.**
The solver *descends* the cold seed (residual ‖R‖ 1.96 → ~0.43) but **floors at
~0.3–0.47 and stalls** — identically at t_f factors 1.12× (k=10 switches),
1.14× (k=24), and 1.25× (k=47), i.e. even for well-separated switches with no
terminal cluster. So the wall is **not** switch clustering; it appears to be the
**small convergence basin of 40-rev shooting seeded from the direct (KKT-dual)
costates**. The dual→costate recovery is also erratic: at t_f=1.85× the seed
residual blows up to 1.2e5 (a scaling failure in the dual map itself).

**2. Min-time all-burn anchor DOES converge — but it's near-singular.**
As a continuation anchor I built the k=0 (no-switch, all-burn) **minimum-time**
solution: min-time indirect costates + λ_t0 = −H_t(0), with τ_f found by
integrating the hard burn to t = t_f,min. This is my first end-to-end rendezvous
convergence: ‖R‖ → 2.3e-7, S<0 everywhere (all-burn confirmed). But ‖R‖ floors
at ~2e-7 because the min-time point is genuinely near-singular (independent
fsolve on the min-time reference also flags "locally singular").

**3. Continuation in t_f off the min-time anchor fails — I believe it's a fold.**
Naive t_f-stepping: even a **0.1%** t_f increase throws the seed residual to
0.30 and the near-singular Jacobian floors the re-solve at ~0.02–0.09 without
converging; max S moves *away* from 0. So I built **pseudo-arclength
continuation** (Keller): state x = [λ0(8); factor], arclength constraint,
predictor + Newton corrector on the 9×9 [R; arclength], complex-step dR/dx,
scaled metric. The extended 8×9 Jacobian at the anchor is **full rank 8**
(singular values 1.7e6 … 8.3e-4, cond ~2e9) — a clean 1-D tangent, so the fold
looks well-posed. **But**: across a 44-step march the **factor stays pinned at
1.0000 to 4+ decimals** while λ0 changes a lot, and max S **wanders
non-monotonically** (−1.48 → −3.23 → … → +0.04) instead of rising monotonically.
My reading: the branch is **effectively vertical** in (t_f, λ0) at min-time AND
carries a **near-null costate gauge** (smallest scaled singular value 8.3e-4);
the corrector — forced to a loose tolerance (‖R‖<1e-5) because the min-time seed
is only ~1e-6 accurate — lets λ0 drift along the gauge at fixed t_f. So
pseudo-arclength is mechanically correct but cannot advance t_f off this
pathological anchor. The apparent "switch birth" at factor=1.0000 is a gauge
artifact, not a physical band result.

## Questions (please be specific and concrete)

1. **Is the fold + costate-gauge diagnosis right?** Given the min-time point is a
   time-optimal solution (H≡0, free-final-time transversality), is a
   near-degenerate costate direction *expected* there, and is that the true
   obstruction — or am I mis-attributing a *seed-accuracy* problem (my ~1e-6
   min-time costates) to a *geometric* fold that isn't really there?

2. **Anchor choice.** My direct side uses **minimum-energy** as its homotopy root
   (Bertrand–Épénoy energy→fuel), not min-time. Would anchoring the indirect
   continuation at a **min-energy** solution (smooth, non-bang-bang, non-singular)
   and homotoping the throttle toward bang-bang be strictly better than my
   min-time anchor? If so, what's the cleanest continuation parameter and how do
   I keep the switch structure from collapsing as ε→0 (my prior ε-smoothed
   attempt slid to the many-switch global basin)?

3. **Gauge regularization.** If I keep the min-time anchor, what's the standard,
   minimal fix for the near-null λ0 direction — a normalization/phase condition
   (‖λ0‖=1 pinned, or a Poincaré-style phase constraint), truncating the null
   direction, or folding the gauge into the arclength system? Concretely, how
   should the augmented [R; arclength] system be modified so the corrector lands
   on unique branch points and the tangent's t_f-component becomes recoverable?

4. **The cold-seed basin (independent of continuation).** For a 40-rev indirect
   spiral, is a KKT-dual→costate map expected to land outside the shooting basin,
   and what's the more reliable warm-start? (Continuation in number of revs? in
   thrust magnitude / a_T? A single-arc adjoint sweep to polish the dual costates
   before multiple shooting? Something else?) Which of these is worth building
   first?

5. **Sanity check on abandoning ε entirely.** Was hard-throttling (no smoothing)
   the right call, or does the min-fuel indirect problem over 40 revs essentially
   *require* an ε-continuation (hyperbolic-tangent throttle smoothing) to have any
   basin at all — accepting the crawl as the price of admission?

If you think the whole IFS framing is the wrong tool for this problem, say that
too, and name what you'd do instead.
