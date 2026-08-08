### 1

**SEVERITY:** CRITICAL  
**LOCATION:** `run_ballistic_target.m`, `minEnergySolve`: residual `gamBoR - gamStarR`, with `gamStarR = pi/4 - angReqR/4`  
**WHY IT IS WRONG:** The classical condition uses the **free-flight central angle from burnout to impact**, with both endpoints on the same reference radius. The code substitutes pad-to-target angle even though burnout occurs downrange and about 80 km above the impact sphere. Matching this angle therefore does not minimize the energy of the actual coast arc. The reported agreement with `gammaStar` verifies an inapplicable residual, not minimum energy.  
**CORRECTED EXPRESSION:**
\[
\min_{\theta_{\rm loft},\,f_c}\;
\epsilon_{BO}=\frac{V_{BO}^{2}}{2}-\frac{\mu}{r_{BO}}
\]
subject to
\[
R(\theta_{\rm loft},f_c)-R_{\rm req}=0.
\]
If a classical diagnostic is retained, compute the actual coast angle
\[
\Lambda_{\rm coast}
=\operatorname{greatCircle}
(\phi_{BO},\lambda_{BO},\phi_I,\lambda_I)
\]
and compare against a formulation derived for \(r_{BO}\ne r_I\); do not use
\[
\gamma^*=\frac{\pi}{4}-\frac{\Lambda_{\rm pad-target}}{4}.
\]

### 2

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, mode name and description `'minimum-energy'`  
**WHY IT IS WRONG:** For a finite atmospheric boost, satisfying the vacuum equal-radius condition on burnout gamma is not equivalent to minimum booster energy, minimum propellant, minimum burnout mechanical energy, or minimum launch mass. With a fixed full-burn booster there is no energy optimization at all; with cutoff available there is, but its objective must be stated and minimized explicitly. The present mode is defensible only as “classical-gamma-matched,” not as minimum energy.  
**CORRECTED EXPRESSION:**
\[
\text{minimum-propellant:}\qquad
\min_{\theta_{\rm loft},\,f_c} f_c
\quad\text{subject to}\quad
R(\theta_{\rm loft},f_c)=R_{\rm req},
\]
or
\[
\text{minimum-burnout-energy:}\qquad
\min_{\theta_{\rm loft},\,f_c}
\left(\frac{V_{BO}^{2}}{2}-\frac{\mu}{r_{BO}}\right)
\quad\text{subject to}\quad
R(\theta_{\rm loft},f_c)=R_{\rm req}.
\]
If neither objective is solved, rename the mode to `classical-gamma-matched`.

### 3

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `maxRangeLoft` returning `loftStarR = 0.5.*(aL + bL)`, followed by branch brackets `[loftLoR,loftStarR]` and `[loftStarR,loftHiR]`  
**WHY IT IS WRONG:** Golden section establishes only that the true maximizer remains inside `[aL,bL]`; its midpoint is not guaranteed to lie at the maximizer. Unless the midpoint happens to be exact, one branch bracket crosses the true maximum and is therefore not monotone. A bisection can then select an unintended root or reject a reachable target. The claimed one-sided bracketing is not certified.  
**CORRECTED EXPRESSION:**
\[
\theta_{\max}\in[a_L,b_L].
\]
Use certified one-sided intervals
\[
[\theta_{\min},a_L]
\quad\text{and}\quad
[b_L,\theta_{\max,\mathrm{user}}],
\]
and classify targets lying in the unresolved top band as coalesced:
\[
R_{\rm req}\ge
\min\!\left(R(a_L),R(b_L)\right)-\mathrm{tol}_R.
\]
Alternatively refine `[aL,bL]` until its induced range uncertainty is below the branch-separation tolerance.

### 4

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `maxRangeLoft`: `nTurn` is only printed as a caution after both branch solves  
**WHY IT IS WRONG:** One coarse sign change does not establish global unimodality; multiple sign changes explicitly disprove the assumption. Nevertheless the code continues with a local golden-section maximum and performs bisections that require monotonic branches. It can report a converged but wrong branch before merely printing a caution. Flat samples and unresolved extrema between 10.4-degree scan intervals make even `nTurn == 1` weak evidence.  
**CORRECTED EXPRESSION:**
```matlab
if nTurn ~= 1 || any(diff(scanRngM(1:kB)) <= 0) || ...
                 any(diff(scanRngM(kB:end)) >= 0)
    error('coorbital:runBallisticTarget:rangeNotUnimodal', ...
          'The sampled range curve does not certify two monotone branches.');
end
```
A production implementation should adaptively refine each interval until the monotonicity evidence is resolved.

### 5

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `measureBranch(hApoM,tFlyS,hApoStar,tFlyStar)`  
**WHY IT IS WRONG:** Apogee and total flight time are not mathematical branch invariants for a finite powered, atmospheric trajectory. Drag, lift, burnout altitude and boost duration can make either nonmonotone with commanded loft. Near maximum, both differences vanish quadratically and their signs are dominated by maximum-search and integration errors. Agreement between two correlated observables does not prove the root lies on the requested side.  
**CORRECTED EXPRESSION:**
```matlab
if pick.loftR < aL - tolClassR
    flownName = 'depressed';
elseif pick.loftR > bL + tolClassR
    flownName = 'lofted';
else
    flownName = 'coalesced/indeterminate';
end
```
Branch identity must come from the root’s location relative to the certified maximizer interval, with apogee and time retained only as descriptive outputs.

### 6

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `minEnergySolve`: inner bracket `[cutFracMin,1]` and claim that range rises monotonically with cutoff  
**WHY IT IS WRONG:** More burn does not guarantee more impact range at fixed loft. Added speed can move the trajectory through its own range maximum; added burn time also changes burnout position, altitude, gamma, gravity loss and drag loss. The outer full-burn loft bracket proves only that the endpoint at `cutFrac = 1` reaches or exceeds the target—it proves neither that `cutFracMin` undershoots nor that the interval contains a unique monotone root. Near the full-burn maximum the Jacobian becomes poorly conditioned and the inner root may disappear or bifurcate.  
**CORRECTED EXPRESSION:**
\[
F_1(\theta,f)=R(\theta,f)-R_{\rm req}=0
\]
must be bracketed using sampled sign changes in \(f\):
\[
F_1(\theta,f_i)\,F_1(\theta,f_{i+1})\le0,
\]
with monotonicity verified on the selected interval. If multiple cutoff roots exist, select the one minimizing the declared objective; do not assume `[cutFracMin,1]` contains one unique root.

### 7

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `minEnergySolve`: outer bisection on `gammaAtRange(loftR)`  
**WHY IT IS WRONG:** The outer residual is not demonstrated to be monotone; nine samples from one placeholder case are not a property of the model family. Worse, each inner solve is only accurate to `tolRngM`, so `gammaAtRange` is a noisy, potentially discontinuous map when the selected cutoff root changes. Bisection can terminate on inner-solve noise rather than on the requested gamma tolerance, particularly where \(\partial R/\partial f\) or \(\partial R/\partial\theta\) approaches zero.  
**CORRECTED EXPRESSION:**
\[
\mathbf F(\theta,f)=
\begin{bmatrix}
R(\theta,f)-R_{\rm req}\\[2mm]
\partial \epsilon_{BO}(\theta,f)/\partial s
\end{bmatrix}
=\mathbf0
\]
for a true constrained energy optimum, solved with a bounded trust-region method and a finite-difference Jacobian. At minimum, require
\[
\left|\frac{\partial R}{\partial f}\right|>\eta_R,\qquad
\operatorname{sign}\Delta\gamma_{BO}\ \text{constant over the outer bracket},
\]
and make inner error small enough that its propagated gamma uncertainty is below `tolGamR`.

### 8

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `atMaxRange = pick.loftR == loftStarR`  
**WHY IT IS WRONG:** This is the reported exact-equality defect. Equality depends on whether `rangeSolve` returns the exact endpoint bit pattern. A trajectory numerically indistinguishable from the maximum can be classified as an ordinary branch, while a copied endpoint can be classified as coalesced despite appreciable uncertainty in the actual maximizer.  
**CORRECTED EXPRESSION:**
```matlab
atMaxRange = pick.rngM >= rngMax - tolRngM && ...
             pick.loftR >= aL - tolLoftR && ...
             pick.loftR <= bL + tolLoftR;
```
Here `[aL,bL]` must be the retained final maximizer bracket, not only its midpoint.

### 9

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, `flyLoft`: completion test `nPhRn < 3 || abs(hEndM - cfg.hStop) > 1e-3`  
**WHY IT IS WRONG:** Presence of three phase labels plus a final altitude does not prove that boost ended by exhaustion/cutoff or that phase 2 ended at apogee. A horizon-truncated coast can still start phase 3 and later hit the altitude event, passing this test and feeding a physically malformed trajectory into every root solve. The later `why*` diagnostics occur only after selection and cannot protect the optimizer.  
**CORRECTED EXPRESSION:**
```matlab
kBO  = find(traj.phaseIdx == 1,1,'last');
kApo = find(traj.phaseIdx == 2,1,'last');

boostOK = (cutFrac < 1 && abs(traj.t(kBO)-cutFrac*cfg.tBurn) <= tolTime) || ...
          (cutFrac == 1 && abs(traj.x(kBO,7)-cfg.mBurnout) <= tolMass);
coastOK = abs(traj.x(kApo,5)) <= tolGamma && ...
          traj.t(kApo)-traj.t(kBO) < cfg.tMaxCoast-tolTime;
impactOK = abs(hEndM-cfg.hStop) <= tolAltitude && ...
           traj.x(end,5) < 0;

if ~(boostOK && coastOK && impactOK)
    error(...);
end
```

### 10

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic.m`, drag attribution: `dragDiff = sDrgB - sVac`, while `envVac.grav = @sphereGrav` but `envDrg.grav = env.grav`  
**WHY IT IS WRONG:** With an overridden gravity model, the alleged drag effect includes the gravity-model difference. The summary nevertheless states that the sign is a drag-only physical check. The same contamination affects `liftDiff` if rotation or other environment settings differ.  
**CORRECTED EXPRESSION:**
```matlab
envDrg = envVac;
envDrg.atmos = env.atmos;
envDrg.aero = @(alphaArg,machArg,vehArg) ...
    dragOnly(aeroFn,alphaArg,machArg,vehArg);

envFullCheck = envDrg;
envFullCheck.aero = aeroFn;

dragDiff = sDrgB - sVacProp;
liftDiff = sFullCheckB - sDrgB;
```
All three comparison trajectories must use identical gravity and rotation settings.

### 11

**SEVERITY:** MINOR  
**LOCATION:** `run_ballistic_target.m`, `maxRangeLoft`: loop condition `nGolden < maxStep` followed by unconditional return  
**WHY IT IS WRONG:** If an extremely tight valid positive `tolLoftR` reaches 200 iterations, the function silently returns an unconverged maximizer while reporting its residual width. This can invalidate both branch brackets without a failure.  
**CORRECTED EXPRESSION:**
```matlab
if (bL - aL) > tolLoftR
    error('coorbital:runBallisticTarget:maxRangeNoConvergence', ...
          'Golden-section search reached %d steps with width %.16g rad, tolerance %.16g rad.', ...
          maxStep,bL-aL,tolLoftR);
end
```

### 12

**SEVERITY:** MAJOR  
**LOCATION:** `run_ballistic_target.m`, shipped `alphaMax = 12` versus `run_ballistic.m` value `6`  
**WHY IT IS WRONG:** A control/structural/aerothermal limit cannot be chosen to make a desired branch reachable. Changing it changes the vehicle definition, boost loads, losses and range envelope; it is not a targeting degree of freedom unless the vehicle is certified for that authority. Shipping two entry points for nominally the same vehicle with different limits makes their performance non-comparable and invites the appearance that the envelope was tuned to the demonstration.  
**CORRECTED EXPRESSION:**
```matlab
alphaMaxR = min(alphaCommandLimitR, ...
                vehicleAllowableAlpha(mach,qbar,loadFactor,veh));
```
Use one vehicle-qualified limit in both scripts. If 6 degrees is the qualified value, report the depressed demonstration target as unreachable; if 12 degrees is qualified, use 12 degrees consistently and document the qualification basis rather than the range benefit.
