### 1
**SEVERITY:** CRITICAL  
**LOCATION:** `missiles/+coorbital/+guide/pitchProgram.m` — `alpha = theta - gamma` while `sigma` may be nonzero  
**WHY IT IS WRONG:** This relation holds only when the thrust/body-axis normal lies in the vertical plane (`sigma = 0`). With the EOM’s bank convention, the thrust direction is  
\[
\hat b=\cos\alpha\,\hat v+\sin\alpha(\cos\sigma\,\hat e_\gamma+\sin\sigma\,\hat e_\psi),
\]
so its radial component is  
\[
\sin\theta=\sin\gamma\cos\alpha+\cos\gamma\cos\sigma\sin\alpha.
\]
The current guidance does not command the documented pitch attitude when banked; the self-demo explicitly combines pitch commands with ±45° bank.  
**CORRECTED EXPRESSION:**  
\[
R=\sqrt{\sin^2\gamma+\cos^2\gamma\cos^2\sigma},\qquad
\phi=\operatorname{atan2}(\sin\gamma,\cos\gamma\cos\sigma),
\]
\[
\alpha=\asin\!\left(\frac{\sin\theta}{R}\right)-\phi,
\]
using the branch continuous near the current command, and rejecting `abs(sin(theta)) > R`. Alternatively, require and enforce `sigma == 0` for this guidance law.

### 2
**SEVERITY:** CRITICAL  
**LOCATION:** `missiles/+coorbital/+prop/phaseRun.m` — single `veh` argument used by every phase: `ph.eom(...,veh,env)`  
**WHY IT IS WRONG:** The documented staging repair is impossible through this interface. A link can change `x(7)`, but cannot replace booster geometry, aerodynamics and `veh.mass`; the following `massConstant` phase either aborts or, without that wrapper, flies the wrong vehicle. This defeats the stated mass contract at exactly the phase junction where it matters.  
**CORRECTED EXPRESSION:** Each phase must carry its complete vehicle:
```matlab
phaseVeh = ph.veh;
odeF = @(t,x) ph.eom(t - ph.tspan(1),x, ...
    ph.guide(t - ph.tspan(1),x),phaseVeh,env);
```
The post-separation phase must supply a payload vehicle whose `mass`, `Sref` and aerodynamic data correspond to the linked state.

### 3
**SEVERITY:** CRITICAL  
**LOCATION:** `missiles/+coorbital/+eom/massConstant.m` — `if abs(x(7) - mVeh) > tol`  
**WHY IT IS WRONG:** MATLAB comparisons with `NaN` are false. Therefore `x(7)=NaN` or `veh.mass=NaN` passes the guard and feeds NaN into the EOM. Equal infinities can also pass because `Inf-Inf` is NaN. Zero and negative agreeing masses pass and cause division by zero or inverted aerodynamic acceleration.  
**CORRECTED EXPRESSION:**
```matlab
mState = x(7);
mVeh = veh.mass;
if ~isreal(mState) || ~isscalar(mState) || ~isfinite(mState) || mState <= 0 || ...
   ~isreal(mVeh) || ~isscalar(mVeh) || ~isfinite(mVeh) || mVeh <= 0
    error('coorbital:massConstant:invalidMass', ...
        'State and vehicle masses must be finite positive real scalars.');
end
tol = max(1e-9,1e-9*max(abs(mState),abs(mVeh)));
if abs(mState-mVeh) > tol
    error('coorbital:massConstant:massMismatch',...);
end
```

### 4
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+prop/phaseRun.m` — EOM, guide and event receive raw `t` from `ph.tspan`  
**WHY IT IS WRONG:** Their contracts call this time “since phase start,” but a phase with `tspan=[10 50]` receives 10–50 rather than 0–40. Time-dependent thrust and guidance therefore depend on the arbitrary numerical origin of `tspan`. The output clock is shifted separately, creating inconsistent clocks.  
**CORRECTED EXPRESSION:**
```matlab
tau = @(t) t - ph.tspan(1);
odeF = @(t,x) ph.eom(tau(t),x,ph.guide(tau(t),x),phaseVeh,env);
eventF = @(t,x) ph.terminate(tau(t),x);
opts = odeset('RelTol',relTol,'AbsTol',absTol,'Events',eventF);
```
Control-history reconstruction must likewise call `ph.guide(tk(kt)-ph.tspan(1),...)`.

### 5
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+prop/phaseRun.m` — `[tk,xk] = ode45(...)` with no event-status check  
**WHY IT IS WRONG:** The code cannot distinguish event termination from exhaustion of `tspan`. A failed burnout or altitude event silently advances the chain and may execute staging with live propellant or begin descent above the intended interface.  
**CORRECTED EXPRESSION:**
```matlab
[tk,xk,te,xe,ie] = ode45(odeF,ph.tspan,xCurr,opts);
if ph.requireEvent && isempty(ie)
    error('coorbital:phaseRun:requiredEventMissing', ...
        'Phase %d reached tspan(end) without its required terminal event.',kp);
end
```
The phase schema must explicitly distinguish time-limited phases from phases requiring an event.

### 6
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+prop/phaseRun.m` — scalar `AbsTol = 1e-10` for all state components  
**WHY IT IS WRONG:** Radius, angles, speed and mass have incompatible scales and units. A scalar tolerance gives no meaningful physical error budget; `RelTol=1e-10` still permits roughly \(6\times10^{-4}\) m local radius error while demanding \(10^{-10}\) kg near zero mass. No `MaxStep` or phase-break handling is provided for short burns, guidance knots or narrow event regions.  
**CORRECTED EXPRESSION:** Use a state-scaled vector, for example:
```matlab
absTol = [1e-3; 1e-11; 1e-11; 1e-6; 1e-11; 1e-11; 1e-6];
opts = odeset('RelTol',1e-10,'AbsTol',absTol, ...
              'MaxStep',ph.maxStep,'Events',eventF);
```
The seventh entry must be omitted for six-state chains, and scheduled discontinuities must be phase boundaries or explicit integration breakpoints.

### 7
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+prop/phaseRun.m` — repeated `ph.guide(...)` calls during sizing and history reconstruction  
**WHY IT IS WRONG:** A stateful guide receives not one extra evaluation, but one sizing call plus one call per returned output sample after all solver-stage calls. Its reconstructed `traj.u` can therefore differ from the controls actually integrated and can alter later controller state. The documentation materially understates this.  
**CORRECTED EXPRESSION:** Require
\[
u(t,x)=\text{a pure, deterministic function of only }(t,x).
\]
Stateful guidance must carry controller state in the integrated state vector; `traj.u` must then be reconstructed from that state without mutating hidden state.

### 8
**SEVERITY:** MAJOR  
**LOCATION:** `glide3DOF.m` and `boost3DOF.m` — guards using comparisons such as `V < 1`, `m <= 0`, and `abs(cos(...)) < 1e-8`  
**WHY IT IS WRONG:** NaN values pass every guard. There is no validation of finite positive `r`, finite controls, finite atmosphere/gravity/aerodynamic outputs, positive sound speed, or finite derivatives. `r=0`, `aSnd=0`, invalid model output or NaN latitude can therefore produce silent NaNs. The arbitrary `V<1` rejection also makes launch from rest unsupported without resolving the underlying \(1/V\) coordinate singularity.  
**CORRECTED EXPRESSION:** Before divisions, require all inputs and model outputs to be finite real scalars, with
\[
r>0,\quad V>V_{\min},\quad m>0,\quad a_{\rm snd}>0,\quad \rho\ge0,\quad gr\ge0.
\]
After assembly:
```matlab
if any(~isfinite(xdot)) || ~isreal(xdot)
    error('coorbital:eom:nonFiniteDerivative', ...
        'The EOM produced a non-finite derivative.');
end
```
Launch through \(V=0\), polar passage and vertical flight require Cartesian velocity equations, not threshold guards.

### 9
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+prop/constThrust.m` — `mdot = ... .* (T > 0)`  
**WHY IT IS WRONG:** A choked motor does not stop burning merely because ambient back-pressure makes the simplified net-thrust expression nonpositive. Zeroing both outputs can leave mass frozen above burnout, after which `phaseRun` silently reaches `tspan(end)`. The condition indicates that the nozzle model is invalid, not that combustion ceased.  
**CORRECTED EXPRESSION:**
```matlab
Traw = veh.thrustVac - veh.Aexit.*P;
if any(Traw <= 0)
    error('coorbital:constThrust:invalidBackPressure', ...
        'Back pressure is outside the valid domain of the nozzle model.');
end
T = Traw;
mdot = veh.thrustVac./(veh.Isp.*c.g0);
```
A separated-flow model may replace the error, but mass flow must remain tied to commanded motor operation.

### 10
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+util/rangeSolve.m` — no validation of `fTarget`  
**WHY IT IS WRONG:** A NaN target makes all convergence and sign comparisons false and drives the bisection to an arbitrary failure result with NaN residual instead of a loud input error. Complex or nonscalar targets are similarly unsupported but unchecked.  
**CORRECTED EXPRESSION:**
```matlab
if ~isscalar(fTarget) || ~isreal(fTarget) || ~isfinite(fTarget)
    error('coorbital:rangeSolve:badTarget', ...
        'fTarget must be a finite real scalar; got %s.',mat2str(fTarget));
end
```

### 11
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+util/rangeSolve.m` — failure returns the final midpoint and claims it is the closest value reached  
**WHY IT IS WRONG:** The final midpoint need not minimize `abs(fRange(x)-fTarget)` for nonlinear, discontinuous or noisy functions. No best-so-far value is retained, contradicting the output contract. The failure message also reports the original endpoints rather than the live collapsed bracket.  
**CORRECTED EXPRESSION:** Initialize and update:
```matlab
[bestErr,kBest] = min([abs(fLo-fTarget),abs(fHi-fTarget)]);
bestX = [xLo,xHi]; bestX = bestX(kBest);
bestF = [fLo,fHi]; bestF = bestF(kBest);

if abs(fMid-fTarget) < bestErr
    bestErr = abs(fMid-fTarget);
    bestX = xMid;
    bestF = fMid;
end
```
On failure return `xSol=bestX; fAch=bestF`, and report `aLo`, `aHi` and their cached achieved values.

### 12
**SEVERITY:** MINOR  
**LOCATION:** `missiles/+coorbital/+util/rangeSolve.m` — `xMid = 0.5.*(aLo + aHi)` and collapse check after evaluation  
**WHY IT IS WRONG:** The sum can overflow for finite same-sign endpoints. Once the bracket reaches machine spacing, `xMid` may equal an endpoint, causing a duplicate expensive evaluation despite the stated guarantee. Two hundred iterations are not sufficient for every representable finite interval.  
**CORRECTED EXPRESSION:**
```matlab
xMid = aLo/2 + aHi/2;
if xMid == aLo || xMid == aHi
    stop = true;
    break;
end
```
Derive the iteration bound from the initial width and endpoint spacing, or continue until midpoint stagnation rather than using a fixed 200-step claim.

### 13
**SEVERITY:** MAJOR  
**LOCATION:** `missiles/+coorbital/+util/greatCircleBearing.m` — deliberate return of plausible bearings for coincident, antipodal and polar cases  
**WHY IT IS WRONG:** These are undefined targeting azimuths, not harmless conventions. Returning zero or a rounding-dependent direction can launch a valid-looking trajectory toward an arbitrary course. Near-antipodal bearings are also catastrophically ill-conditioned.  
**CORRECTED EXPRESSION:** Form
\[
y=\sin(\Delta\lambda)\cos\phi_2,\qquad
x=\cos\phi_1\sin\phi_2-\sin\phi_1\cos\phi_2\cos(\Delta\lambda),
\]
and reject when
\[
\operatorname{hypot}(x,y)\le \tau
\]
or `abs(cos(lat1)) <= tau`, using a documented angular tolerance. Only otherwise compute
\[
\psi_0=\operatorname{mod}(\operatorname{atan2}(y,x),2\pi).
\]
