**SEVERITY: MAJOR**  
**QUOTE:** `T = max(bst.thrustVac - bst.Aexit*P, 0)` with `mdot = bst.thrustVac/(bst.Isp*c.g0)`  
**WHY IT IS WRONG:** Clamping a negative pressure-corrected thrust to zero while continuing to consume propellant violates the stated fixed-flow nozzle model. Fixed `mdot` is physically reasonable for a choked, fixed-chamber-pressure engine; the unilateral force clamp is not.  
**CORRECTED EXPRESSION:** Either retain the model consistently:
\[
T=T_{\rm vac}-A_eP,\qquad \dot m={T_{\rm vac}\over I_{\rm sp,vac}g_0},
\]
or, if zero thrust denotes commanded/shutdown behavior:
\[
T=\max(T_{\rm vac}-A_eP,0),\qquad
\dot m=\begin{cases}
T_{\rm vac}/(I_{\rm sp,vac}g_0),&T>0\\
0,&T=0.
\end{cases}
\]
`Isp` must explicitly mean vacuum \(I_{\rm sp,vac}\); then deriving fixed \(\dot m\) from vacuum thrust is self-consistent.

**SEVERITY: MAJOR**  
**QUOTE:** “With `T = 0` and `mdot = 0` … [reduction] is the cheapest and strongest check available”  
**WHY IT IS WRONG:** This reduction cannot validate any term that vanishes with thrust or mass flow: all thrust projections/signs, `dm/dt`, use of state mass versus `veh.mass` when \(x_7=veh.mass\), or thrust–lift coupling. It also cannot find an error shared by `boost3DOF` and `glide3DOF`.  
**CORRECTED EXPRESSION:** Retain the reduction test, but add an independent force-increment test at \(T>0\), \(\alpha\ne0\), \(\sigma\ne0\), with identical states and all non-thrust models unchanged:
\[
\Delta\dot V={T\cos\alpha\over m},\qquad
\Delta\dot\gamma={T\sin\alpha\cos\sigma\over mV},
\]
\[
\Delta\dot\psi={T\sin\alpha\sin\sigma\over mV\cos\gamma},\qquad
\dot m=-\dot m_{\rm prop}.
\]
Use \(m\ne veh.mass\). This directly checks all new thrust-component signs and denominators.

**SEVERITY: MAJOR**  
**QUOTE:** “Three checks, in increasing strength.”  
**WHY IT IS WRONG:** The listed Tsiolkovsky test uses `alpha = 0`; therefore it cannot detect replacing \(T\cos\alpha\) with \(T\), nor any wrong sign/projection in the normal or banked thrust components. The plan notices this only during an intentional mutation, rather than making it a required validation from the start.  
**CORRECTED EXPRESSION:** Make the nonzero-\(\alpha\), nonzero-\(\sigma\) force-increment test above a required Task 6 check, not a conditional addition. Tsiolkovsky is exact only under the additional assumptions \(D=0\), constant \(T\), constant \(\dot m\), \(T=I_{\rm sp}g_0\dot m\), and comparison of planet-relative speed:
\[
V_f-V_0=I_{\rm sp}g_0\ln\!\left({m_0\over m_f}\right).
\]

**SEVERITY: MAJOR**  
**QUOTE:** “A single-stage solid-equivalent booster sized to loft the existing 900 kg glide vehicle. Fields: `massWet`, `massDry`…”  
**WHY IT IS WRONG:** The plan never defines whether these masses include the 900 kg payload, retained dry booster hardware, and/or a stage-separation mass drop. That ambiguity corrupts the acceleration denominator, liftoff T/W, burnout event, rocket-equation mass ratio, and coast initial mass. A continuous seven-state chain cannot silently discard dry-stage mass.  
**CORRECTED EXPRESSION:** Define the state mass explicitly as total carried mass:
\[
m(0)=m_{\rm payload}+m_{\rm booster,dry}+m_{\rm propellant},
\qquad
m_{\rm burnout}=m_{\rm payload}+m_{\rm booster,dry},
\]
\[
t_{\rm burn}={m_{\rm propellant}\over\dot m}
={m(0)-m_{\rm burnout}\over\dot m}.
\]
If stage separation is intended, add an explicit instantaneous junction map
\[
m^+=m^- - m_{\rm jettison},
\]
and do not claim a state-continuous seven-state chain across that junction.

**SEVERITY: MAJOR**  
**QUOTE:** “For a vacuum ballistic trajectory over a spherical Earth with burnout speed `V`, flight-path angle `gamma`, at radius `r`, the range angle has a closed form. Compute it independently and compare against the propagated ground range.”  
**WHY IT IS WRONG:** No formula is specified, and comparing a closed-form *vacuum* coast with a propagated trajectory that includes drag and a terminal descent is not an analytic cross-check. It is only an uncalibrated model-discrepancy observation. The analytic reference must use the same burnout state and a separately propagated vacuum, central-gravity coast.  
**CORRECTED EXPRESSION:** For \(\mu\), impact radius \(R_I\), burnout \(r_0,V_0,\gamma_0\), define
\[
h=r_0V_0\cos\gamma_0,\qquad
p={h^2\over\mu},\qquad
\epsilon={V_0^2\over2}-{\mu\over r_0},
\]
\[
e=\sqrt{1+{2\epsilon h^2\over\mu^2}},\qquad
f_0=\operatorname{atan2}\!\left({hV_0\sin\gamma_0\over\mu},
{p\over r_0}-1\right),
\]
\[
q={p/R_I-1\over e},\qquad
f_I=2\pi-\arccos(q),
\]
choosing the forward, descending intersection and requiring \(|q|\le1\). Then
\[
\Delta\chi=\operatorname{mod}(f_I-f_0,\,2\pi),\qquad
s_{\rm vac}=R_I\Delta\chi.
\]
Validate this against a separately integrated vacuum trajectory with \(D=L=0\), no rotation, and \(\mathbf g=-\mu\mathbf r/r^3\); report atmospheric/full-trajectory discrepancy separately.
