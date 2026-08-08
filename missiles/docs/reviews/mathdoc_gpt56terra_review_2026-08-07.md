- **SEVERITY: MAJOR**  
  **LOCATION: Section 4.1, Eqs. (12)–(14)**  
  **WHY IT IS WRONG:** The derivation calls \(\boldsymbol{\omega}_f\) the horizon basis angular velocity “relative to inertial space,” but Eq. (14) includes only motion over the spherical surface. With rotation enabled it must include Earth rotation if inertial; the later use of rotating-frame pseudo-forces instead requires a basis rate relative to the rotating frame.  
  **CORRECTION:** Define \(\boldsymbol{\omega}_f\) as the horizon-basis angular velocity **relative to the planet-fixed rotating frame**, and retain Eq. (14). Alternatively use its inertial rate \(\boldsymbol{\Omega}+\boldsymbol{\omega}_f\) and do not add rotating-frame pseudo-forces.

- **SEVERITY: MAJOR**  
  **LOCATION: Section 8.3, Eqs. (Sänger range development)**  
  **WHY IT IS WRONG:** The closed-form range derivation silently treats \(r\), \(g\), and \(V_c\) as constants, although it is presented after an altitude-varying equilibrium condition. Without a constant-radius approximation, \(V_c^2/g=r\) varies along the path and Eq. (Sänger) does not follow. Also, \(V_c^2=gr\) is dimensionally/symbolically ambiguous: \(g\) was not defined.  
  **CORRECTION:** State explicitly before the energy equation: “Assume a shallow, near-constant-altitude glide, so \(r\simeq r_E\), \(g\simeq \mu/r_E^2\), and \(V_c^2\simeq g r_E\) are constants.” Write \(V_c^2=g r_E\), then derive \(R=(L/D)r_E\ln[1/(1-(V/V_c)^2)]/2\).

- **SEVERITY: MAJOR**  
  **LOCATION: Section 8.3, Physical reading after Eq. (Sänger)**  
  **WHY IT IS WRONG:** “The last increment of speed buys very little range” is the opposite of the formula. As \(V\to V_c^-\), \(\mathrm{d}R/\mathrm{d}V\) diverges; each increment closer to circular speed buys *more* range in this idealized model.  
  **CORRECTION:** Replace it with: “Conversely, the final increments toward \(V_c\) buy increasingly large range in this ideal model; this divergence signals the breakdown of the constant-altitude, fixed-\(L/D\), atmospheric-glide assumptions near orbital speed.”

- **SEVERITY: MAJOR**  
  **LOCATION: Section 9.3, Eq. (cross-track)**  
  **WHY IT IS WRONG:** The cross-track formula requires the initial bearing from launch to the **impact point**, not the vehicle heading \(\psi\) from the terminal flown state. A heading at impact is tangent at impact and generally differs from the initial great-circle bearing; using it invalidates the stated formula.  
  **CORRECTION:** Define \(\theta_{13}=\operatorname{bearing}(\text{launch},\text{impact})\) and \(\theta_{12}=\operatorname{bearing}(\text{launch},\text{target})=\psi_0\), then use  
  \[
  x_{\rm track}=r_I\arcsin\!\left[\sin(\Delta\theta_{13})\sin(\theta_{13}-\theta_{12})\right],
  \]
  with the chosen sign convention documented.

- **SEVERITY: MINOR**  
  **LOCATION: Section 5.4, “Numbers for the shipped placeholder booster”**  
  **WHY IT IS WRONG:** The stated “liftoff \(T/W=2.9899\)” uses vacuum thrust, despite immediately calculating delivered sea-level thrust as \(840112.6\) N. At liftoff under this model, delivered thrust must be used.  
  **CORRECTION:** Report vacuum \(T/W=2.9899\) separately, and state liftoff delivered \(T/W=840112.6/(32400\,g_0)=2.6441\).

- **SEVERITY: MINOR**  
  **LOCATION: Section 8.7, “Three degenerate cases”**  
  **WHY IT IS WRONG:** “At the poles every direction is south” is false at the south pole: every horizontal direction from the south pole is north.  
  **CORRECTION:** Write: “At the north pole every departure direction is south, and at the south pole every departure direction is north; longitude and a navigational bearing remain coordinate-dependent there.”
