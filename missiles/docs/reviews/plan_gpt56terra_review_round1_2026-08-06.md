## Findings

- **CRITICAL — latitudinal-gravity sign:** `+ gLat.*sin(gamma).*cos(psi)./V` is wrong for `gLat` defined northward.  
  **Correct:** `- gLat.*sin(gamma).*cos(psi)./V`. The northward gravity projection onto the positive flight-path-normal direction is negative.

- **MAJOR — equilibrium-glide balance:** `0.5 rho V^2 S CL = m g - m V^2 / r` omits the flight-path-angle factor.  
  **Correct:** `0.5*rho*V^2*S*CL = m*(g - V^2/r)*cos(gamma)`, hence `V^2 = m*g*cos(gamma)/(0.5*rho*S*CL + m*cos(gamma)/r)`. At \(-0.5^\circ\) this omission is tiny, but a fixed-\(\alpha\), open-loop trajectory is not guaranteed to settle onto this instantaneous \(\dot\gamma=0\) curve; the stated 5%/20% criterion is not a reliable validation reference.

- **MAJOR — Allen–Eggers measured quantity:** `accel = abs(diff(vArc)./diff(tArc));` is net longitudinal deceleration, not the aerodynamic deceleration in the Allen–Eggers result.  
  **Correct:** `aDrag = -Vdot - gr.*sin(gamma)` (or directly `0.5*rho.*V.^2*S*CD/m`). The quoted \(a_{\max}=V_e^2|\sin\gamma_e|/(2eH)\) and ballistic-coefficient independence are correct only under its flat, constant-\(\gamma\), gravity-free assumptions. At \(30^\circ\), 10% is defensible, but only after comparing drag deceleration.

- **MAJOR — polar singularities:** `londot = .../(r.*cos(lat))` and `...sin(psi).*tan(lat)./r` are unguarded at `cos(lat)=0`. No finite corrected expression exists in longitude/latitude coordinates at a pole; require `abs(cos(lat)) > tolerance` or switch to Cartesian/local coordinates.

- **MAJOR — vertical/zero-speed singularities:** `aLift.*sin(sigma)./(V.*cos(gamma))`, `...gLat.../(V.*cos(gamma))`, and the rotating term `cos(psi).*tan(gamma)` are singular as `cos(gamma)->0`; all dynamics with `/V` are singular as `V->0`. These are coordinate/model-domain singularities, so enforce `V > Vmin` and `abs(cos(gamma)) > tolerance` or use a nonsingular velocity-vector representation.

- **MINOR — isothermal atmosphere consistency:** `rho = rho0*exp(-h/H)` with fixed `T=250 K` is not exactly hydrostatic under `gr=mu/r^2`. The hydrostatic isothermal spherical expression is `rho(r)=rho0*exp((mu/(R*T))*(1/r - 1/rE))`. The 120-km model has unrealistic thermospheric density/temperature and therefore sound speed/Mach, but it does **not** invalidate either analytic check because both references use the imposed exponential density law.

The rotating-Earth Coriolis and centrifugal terms shown are otherwise self-consistent for planet-relative velocity and heading clockwise from north; no missing rotating-frame term was found.
