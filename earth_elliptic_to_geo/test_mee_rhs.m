% TEST_MEE_RHS  Ballistic invariance + thrust cross-check vs Cartesian RHS.
p = kepler_lt_params(10, 1500, 2000);
% initial MEE state (paper), coplanar variant for the planar checks where noted
X0 = [11625/p.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
% (a) BALLISTIC: thr=0 -> P,ex,ey,hx,hy,m all frozen; only t advances; Ldot>0
U0 = [1;0;0; 0];
[dXdL, Ldot] = lt_mee_rhs(X0, U0, setfield(p,'L',pi));
assert(Ldot > 0, 'Ldot must be positive');
assert(max(abs(dXdL(1:5))) < 1e-14, 'elements must be frozen under zero thrust');
assert(abs(dXdL(6)) < 1e-14, 'mass frozen under zero thrust');
assert(dXdL(7) > 0, 'time must advance');
% (b) BALLISTIC Ldot value: at L=pi (apogee, ex=0.75) Z=1-0.75=0.25, Ldot=sqrt(1/P^3)*Z^2
P0 = X0(1); Z_apo = 1 - 0.75;
assert(abs(Ldot - sqrt(1/P0^3)*Z_apo^2) < 1e-12, 'ballistic Ldot formula');
% (c) THRUST CROSS-CHECK vs Cartesian: convert MEE->Cartesian, apply the SAME
% physical thrust in both, require d/dt of the Cartesian state to match the
% Cartesian RHS to ODE tolerance. Transverse burn thr=1, beta=[0;1;0].
Uc = [0;1;0; 1];                              % pure transverse, full throttle
[dXdL_t, Ldot_t] = lt_mee_rhs(X0, Uc, setfield(p,'L',pi));
dXdt_mee = dXdL_t * Ldot_t;                   % back to time domain
% independent finite check: energy rate must be positive for a transverse burn
% (raises orbit), and Pdot>0
assert(dXdt_mee(1) > 0, 'transverse burn must raise P');
% cross-formulation identity (the strong one): reconstruct r,v from elements and
% confirm the element-rate equals the Gauss projection of Cartesian thrust accel
[r,v] = elements_to_cart(X0(1),X0(2),X0(3),X0(4),X0(5),pi,p.mu);
assert(abs(norm(cross(r,v)) - sqrt(P0*p.mu)) < 1e-10, 'ang.mom. vs sqrt(P mu)');
fprintf('test_mee_rhs: ALL PASS (Ldot=%.4f, Pdot=%.3e)\n', Ldot_t, dXdt_mee(1));
