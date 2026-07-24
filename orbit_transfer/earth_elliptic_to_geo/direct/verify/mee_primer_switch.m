function [primerVec, S, info] = mee_primer_switch(X, U, lam, sigma, dL, par)
% MEE_PRIMER_SWITCH  MEE/sigma-domain PMP primer vector + switching function.
%
% DERIVATION (eps=0, fuel objective; this is the only case exercised -- the
% homotopy's eps>0 running cost is quadratic in thr and the bang-bang
% monotonicity argument below does not apply there).
%
% casadi_lt_mee.m collocates, in sigma in [0,1], the state X=[P;ex;ey;hx;hy;
% m;t] with dX/dsigma = DeltaL * dXdL(X,U), dXdL = dXdt/Ldot(X,U) (lt_mee_rhs).
% DeltaL is a single free SCALAR (the total true-longitude span), not a
% per-node slack row -- so, unlike casadi_lt_2body's cScale, it cannot couple
% asymmetrically into any one node's stationarity condition; it is a common
% factor across every node and every state row (see below).
%
% Every row of dX/dsigma carries the SAME factor DeltaL/Ldot(X,U) (because
% lt_mee_rhs computes dXdL = dXdt/Ldot elementwise, uniformly across all 7
% rows), and casadi_lt_mee's own running-cost integrand is w = (DeltaL/Ldot)*
% (thr - eps*thr*(1-thr)) -- i.e. it carries that SAME factor too. So the
% sigma-domain Hamiltonian factors cleanly:
%
%   H_sigma(X,U,lam) = lam' * (dX/dsigma) + w
%                     = (DeltaL/Ldot(X,U)) * G(X,U,lam),   with
%   G(X,U,lam) = lam(1:5)'*dXdt(1:5) + lam(6)*mdot + lam(7)*1 + (thr - eps*thr*(1-thr))
%
% G is exactly the classical time-domain Hamiltonian bracket (Haberkorn-
% Martinon-Gergaud 2004 eq. for H, their p/pm playing the role of our lam)
% -- EXCEPT here Ldot(X,U) itself depends on control, through the RTN-normal
% thrust component w=thr*beta_3 (lt_mee_rhs's Ldot = Ldot0(X) + Tmax/m*sqrt(P/
% mu)*(hterm/Z)*w), because L (not t) is this transcription's independent
% variable. That is the one genuine structural difference from a plain
% free-horizon reparametrization: minimizing H_sigma over U is NOT simply
% minimizing G, because the positive scalar DeltaL/Ldot(X,U) also depends on
% U. We do NOT drop this coupling -- both stationarity conditions below are
% solved keeping it, and its size is reported (info.KLoverLdot0) for honesty.
%
% Write dXdt(1:5) = (Tmax/m)*B(X)*[q;s;w], q=thr*beta_R, s=thr*beta_T,
% w=thr*beta_N (beta unit RTN thrust direction; B is the Gauss matrix
% implicit in lt_mee_rhs, extracted numerically below, NOT hand-transcribed,
% by exploiting exact linearity in (q,s,w): B(:,j) = (m/Tmax)*dXdt(1:5)
% evaluated at beta=e_j, thr=1). Define the "MEE primer" p_el = B(X)'*
% lam(1:5) (3-vector, RTN components) and K_L(X) = Ldot(X,beta=e3,thr=1) -
% Ldot(X,thr=0) (the w-coefficient in Ldot, extracted the same way, not
% transcribed).
%
% (1) BETA stationarity (Lagrange multiplier mu for ||beta||^2=1, at fixed
%     thr): minimizing H_sigma over the unit sphere,
%       grad_beta H_sigma = DeltaL*[ (Tmax/m)*Ldot*p_el - G*K_L*e3 ] / Ldot^2
%                          = 2*mu*beta
%     so the STATIONARY thrust direction satisfies
%       beta* is parallel to  (Tmax/m)*Ldot*p_el - G*K_L*e3  =: primerVec
%     The MINIMIZING branch (H_sigma has a "+thr*(Tmax/m)*(p_el.beta)" term,
%     so H_sigma is reduced by pointing beta OPPOSITE p_el, exactly mirroring
%     the paper's u = -B'p/||B'p|| and this repo's own established primer
%     convention "thrust || -lam_v", verify_pmp_2body.m) is
%       beta* = -primerVec / ||primerVec||.
%     A single GLOBAL sign flip (IPOPT/CasADi Lagrangian-sign convention on
%     opti.dual, f +/- lam'*g -- ambiguous, standard, and NOT specific to any
%     anomaly) is still resolved empirically by the caller (verify_pmp_mee.m),
%     exactly as verify_pmp_2body.m already does; this function reports the
%     RAW (unflipped) primerVec and lets the caller apply one global flip to
%     lam (hence to primerVec and S together, consistently, since both are
%     linear in lam).
%
% (2) THR stationarity (fixed beta): H_sigma(thr) = DeltaL*(thr*C1+lam(7)) /
%     (Ldot0+thr*K_L*beta_3), C1 = (Tmax/m)*(p_el.beta) - (Tmax/c)*lam(6) + 1
%     (linear-fractional in thr since G and Ldot are both affine in thr at
%     fixed beta). Its thr-derivative reduces EXACTLY (the thr-dependent
%     terms cancel in the quotient-rule numerator) to
%       dH_sigma/dthr  proportional-to  S := C1*Ldot0 - lam(7)*K_L*beta_3
%     with S<0 => thr*=1 (H_sigma decreasing in thr, push to the upper
%     bound), S>0 => thr*=0 -- the MEE analog of the paper's switching
%     function psi (S<0 <=> full thrust), up to the positive overall factor
%     Ldot0*DeltaL/Ldot^2 that does not affect its sign.
%
% Both (1) and (2) are evaluated AT the solver's own (X,U) -- a necessary-
% condition CHECK on an already-converged trajectory, not a re-optimization
% (mirrors how verify_pmp_2body.m evaluates its S at the solver's own state).
%
% LUNAR-AWARE AMENDMENT (task B, 2026-07-23). The B(X)/pel extraction above
% ("exploiting exact linearity in (q,s,w)") implicitly assumed dXdt(1:5) is
% PURELY control-affine, i.e. dXdt(1:5)|_{thr=0} == 0 -- true for the plain
% 2-body RHS (test_mee_rhs.m's ballistic invariance: P,ex,ey,hx,hy freeze
% under zero thrust) but FALSE once lt_mee_rhs's par.pert (lunar third body)
% is active: the lunar direct+indirect acceleration is a control-INDEPENDENT
% forcing term, so dXdt(1:5) = A0(X,t) + (Tmax/m)*B(X)*[q;s;w] with A0(X,t)
% (the "ballistic bracket", zero unless par.pert.gain>0) the zero-throttle
% element-rate probed the SAME way Ldot0 already is (lt_mee_rhs(Xk,[e3;0],
% parK), never re-derived). Both consumers of the unit-thrust probes are
% amended to subtract this baseline before use, keeping the two derivations
% below EXACT rather than approximate:
%   (1') B(X) extraction: Bcol_j = (m/Tmax)*(dXdt(1:5)|_{beta=ej,thr=1} - A0),
%        not the raw probe -- A0 is CONTROL-INDEPENDENT so it cancels out of
%        every column identically and must not leak into p_el = B(X)'*lam_el.
%        (K_L = Le3-Ld0 needs NO such fix: the pert acceleration term is
%        control-independent there too, and a DIFFERENCE of two probes at
%        beta=e3,thr=1 vs thr=0 already cancels it exactly, whether or not
%        A0 is subtracted -- verified algebraically and by test 3/4 below.)
%   (2') THR stationarity: at fixed beta, G(X,U,lam) = G0 + thr*C1 where
%        G0 = lam(1:5)'*A0 + lam(7) (NOT just lam(7) -- the ballistic
%        bracket's projection onto the costate is a real, generally nonzero,
%        thr=0 offset once A0!=0). The quotient-rule cancellation that
%        isolates S is otherwise unchanged: S = C1*Ldot0 - G0*K_L*beta_3.
% Both amendments reduce IDENTICALLY to the pre-existing formulas when A0=0
% (pert absent or gain=0) -- G0 -> lam(7), Bcol_j unchanged -- so this file
% is BYTE-IDENTICAL on every already-certified 2-body artifact (A0 is exact
% floating-point zero there, not merely small).
%
% INPUTS:
%   X     - state trajectory [P;ex;ey;hx;hy;m;t] [7x(N+1)]
%   U     - control [beta(3);thr] [4x(N+1)]
%   lam   - nodal costate (mee_dual_to_costate output) [7x(N+1)]
%   sigma - node grid [(N+1)x1]
%   dL    - converged DeltaL (total true-longitude span) [scalar]
%   par   - kepler_lt_params struct (.Tmax .c .mu used; par.L is overwritten
%           per node internally). Optional par.pert (lunar_params struct,
%           CR3BP campaign): opt-in, forwarded verbatim to every lt_mee_rhs
%           probe below (parK = par so parK.pert carries automatically) --
%           absent/empty/gain==0 takes the untouched nominal path (see
%           amendment above and core/lt_mee_rhs.m's own pertOn branch).
%
% OUTPUTS:
%   primerVec - raw (unflipped) primer vector, RTN components [3x(N+1)]
%   S         - raw (unflipped) switching function [1x(N+1)]
%   info      - struct: .pel [3x(N+1)] (B(X)'*lam_el, unflipped, baseline-
%               subtracted), .Ldot0 [1x(N+1)], .Ldot [1x(N+1)] (actual, at
%               the solver's U), .KL [1x(N+1)], .G [1x(N+1)], .G0 [1x(N+1)]
%               (G at thr=0, the S-formula amendment above), .C1 [1x(N+1)],
%               .KLoverLdot0 [1x(N+1)] (coupling-strength diagnostic, should
%               be << 1 for "low thrust" to justify treating Ldot's
%               control-dependence as a small correction; NOT assumed small
%               anywhere in the code above, reported for honesty only),
%               .A0 [5x(N+1)] (the ballistic bracket, zero unless par.pert
%               is active), .A0overLdot0 [1x(N+1)] (|A0|/Ldot0 coupling-
%               strength diagnostic for the lunar term, same honesty-only
%               reporting convention as KLoverLdot0)
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, p.7 (H, psi, B'p, the
%       t-domain analog this L/sigma-domain derivation specializes to when
%       K_L -> 0).
%   [2] earth_elliptic_to_geo/lt_mee_rhs.m (dXdL, Ldot -- the RHS this file
%       probes, never hand-transcribed; par.pert branch is the lunar-aware
%       amendment's source of A0).
%   [3] earth_elliptic_to_geo/mee_dual_to_costate.m (the lam this consumes).
%   [4] earth_elliptic_to_geo/process/DESIGN_dual_map.md (Campaign-B context; no
%       cScale-class anomaly here -- see mee_dual_to_costate.m header).
%   [5] earth_elliptic_to_geo_CR3BP/direct/solve_cr3bp_minfuel.m (recorded the
%       CAVEAT this amendment closes: "the zero-throttle ballistic dXdt would
%       have to be subtracted out of the costate/primer bracket").
Nn = size(X, 2);
e1 = [1;0;0];  e2 = [0;1;0];  e3 = [0;0;1];

Ldot0 = zeros(1, Nn);  Ldot = zeros(1, Nn);  KL = zeros(1, Nn);
pel   = zeros(3, Nn);  G    = zeros(1, Nn);  C1 = zeros(1, Nn);
primerVec = zeros(3, Nn);  S = zeros(1, Nn);
A0all = zeros(5, Nn);  G0 = zeros(1, Nn);   % lunar-aware amendment (see header)

for k = 1:Nn
    parK = par;  parK.L = pi + sigma(k) * dL;
    Xk = X(:, k);  Uk = U(:, k);  mK = Xk(6);

    [dXdL0,   Ld0] = lt_mee_rhs(Xk, [e3; 0], parK);   % thr=0 probe -> Ldot0(X)
                                                       % AND the ballistic
                                                       % element baseline A0
    [dXdLe1,  Le1] = lt_mee_rhs(Xk, [e1; 1], parK);   % beta=e1(radial),  thr=1
    [dXdLe2,  Le2] = lt_mee_rhs(Xk, [e2; 1], parK);   % beta=e2(transv.), thr=1
    [dXdLe3,  Le3] = lt_mee_rhs(Xk, [e3; 1], parK);   % beta=e3(normal),  thr=1
    [dXdLk,   Ldk] = lt_mee_rhs(Xk, Uk,      parK);   % actual solver control

    Ldot0(k) = Ld0;
    KL(k)    = Le3 - Ld0;                              % w-coefficient in Ldot
    Ldot(k)  = Ldk;

    % Ballistic (control-independent) element-rate baseline: A0 = dXdt(1:5)
    % at thr=0. EXACTLY zero in the pert-absent 2-body case (test_mee_rhs.m's
    % ballistic-invariance check) -- nonzero only under par.pert (lunar
    % drift). Must be subtracted from each unit-thrust probe before the
    % difference isolates the control-affine B(X) columns (header amendment).
    A0 = dXdL0(1:5) * Ld0;
    A0all(:, k) = A0;

    Bcol1 = (dXdLe1(1:5) * Le1 - A0) * (mK / par.Tmax);
    Bcol2 = (dXdLe2(1:5) * Le2 - A0) * (mK / par.Tmax);
    Bcol3 = (dXdLe3(1:5) * Le3 - A0) * (mK / par.Tmax);
    lamEl = lam(1:5, k);
    pel(:, k) = [Bcol1.' * lamEl; Bcol2.' * lamEl; Bcol3.' * lamEl];

    dXdtk = dXdLk * Ldk;                               % = dXdt(X,Uk) exactly
    G(k)  = lam(:, k).' * dXdtk + Uk(4);                % eps=0: ell(thr)=thr
    G0(k) = lamEl.' * A0 + lam(7, k);                   % G at thr=0 (header amendment)

    primerVec(:, k) = (par.Tmax / mK) * Ldk * pel(:, k) - G(k) * KL(k) * e3;

    C1(k) = (par.Tmax / mK) * (pel(:, k).' * Uk(1:3)) - (par.Tmax / par.c) * lam(6, k) + 1;
    S(k)  = C1(k) * Ld0 - G0(k) * KL(k) * Uk(3);
end

info = struct('pel', pel, 'Ldot0', Ldot0, 'Ldot', Ldot, 'KL', KL, 'G', G, ...
              'G0', G0, 'C1', C1, 'KLoverLdot0', KL ./ Ldot0, ...
              'A0', A0all, 'A0overLdot0', sqrt(sum(A0all.^2, 1)) ./ Ldot0);
end
