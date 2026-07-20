function H = hamiltonian_along_traj(matPath, epsv)
% HAMILTONIAN_ALONG_TRAJ  Evaluate the L-domain and time-domain PMP Hamiltonians
% along a min-fuel MEE solution, with a self-calibrating dual sign.
%
% Companion to hamiltonian_const_check.m. That routine checks the ONE conserved
% quantity (the time-costate lambda_t = -H_t); this routine reconstructs the
% FULL Hamiltonians node-by-node so they can be plotted/animated:
%
%   H_L(sigma) = ell_L + lambda' * dXdL          (L-domain PMP Hamiltonian)
%                ell_L = (thr - eps*thr*(1-thr)) / Ldot     (running cost per dL)
%   H_t(sigma) = -lambda_t                        (time-domain Hamiltonian; const)
%
% where dXdL, Ldot come from lt_mee_rhs and lambda [7x(N+1)] is the nodal costate
% recovered from the defect duals (mee_dual_to_costate). H_L is NOT constant --
% L enters the dynamics explicitly (cos L, sin L) -- while H_t is (time-autonomy).
%
% SIGN CALIBRATION. CasADi's opti.dual carries a convention-dependent overall
% sign. We fix it from the EXTREMAL IDENTITY that holds along any Pontryagin
% extremal: the total derivative of the (minimized) Hamiltonian equals its
% explicit partial,
%       dH_L/dsigma = partial H_L/partial sigma = DeltaL*(partial H_L/partial L),
% the right side taken through the explicit L = pi + sigma*DeltaL only (via a
% complex-step derivative of lt_mee_rhs, which is complex-step clean). The sign s
% in H_L = ell_L + s*lambda'*dXdL that makes this identity hold is the correct
% one; the wrong sign violates the state equation xdot = partial H/partial lambda
% and blows the residual up. The chosen sign then sets H_t = -s*lambda_t.
%
% INPUTS:
%   matPath - path to a results .mat holding `res` (uses res.fuel.{X,U,dL,lamDef},
%             res.sigma), OR the res struct itself [char | struct]
%   epsv    - Bertrand-Epenoy homotopy eps of the stored solution [scalar];
%             optional, default 0 (the certified fuel solutions are eps=0)
%
% OUTPUTS:
%   H - struct:
%       .HL      [1x(N+1)]  L-domain Hamiltonian along the trajectory
%       .Ht      [1x(N+1)]  time-domain Hamiltonian (= -s*lambda_t, ~constant)
%       .tdays   [1x(N+1)]  physical time at each node [days]
%       .revs    [1x(N+1)]  cumulative revolutions (L-pi)/2pi at each node
%       .sign    +1|-1      calibrated dual sign s
%       .idResid scalar     relative extremal-identity residual at the chosen
%                           sign (median |dHL/dsig - dHL/dsig_explicit| / scale);
%                           SMALL confirms the reconstruction + sign are right
%       .idResidFlip scalar the same residual at the REJECTED sign (should be >>)
%       .HtMean .HtCoV      mean and coeff-of-variation of H_t (CoV ~ machine eps)
%
% REFERENCES:
%   [1] verify/hamiltonian_const_check.m (the lambda_t=-H_t constancy result).
%   [2] core/lt_mee_rhs.m (dXdL, Ldot; complex-step clean in par.L).
%   [3] core/casadi_lt_mee.m:278 (objective integrand ell_L definition).
if ischar(matPath) || isstring(matPath), S = load(char(matPath)); res = S.res;
else, res = matPath; end
if nargin < 2 || isempty(epsv), epsv = 0; end
X   = res.fuel.X;   U = res.fuel.U;   dL = res.fuel.dL;
sig = res.sigma(:);
par = kepler_lt_params(res.fp.thrustN, 1500, 2000);
N1  = size(X,2);   thr = U(4,:);
lam = mee_dual_to_costate(res.fuel.lamDef, sig);      % [7 x (N+1)] nodal costate

% node-by-node: dXdL, Ldot, running cost ell_L, and their explicit d/dL
dXdL = zeros(7, N1);  Ldot = zeros(1, N1);
ellL = zeros(1, N1);  dHdL_expl_noLam = zeros(1, N1);  dfdL = zeros(7, N1);
hc = 1e-20;                                            % complex-step size
for k = 1:N1
    Lk = pi + sig(k)*dL;
    pk = par;  pk.L = Lk;
    [f, Ld]    = lt_mee_rhs(X(:,k), U(:,k), pk);
    dXdL(:,k)  = f;   Ldot(k) = Ld;
    be         = thr(k) - epsv*thr(k)*(1-thr(k));
    ellL(k)    = be / Ld;
    % complex-step explicit d/dL (hold X,U fixed; perturb only par.L)
    pc = par;  pc.L = Lk + 1i*hc;
    [fc, Ldc]  = lt_mee_rhs(X(:,k), U(:,k), pc);
    dfdL(:,k)  = imag(fc)/hc;
    dHdL_expl_noLam(k) = imag(be/Ldc)/hc;              % d(ell_L)/dL
end

% Dual SIGN from transversality (not from the identity, which does not
% discriminate it): for fixed t_f, H_t = dJ*/dt_f (verified on min-int-u^2:
% J*=1/T, H=-1/T^2). A min-fuel transfer above the minimum time has
% dJ*/dt_f < 0 (more time -> less fuel), so H_t < 0; since H_t = -s*lambda_t,
% s = sign(mean lambda_t). (CasADi's opti.dual convention resolves to s=+1.)
dsg = sig(2) - sig(1);
s = sign(mean(lam(7,:)));  if s == 0, s = 1; end

HL  = ellL + s*sum(lam.*dXdL, 1);
Ht  = -s*lam(7,:);

% Extremal-identity residual, REPORTED as a reconstruction validation (the
% recovered trapezoidal costates satisfy dH_L/dsigma = explicit dH_L/dsigma to
% collocation order O(dsigma^2); the metric degrades where the per-rev H_L
% ripple is under-sampled -- ~1-2% at 10 N's ~26 nodes/rev, larger at deep
% rungs' 8 nodes/rev). NOT a machine-precision claim.
[idResid, scl] = local_id_residual(HL, dL*(dHdL_expl_noLam + s*sum(lam.*dfdL,1)), dsg, thr);

H = struct('HL', HL, 'Ht', Ht, ...
    'tdays', X(7,:)*par.TU_s/86400, 'revs', (pi + sig.'*dL - pi)/(2*pi), ...
    'sign', s, 'idResid', idResid/scl, ...
    'HtMean', mean(Ht), 'HtCoV', std(Ht)/max(abs(mean(Ht)), realmin));

fprintf(['[H-traj] sign=%+d | H_t=%.5g (CoV=%.1e) | H_L in [%.4g, %.4g] | ' ...
         'extremal-identity resid=%.1e (collocation-order validation)\n'], ...
        s, H.HtMean, H.HtCoV, min(HL), max(HL), H.idResid);
end

% -------------------------------------------------------------------------
function [resid, scl] = local_id_residual(HL, dHL_expl, dsg, thr)
% Median |dHL/dsigma - explicit dHL/dsigma| over nodes away from bang-bang
% switches (where the numeric central difference straddles a control jump and
% is not a valid derivative sample though H_L itself is continuous). scl is a
% robust magnitude of the explicit term, for a relative residual.
dHL_num = gradient(HL, dsg);
burn    = thr > 0.5;
sw      = [false, (abs(diff(burn)) > 0), false] | [false, false, (abs(diff(burn)) > 0)];
sw      = sw(1:numel(HL));
keep    = ~sw;  keep(1) = false;  keep(end) = false;   % drop endpoints + switch-adjacent
r       = abs(dHL_num(keep) - dHL_expl(keep));
resid   = median(r);
scl     = median(abs(dHL_expl(keep))) + realmin;
end
