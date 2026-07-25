function D = diag_t1_beta(matPath)
% DIAG_T1_BETA  Campaign-B test T1, done properly: full-Lagrangian beta-
% stationarity from the RAW nlpsol multipliers, compared term-by-term against
% the analytic primer vector -- the decisive discriminator for the standing
% PMP primer anomaly (10-60 deg misalignment, eccentricity-correlated).
%
% WHY THIS TEST. process/DESIGN_dual_map.md's T1 was specified but never built:
% it asks whether the returned duals assemble into a valid KKT stationarity
% point at all. Two things have changed since that spec was written:
%   (a) the production transcription is now MEE/L-domain (casadi_lt_mee), not
%       the cScale-augmented Cartesian casadi_lt_2body the spec was written
%       against -- and the anomaly REPRODUCES in MEE, which has none of the
%       structural features the spec's hypotheses blamed (no per-node slack
%       state; DeltaL is one free scalar). So the spec's suspect list is
%       already refuted by the MEE evidence.
%   (b) verify/sosc/sosc_recover_kkt.m already recovers the raw low-level
%       multiplier vector opti.lam_g (NOT opti.dual) plus grad_f and the full
%       constraint Jacobian -- exactly the tooling the spec proposed building.
%
% So T1 is now cheap, and it can be sharpened from a yes/no into a term-naming
% test. The analytic beta-stationarity of the sigma-domain Hamiltonian
% (mee_primer_switch.m header, derivation re-checked independently) predicts
%
%   dLagrangian/dbeta_k  =  wq(k) * (DeltaL*thr_k/Ldot_k^2) * primerVec_k
%
% with wq(k) the trapezoid nodal weight. The LHS is computed here by AD from
% the real NLP using the COMPLETE multiplier set (defect + betaNorm + Ldot
% guard + every box), the RHS from mee_primer_switch's reconstruction. The
% tangential (radial-projected-out) parts must agree: the ||beta||=1
% multiplier is purely radial and drops out of both sides.
%
%   - LHS tangential ~ 0 and RHS tangential ~ 0  => no anomaly on this row.
%   - LHS tangential ~ 0 but RHS tangential large => the duals ARE a valid KKT
%     certificate and the analytic primer is MISSING A TERM; the residual
%     LHS-RHS then names which one (the Ldot-guard multiplier contribution is
%     reported separately for exactly this reason -- it is control-dependent,
%     so unlike every box bound it CAN rotate the primer).
%   - LHS tangential large => the recovered multipliers are not a stationarity
%     point; escalate per DESIGN_dual_map.md sec 6.
%
% Also reported: opti.dual vs raw lam_g agreement on the defect rows (the
% spec's suspected extraction-convention corruption, now directly measurable),
% and the correlation of primer misalignment with orbital eccentricity (the
% characterization the anomaly has always been described by, never measured).
%
% INPUTS:
%   matPath - certified fuel row to test [char]; default results/MEE_M2_10N.mat
%
% OUTPUTS:
%   D - struct of diagnostics (also saved to results/dual_anomaly/):
%       .statResidP .statResidM .sLag  - KKT stationarity ||grad_x L||_inf for
%                                        L = f +/- lam'g, and the chosen sign
%       .tanAD_med .tanAD_max          - tangential |dL/dbeta| on burn nodes,
%                                        AD / full multiplier set (normalized)
%       .tanAN_med .tanAN_max          - same from the analytic primerVec
%       .primerDegMed                  - verify_pmp_mee's own primer angle
%       .lamLdotMaxAbs .lamLdotActiveN - Ldot-guard multiplier magnitude/count
%       .dualVsLamgRelMax              - max rel diff opti.dual vs raw lam_g
%                                        on the defect rows (global sign
%                                        resolved), the extraction check
%       .eccCorr                       - corr(primer angle, eccentricity)
%       .missRelMed                    - median |AD - analytic| / |AD| after
%                                        removing the Ldot-guard contribution
%
% REFERENCES:
%   [1] process/DESIGN_dual_map.md (Campaign B; T1 spec this implements).
%   [2] verify/mee_primer_switch.m (the analytic primerVec being tested).
%   [3] verify/sosc/sosc_recover_kkt.m (raw lam_g recovery precedent).
import casadi.*
if nargin < 1 || isempty(matPath)
    matPath = fullfile(module_root(), 'results', 'MEE_M2_10N.mat');
end

saved = sosc_load_row(matPath);
par   = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
opts  = struct('par', par, 'mode', 'fixedtf', 'eps', 0, ...
    'tfTarget', saved.tfTarget, 'x0', saved.X(:,1), 'xf', saved.xf, ...
    'maxIter', saved.maxIter, 'warmTight', true, 'printLevel', 0, ...
    'returnModel', true);

fprintf('=== diag_t1_beta: %s ===\n', saved.tag);
o = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts);
fprintf('warm re-solve: %s | maxDefect %.3e | drift %.3e\n', ...
    o.ipoptStatus, o.maxDefect, max(abs(o.X(:) - saved.X(:))));

opti = o.model.opti;  sol = opti.debug;
xsym = opti.x;  gsym = opti.g;  fsym = opti.f;  lsym = opti.lam_g;
xv   = full(sol.value(xsym));
lamv = full(sol.value(lsym));

Xv = o.X;  Uv = o.U;  dLv = o.dL;  N1 = size(Uv, 2);  Nseg = N1 - 1;
sg = saved.sigma(:);  dsig = diff(sg).';
burn = Uv(4,:) > 0.5;

% --- variable layout (Opti stacks in declaration order: X, then U, then dL).
% ASSERTED against the primal, not assumed.
nX   = 7 * N1;
idxU = nX + (1:4*N1);
assert(max(abs(xv(idxU) - Uv(:))) < 1e-10, ...
    'diag_t1_beta: Opti variable-layout assumption failed (X then U)');

% --- KKT stationarity with the FULL raw multiplier set -----------------------
Fk = Function('Fk', {xsym, lsym}, {gradient(fsym, xsym), jacobian(gsym, xsym)});
[gfD, AD_] = Fk(xv, lamv);
gf = full(gfD);  A = sparse(AD_);
gLp = gf + A.' * lamv;
gLm = gf - A.' * lamv;
D.statResidP = norm(gLp, inf);
D.statResidM = norm(gLm, inf);
if D.statResidP <= D.statResidM, sLag = +1; gL = gLp; else, sLag = -1; gL = gLm; end
D.sLag = sLag;
fprintf(['KKT stationarity ||grad_x L||_inf : (f+lam''g) %.3e | (f-lam''g) ' ...
         '%.3e -> convention s=%+d\n'], D.statResidP, D.statResidM, sLag);

% --- LHS: tangential dL/dbeta from AD, complete multiplier set ---------------
gBeta = zeros(3, N1);
for k = 1:N1
    gBeta(:,k) = gL(nX + (k-1)*4 + (1:3));
end
tanAD = zeros(1, N1);  nrmAD = zeros(1, N1);
for k = 1:N1
    b = Uv(1:3,k) / norm(Uv(1:3,k));
    v = gBeta(:,k);
    tanAD(k) = norm(v - (v.'*b)*b);
    nrmAD(k) = norm(v);
end
D.tanAD_med = median(tanAD(burn));
D.tanAD_max = max(tanAD(burn));
fprintf(['T1 LHS  (AD, all multipliers) tangential |dL/dbeta| on burns: ' ...
         'median %.3e  max %.3e   (|dL/dbeta| median %.3e)\n'], ...
    D.tanAD_med, D.tanAD_max, median(nrmAD(burn)));

% --- RHS: the analytic primer ------------------------------------------------
lamBar = mee_dual_to_costate(o.lamDef, sg);
ver    = verify_pmp_mee(o, par, sg, struct('eps', 0));
D.primerDegMed = ver.primerMedianDeg;

% trapezoid nodal weights (one-sided at the ends), and the analytic prediction
wq = zeros(1, N1);
wq(1) = dsig(1)/2;  wq(N1) = dsig(Nseg)/2;
for k = 2:N1-1, wq(k) = (dsig(k-1) + dsig(k))/2; end

bestMiss = inf;  bestSgn = 1;
for sgn = [1 -1]
    [pv, ~, info] = mee_primer_switch(Xv, Uv, sgn*lamBar, sg, dLv, par);
    predk = zeros(3, N1);
    for k = 1:N1
        predk(:,k) = wq(k) * (dLv * Uv(4,k) / info.Ldot(k)^2) * pv(:,k);
    end
    miss = median(vecnorm(predk(:,burn) - gBeta(:,burn)) ./ ...
                  max(vecnorm(gBeta(:,burn)), 1e-30));
    if miss < bestMiss, bestMiss = miss; bestSgn = sgn; bestPred = predk; bestPv = pv; end
end
tanAN = zeros(1, N1);
for k = 1:N1
    b = Uv(1:3,k) / norm(Uv(1:3,k));
    v = bestPv(:,k);
    tanAN(k) = norm(v - (v.'*b)*b) / max(norm(v), 1e-30);
end
D.tanAN_med = median(tanAN(burn));
D.tanAN_max = max(tanAN(burn));
D.lamSign   = bestSgn;
fprintf(['T1 RHS  (analytic primerVec) RELATIVE tangential residual on burns: ' ...
         'median %.3e  max %.3e   | verify_pmp_mee primer median %.3f deg\n'], ...
    D.tanAN_med, D.tanAN_max, D.primerDegMed);

% --- the term-naming residual: AD minus analytic -----------------------------
D.missRelMed = median(vecnorm(bestPred(:,burn) - gBeta(:,burn)) ./ ...
                      max(vecnorm(gBeta(:,burn)), 1e-30));
fprintf('AD-vs-analytic dL/dbeta relative gap on burns: median %.3e\n', D.missRelMed);

% --- suspect 1: the Ldot guard is the ONLY control-dependent constraint whose
% multiplier can rotate the primer (every box bound is on a state or on beta
% itself; betaNorm is radial). Measure it.
gi = find(strcmp({o.model.creg.label}, 'ldotGuard'), 1);
lamLdot = lamv(o.model.creg(gi).rows);
D.lamLdotMaxAbs   = max(abs(lamLdot));
D.lamLdotActiveN  = nnz(abs(lamLdot) > 1e-8 * max(1, max(abs(lamv))));
fprintf('Ldot-guard multipliers: max|lam| %.3e | # nonzero %d of %d\n', ...
    D.lamLdotMaxAbs, D.lamLdotActiveN, numel(lamLdot));

% --- suspect 2: opti.dual vs raw lam_g on the defect rows (extraction check) --
di = find(strcmp({o.model.creg.label}, 'defect'), 1);
lamDefRaw = reshape(lamv(o.model.creg(di).rows), 7, Nseg);
sc = sign(sum(sum(lamDefRaw .* o.lamDef)));  if sc == 0, sc = 1; end
den = max(abs(o.lamDef(:)));
D.dualVsLamgRelMax = max(abs(sc*lamDefRaw(:) - o.lamDef(:))) / max(den, 1e-30);
fprintf('opti.dual vs raw lam_g on defect rows: max rel diff %.3e (sign %+d)\n', ...
    D.dualVsLamgRelMax, sc);

% --- the never-measured characterization: primer angle vs eccentricity -------
ecc = sqrt(Xv(2,:).^2 + Xv(3,:).^2);
pdb = ver.primerDeg(burn);  eb = ecc(burn);
ok  = isfinite(pdb);
D.eccCorr = corr(eb(ok).', pdb(ok).');
fprintf('corr(primer angle, eccentricity) on burns: %+.3f  (ecc range %.3f-%.3f)\n', ...
    D.eccCorr, min(eb), max(eb));

outFile = fullfile(module_root(), 'results', 'dual_anomaly', ...
                   sprintf('diag_t1_beta_%s.mat', saved.tag));
save(outFile, 'D');
fprintf('saved %s\n', outFile);
end
