function E = pmp_objective_error(R, opts)
% PMP_OBJECTIVE_ERROR  First-order objective error implied by a PMP residual map.
%
% Converts "the mesh is mis-allocated" into "the mesh costs X kg", which is the
% question the whole study opened with and which no earlier measurement could
% answer.
%
% THE ARGUMENT. The discrete solution satisfies its own trapezoidal defect
% constraints exactly. Measured against the EXACT flow, each interval instead
% carries a residual
%
%     d_k = X(:,k+1) - x_exact(k+1 | X(:,k), u)
%
% so the discrete solution is the exact solution of a problem whose defect
% constraints are displaced by d_k. The multiplier of defect k is, by
% definition, the sensitivity of the optimal objective to exactly that
% displacement:  lamDef(:,k) = dJ*/d(defect k). Therefore, to first order,
%
%     Delta J  ~  sum_k  lamDef(:,k)' * d_k
%
% No extra adjoint propagation is required -- the NLP already computed the
% sensitivities, which is the whole point of using the multipliers rather than
% integrating a fresh adjoint backwards (which would hit the same instability
% that forced per-interval propagation in the first place).
%
% THE SIGNED SUM AND THE ABSOLUTE SUM ANSWER DIFFERENT QUESTIONS, and reporting
% only one would hide the more interesting result:
%   sum|lam'd|   the no-cancellation bound -- what the error would be if every
%                interval's contribution pushed the same way.
%   |sum lam'd|  the actual first-order estimate, after cancellation.
% Their RATIO measures how much the per-interval errors cancel. This is a
% direct test of the H2 hypothesis, which argued that alternating-sign errors
% across many switches partially cancel in an integrated quantity. H2 was
% argued and never measured; this measures it.
%
% CONVERSION TO MASS. The objective is the fuel integral, and the mass state
% obeys mdot = -(Tmax/c)*thr, so propellant and objective are proportional:
% (m0 - m_f) = (Tmax/c) * J. The constant is taken from the SOLUTION ITSELF
% rather than re-derived, so the conversion cannot drift from the transcription.
%
% LIMITS, STATED PLAINLY. This is a FIRST-ORDER estimate. It assumes the
% residuals act as small perturbations whose effects superpose linearly. With
% ~200 intervals each missing by ~4e-4 that assumption is worth testing, not
% asserting -- which is what the mesh-redistribution experiment is for: predict
% a mass change here, then change the mesh and see whether it appears.
%
% INPUTS:
%   R    - output of pmp_residual_mee, carrying .dx, .lamDef, .par, .X [struct]
%   opts - struct (optional): .verbose [default true]
%
% OUTPUTS:
%   E - struct:
%       .dJ_signed    sum_k lamDef(:,k)'*d_k        (first-order estimate)
%       .dJ_abs       sum_k |lamDef(:,k)'*d_k|      (no-cancellation bound)
%       .cancelRatio  .dJ_abs / |.dJ_signed|  -- 1 means no cancellation at all
%       .perCell      [1xN] the per-interval contributions
%       .dMass_kg     mass error implied by .dJ_signed [kg]
%       .dMassAbs_kg  mass error implied by .dJ_abs [kg]
%       .mfKg         the solution's own final mass [kg]
%       .relMass      |.dMass_kg| / .mfKg
%
% REFERENCES:
%   [1] Nocedal & Wright, "Numerical Optimization," 2nd ed., sec. 12.9
%       (multipliers as sensitivities of the optimal value).
%   [2] pmp_residual_mee.m (produces .dx and .lamDef).

if nargin < 2, opts = struct(); end
verbose = local_default(opts, 'verbose', true);

assert(isfield(R,'dx') && isfield(R,'lamDef'), 'pmp_objective_error:input', ...
    'R must carry .dx and .lamDef (re-run pmp_residual_mee after the signed-residual change)');

dx = R.dx;  ld = R.lamDef;
N  = min(size(dx,2), size(ld,2));
dx = dx(:,1:N);  ld = ld(:,1:N);
good = all(isfinite(dx), 1) & all(isfinite(ld), 1);

per = nan(1, N);
per(good) = sum(ld(:,good) .* dx(:,good), 1);

E.perCell    = per;
E.dJ_signed  = sum(per(good));
E.dJ_abs     = sum(abs(per(good)));
E.cancelRatio = E.dJ_abs / max(abs(E.dJ_signed), realmin);

% --- objective -> mass, with the constant taken from the solution -----------
m0kg = R.par.m0kg;
mfND = R.X(6,end);
E.mfKg = mfND * m0kg;
propKg = m0kg - E.mfKg;                       % propellant actually used
% J for this solution, in the same units the multipliers price:
%   J = int thr dt, and prop = (Tmax/c)*J, so dProp = (Tmax/c)*dJ.
% Recover (Tmax/c) empirically from the solution rather than re-deriving it.
thr = R.U(4,:);  t = R.X(7,:);
Jsol = trapz(t, thr);
kTc  = propKg / max(Jsol, realmin);           % kg per unit J
E.J_solution = Jsol;
E.kg_per_J   = kTc;
E.dMass_kg    = kTc * E.dJ_signed;
E.dMassAbs_kg = kTc * E.dJ_abs;
E.relMass     = abs(E.dMass_kg) / max(E.mfKg, realmin);

if verbose
    fprintf('\n=== objective error implied by the PMP residual ===\n');
    fprintf('  intervals used            : %d of %d\n', sum(good), N);
    fprintf('  sum |lam''d|  (no cancel)  : %.6e\n', E.dJ_abs);
    fprintf('  |sum lam''d|  (actual)     : %.6e\n', abs(E.dJ_signed));
    fprintf('  cancellation ratio        : %.1f   (1 = none; large = strong)\n', E.cancelRatio);
    fprintf('  ---\n');
    fprintf('  solution J = int thr dt   : %.6f   -> %.6f kg per unit J\n', Jsol, kTc);
    fprintf('  m_f                       : %.6f kg\n', E.mfKg);
    fprintf('  implied mass error        : %+.6f kg  (%.2e relative)\n', E.dMass_kg, E.relMass);
    fprintf('  no-cancellation bound     : %.6f kg\n', E.dMassAbs_kg);
    fprintf('  NOTE: first-order, and the overall SIGN depends on the multiplier\n');
    fprintf('        convention -- the magnitude is the defensible part until the\n');
    fprintf('        mesh-redistribution experiment fixes the sign empirically.\n');
end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
