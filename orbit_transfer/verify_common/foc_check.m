function rep = foc_check(out, sigma, man, opts)
% FOC_CHECK  Generic AD-based first-order optimality (PMP/KKT) gate.
%
% Campaign-agnostic first-order check on a converged direct-collocation NLP
% solved via CasADi Opti. Recomputes the full Lagrangian gradient by AD
% (grad_f + A'*lam, sign auto-resolved against opti.lam_g's actual
% convention -- see LESSONS_DUAL_EXTRACTION.md), then, driven purely by the
% campaign manifest (foc_manifest) and the solver's constraint registry
% (out.model.creg), reports: KKT stationarity, the bang-bang switching-
% function sign law, tangential direction-vector stationarity, nodal
% costates (via foc_dual_to_costate), time-costate behavior, free-mass
% transversality, and singular-arc / regular-switching diagnostics. All
% verdicts are advisory (.pass is a burn-in signal, not a hard gate).
%
% INPUTS:
%   out   - solved-model struct [struct]:
%           .X          state trajectory [nx x (N+1)]
%           .U          control trajectory [nu x (N+1)]
%           .model.opti solved CasADi Opti object
%           .model.creg constraint registry, struct array with fields
%                       .label [char], .rows [1xk] (row range into opti.g);
%                       label 'defect' required, 'thrLo'/'thrHi' required
%                       when man.thrRow is nonempty
%   sigma - node grid, monotonic increasing [(N+1)x1] or [1x(N+1)]
%   man   - campaign manifest from foc_manifest [struct]
%   opts  - optional tolerances [struct], fields (all scalar, all optional):
%           .tolStat  [1e-6] KKT stationarity / tangential-residual tolerance
%           .tolSign  [99]   sign-law pass threshold, percent
%           .tolTrans [1e-3] free-mass transversality relative tolerance
%           .sdotMin  [1e-3] minimum regular-switching relative Sdot
%
% OUTPUTS:
%   rep - report struct [struct], fields:
%         .kktStatInf       max-norm full-Lagrangian stationarity residual
%         .sLag             resolved Lagrangian sign convention (+1 or -1)
%         .dirTanMax        max tangential dL/dbeta residual over burn nodes
%         .dirTanMed        median tangential dL/dbeta residual over burn nodes
%         .signPct          switching-function sign-law agreement, percent
%         .Sd               switching function samples [1 x (N+1)]
%         .lam              nodal costates [nx x (N+1)]
%         .lamTimeCoV       time-costate coefficient of variation
%         .lamTimeEnd       time-costate value at the final node
%         .lamMassEndMapped free-final-mass transversality residual (relative),
%                           Hager-mapped terminal covector -- THE GATED VALUE
%         .singularArcNodes count of >=3-node near-zero switching-function runs
%         .sdotMinRel       minimum relative |Sdot| across detected switches
%         .nSwitches        number of burn/coast sign changes
%         .horizonNote      horizon-condition caveat text [char]
%         .checksRun        cellstr of which checks executed
%         .pass             advisory overall verdict [logical]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md sec 3
%       (opti.lam_g vs opti.dual sign/convention pitfall).
%   [2] earth_elliptic_to_geo/direct/results/dual_anomaly/diag_t1_beta.m
%       (tangential dL/dbeta stationarity diagnostic, the method precedent).
%   [3] verify_common/OPTIMALITY_CERTIFICATION.md sec A2 (generic FOC gate
%       design: manifest-driven routing, advisory verdict).
%   [4] Hager, W.W., "Runge-Kutta Methods in Optimal Control and the
%       Transformed Adjoint System," Numer. Math. 87, 247-282, 2000 (the
%       mapped/transformed terminal covector gated in section (6) below).

import casadi.*
if nargin < 4, opts = struct(); end
gd = @(f,v) fcdef(opts, f, v);
tolStat = gd('tolStat',1e-6); tolSign = gd('tolSign',99);
tolTrans = gd('tolTrans',1e-3); sdotMin = gd('sdotMin',1e-3);

opti = out.model.opti;  creg = out.model.creg;  sol = opti.debug;
X = out.X;  U = out.U;  N1 = size(X,2);  Nseg = N1-1;
nx = man.nx;  nu = man.nu;
assert(size(X,1)==nx && size(U,1)==nu, 'foc_check: X/U shape vs manifest');
sg = sigma(:);

xall   = full(sol.value(opti.x));
lamAll = full(sol.value(opti.lam_g));
nX = nx*N1;
% Variable-layout assert (X block then U block, column-major) -- never assume:
assert(max(abs(xall(1:nX) - X(:))) < 1e-8, 'foc_check: X layout');
assert(max(abs(xall(nX+(1:nu*N1)) - U(:))) < 1e-8, 'foc_check: U layout');
uix = @(rows,k) nX + (k-1)*nu + rows;

% Manifest semantic guards (cheap, catch a wrong manifest outright):
if ~isempty(man.massRow)
    assert(all(diff(X(man.massRow,:)) <= 1e-6), 'foc_check: massRow not nonincreasing -- manifest wrong?');
end
if ~isempty(man.timeRow)
    assert(all(diff(X(man.timeRow,:)) >= -1e-6), 'foc_check: timeRow not nondecreasing -- manifest wrong?');
end

% --- (1) full-Lagrangian KKT stationarity, sign auto-resolved ----------------
Fk = Function('Fk', {opti.x, opti.lam_g}, {gradient(opti.f,opti.x), jacobian(opti.g,opti.x)});
[gfD, AD_] = Fk(xall, lamAll);
gf = full(gfD); A = sparse(AD_);
gP = gf + A.'*lamAll;  gM = gf - A.'*lamAll;
if norm(gP,inf) <= norm(gM,inf), s = +1; gL = gP; else, s = -1; gL = gM; end
rep = struct('kktStatInf', norm(gL,inf), 'sLag', s);
checks = {'kktStat'};

rowsOf = @(lab) [creg(strcmp({creg.label},lab)).rows];
defRows = rowsOf('defect');
assert(~isempty(defRows), 'foc_check: creg must register a ''defect'' group');

% --- (2) burn mask + throttle-block switching function -----------------------
% BANG-BANG APPLICABILITY (finding: verdict was not eps-aware, external review
% 2026-07-25). The sign law, the singular-arc test and the regular-switching
% test are all statements about an AFFINE-in-throttle cost. On an eps>0
% homotopy leg the cost is quadratic in thr, Sd is the derivative of the
% SMOOTHED objective rather than a switching function, an interior throttle can
% legitimately give Sd=0, and the sign of Sd carries no bang-bang implication.
% Previously these three checks ran anyway and were folded into rep.pass --
% e.g. the ELFO energy seed (eps=1) reported "sign law 100%" and an ADVISORY
% PASS, a number with no meaning at that eps. They are now SKIPPED (NaN,
% reported as '--') unless the throttle cost is affine.
%
% Resolution order: opts.eps (if given, affine iff eps==0) overrides
% opts.throttleCostKind, which overrides the manifest's own field.
epsGiven = fcdef(opts, 'eps', []);
if ~isempty(epsGiven)
    thrKind = 'affine';  if epsGiven ~= 0, thrKind = 'quadratic'; end
else
    thrKind = fcdef(opts, 'throttleCostKind', ...
                    fcdef(man, 'throttleCostKind', 'affine'));
end
rep.throttleCostKind = thrKind;
bangBang = ~isempty(man.thrRow) && strcmp(thrKind, 'affine');
rep.bangBangChecksRun = bangBang;

if ~isempty(man.thrRow)
    thr = U(man.thrRow,:);  burn = thr > 0.5;
    lamNB = lamAll;  lamNB([rowsOf('thrLo'), rowsOf('thrHi')]) = 0;
    gNB = gf + s*(A.'*lamNB);           % Lagrangian gradient sans throttle-bound duals
    Sd = zeros(1,N1);
    for k = 1:N1, Sd(k) = gNB(uix(man.thrRow,k)); end
    rep.Sd = Sd;
else
    burn = true(1,N1);  rep.Sd = [];   % min-time: all-burn, no throttle row
end
if bangBang
    rep.signPct = 100*mean((rep.Sd < 0) == burn);   % S<0 <=> full thrust
    checks{end+1} = 'signLaw';
else
    rep.signPct = NaN;
end

% --- (3) minimum condition, direction part: tangential dL/dbeta --------------
if ~isempty(man.dirRows)
    tanAbs = nan(1,N1);
    for k = 1:N1
        b = U(man.dirRows,k);  b = b/norm(b);
        v = gL(uix(man.dirRows,k));
        tanAbs(k) = norm(v - (v.'*b)*b);
    end
    rep.dirTanMax = max(tanAbs(burn));  rep.dirTanMed = median(tanAbs(burn));
    checks{end+1} = 'dirTangential';
else
    rep.dirTanMax = 0;  rep.dirTanMed = 0;
end

% --- (4) nodal costates (sign-resolved) --------------------------------------
assert(numel(defRows) == nx*Nseg, ...
    'foc_check: defect group has %d rows, expected nx*Nseg = %d -- creg mis-registered?', ...
    numel(defRows), nx*Nseg);
LamDef = reshape(lamAll(defRows), nx, Nseg);
% The stationarity combination (one-sided at the endpoints). Endpoint-VALUED
% checks do not use this -- see (6): they read the mapped covector off gL.
rep.lam = foc_dual_to_costate(s*LamDef, sg);
checks{end+1} = 'costates';

% --- (5) time-costate behavior (Hamiltonian conditions in dual form) ---------
if ~isempty(man.timeRow)
    lt = rep.lam(man.timeRow,:);
    rep.lamTimeCoV = std(lt)/max(abs(mean(lt)),1e-30);
    rep.lamTimeEnd = lt(end);
    checks{end+1} = 'lamTime';
else
    rep.lamTimeCoV = NaN;  rep.lamTimeEnd = NaN;
end

% --- (6) transversality: free final mass -> lam_m(tf)=0 ----------------------
% MAPPED TERMINAL COVECTOR (Hager 2000 [4]). The exact discrete object is
%
%   lamHat_f = (h_N/2)*L_x(N+1) + (I - (h_N/2)*f_x(N+1))' * Lambda_N
%
% whose mass component is EXACTLY ZERO under the discrete KKT system whenever
% the final mass is genuinely free (man.massFreeAtTf) with no terminal cost or
% constraint touching it. Rather than re-derive f_x/L_x by hand (forbidden --
% AD from the solved opti.f/opti.g, never a hand derivation), read it directly
% off the assembled Lagrangian gradient at the terminal-node mass index:
% grad(f)'s only dependence on X(:,N+1) is the (h_N/2)*L quadrature term, and
% A's only dependence there is the last defect block's Jacobian, so that one
% entry of gL IS s*lamHat_f_mass. This is also the most robust form -- gL
% already sums any OTHER constraint touching X(massRow,N+1) (an active box
% bound, say), which a hand reconstruction would silently miss.
%
% CONSEQUENCE, stated plainly: this quantity is one ENTRY of the very vector
% whose max-norm is rep.kktStatInf. It is therefore bounded by it (up to the
% normalisation below) and CANNOT fail unless kktStat already has. It is a
% readable projection of the master residual, not an independent test -- the
% same structural situation as the direction check at (3). Both are listed in
% rep.derivedFromKKT and the report groups them under the master residual
% rather than beside the structural checks. Do not present either as
% independent evidence.
%
% Superseded and REMOVED 2026-07-25 (kept only in the record, not the code):
% the raw one-sided interval dual and its linear extrapolation. Both estimated
% lambda_m(t_f) from Lambda_N, which behaves like a sample near
% sigma_end - h_N/2 rather than a value at the endpoint -- an O(h) endpoint
% REPRESENTATION offset, not evidence about the solution. They produced a
% spurious 1.265e-03 "finding" on the earth 2.5 N row that the mapped covector
% resolves to 2.27e-18. History: verify_common/doc/first_order_checks.tex and
% OPTIMALITY_CERTIFICATION.md.
if ~isempty(man.massRow) && man.massFreeAtTf
    scale = max(abs(rep.lam(man.massRow,:)));
    xix   = @(rows,k) (k-1)*nx + rows;
    rep.lamMassEndMapped = abs(gL(xix(man.massRow, N1))) / max(scale,1e-30);
    checks{end+1} = 'transversality';
else
    rep.lamMassEndMapped = NaN;
end

% --- (7) singular arcs + regular switching (Sdot != 0) -----------------------
% MESH NORMALIZATION (finding I1, external review 2026-07-25). Sd_k carries the
% local trapezoid quadrature weight w_k, so Sd ~ w_k * S(sigma_k): it is NOT a
% faithful sampling of S. The old statistic differenced Sd across a switch
% WITHOUT dividing by the step, so it shrank with mesh refinement even for a
% perfectly transversal crossing -- and PSR-refined meshes densify exactly AT
% switches, which is where it bites hardest. Both reviewers judged the readings
% it produced (1.4e-7 tulip, 4.5e-5 ELFO) unable to diagnose grazing.
%
% Corrected statistic: deweight, take a true divided difference in sigma, then
% non-dimensionalize by the total span and a coast-side scale:
%     Shat_k = Sd_k / w_k
%     D_i    = |Shat(k_i+1) - Shat(k_i)| / (sigma(k_i+1) - sigma(k_i))
%     R_i    = (sigma_end - sigma_1) * D_i / Sref
% R is dimensionless and invariant to local refinement for a smooth crossing.
% Sref excludes a one-node neighbourhood of every switch so the scale is not
% contaminated by the crossing itself. The legacy value is retained as
% (The legacy un-normalised statistic was retained during burn-in for
% attribution and has been REMOVED now that the confound is confirmed: it
% differed by 8 orders on the tulip flagship, 1.376e-07 vs 2.701e+01.)
%
% CAVEAT that survives this fix (GPT): no threshold on a SINGLE mesh is
% defensible. Regularity should be asserted only if R stays above the floor on
% two meshes; the sdotMin gate is provisional until recalibrated.
rep.singularArcNodes = NaN;  rep.sdotMinRel = NaN;
rep.nSwitches = 0;  rep.Sdeweighted = [];
if ~isempty(man.thrRow)
    hs = diff(sg(:).');
    w  = [hs(1)/2, (hs(1:end-1) + hs(2:end))/2, hs(end)/2];   % nodal weights [1xN1]
    Shat = rep.Sd ./ w;
    rep.Sdeweighted = Shat;

    swI = find(diff(burn) ~= 0);  rep.nSwitches = numel(swI);

    % coast-side scale, excluding a one-node neighbourhood of each switch
    nearSw = false(1,N1);
    for q = 1:numel(swI)
        nearSw(max(swI(q)-1,1):min(swI(q)+2,N1)) = true;
    end
    cand = abs(Shat(~burn & ~nearSw));
    cand = cand(cand > 0);
    if isempty(cand), cand = abs(Shat(~burn)); cand = cand(cand > 0); end
    if isempty(cand), Sref = max(abs(Shat)); else, Sref = median(cand); end
    if isempty(Sref) || isnan(Sref) || Sref == 0, Sref = max(max(abs(Shat)), 1e-30); end
end
if bangBang
    rep.singularArcNodes = 0;
    nearZ = abs(rep.Sdeweighted) < max(1e-6*Sref, 1e-14);   % deweighted (I1)
    runL = 0;
    for k = 1:N1
        if nearZ(k), runL = runL+1; else, runL = 0; end
        if runL >= 3, rep.singularArcNodes = rep.singularArcNodes + 1; end
    end
    if ~isempty(swI)
        ka = max(swI,1);  kb = min(swI+1,N1);
        D  = abs(rep.Sdeweighted(kb) - rep.Sdeweighted(ka)) ./ max(sg(kb).' - sg(ka).', 1e-30);
        rep.sdotMinRel = min( (sg(end)-sg(1)) * D / Sref );
    end
    checks = [checks, {'singularArc','sdotRegular'}];
end

% --- (8) horizon note (G4): named, honest ------------------------------------
switch man.horizonKind
    case 'fixedtf'
        rep.horizonNote = 'fixed t_f: H=const generally nonzero; constancy via lamTimeCoV';
    case 'freetf-cscale'
        rep.horizonNote = ['free t_f via cScale state: horizon condition enters the ' ...
            'cScale adjoint rows, already inside kktStatInf; value-form H(tf) check ' ...
            'reported via lamTimeEnd (informational, derivation pending)'];
    otherwise
        rep.horizonNote = 'no horizon check applicable';
end
rep.checksRun = checks;
% Which reported lines are PROJECTIONS of rep.kktStatInf rather than
% independent tests (see (3) and (6)). Consumers should not count these as
% separate evidence; foc_report groups them under the master residual.
rep.derivedFromKKT = {'dirTangential','transversality'};

% --- advisory verdict (REPORT-ONLY burn-in) ----------------------------------
% Every ok* is NaN-tolerant: a skipped check (manifest field empty, or the
% bang-bang family skipped because the throttle cost is not affine) must
% neither pass nor fail the verdict -- it drops out.
okSign  = isnan(rep.signPct)          || rep.signPct >= tolSign;
okTrans = isnan(rep.lamMassEndMapped) || rep.lamMassEndMapped <= tolTrans;
okSdot  = isnan(rep.sdotMinRel)       || rep.sdotMinRel > sdotMin;
okSing  = isnan(rep.singularArcNodes) || rep.singularArcNodes == 0;
rep.pass = rep.kktStatInf <= tolStat && rep.dirTanMax <= tolStat && ...
           okSign && okTrans && okSing && okSdot;

% --- (opt) conjugate-point test (second-order, ADVISORY) ---------------------
% Available since migration #5: pass opts.msInfo (an ms_bvp info struct with
% .PHI/.Y/.tGrid, e.g. from ms_tfmin with conjTest/keepSTMs) and opts.msFlow
% (state-rows dynamics handle). Runs the free-time quotiented Jacobi test
% from costate_common. ADVISORY: reported as rep.conj, NOT folded into
% rep.pass -- existing campaign reports stay bit-identical when unused.
if isfield(opts,'msInfo') && ~isempty(fcdef(opts,'msInfo',[]))
    if isempty(which('ms_conjugate_test'))
        addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'costate_common'));
    end
    rep.conj = ms_conjugate_test(opts.msInfo, ...
                   struct('flow', opts.msFlow));
end
end

function v = fcdef(o, f, dflt)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = dflt; end
end
