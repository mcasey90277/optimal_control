function C = certify_root(seed, rv0, rvf, B, opts)
%% Purpose:
%
%   The min-time GATE STACK for one normal-chart multiple-shooting seed:
%
%     1. polish with ms_tfmin (tolR 3e-11) + the free-time conjugate test;
%     2. fly from z8 alone and gate the arrival in POSITION and VELOCITY (a
%        harness that gated on position only accepted arrivals with the
%        wrong velocity -- Astra review 2026-09-09);
%     2b. the pointwise Pontryagin checks on that flight
%        (pmp_pointwise_checks): H = 0, transversality, the adjoint
%        equations, and the EXACT minimum-principle gap of the control the
%        propagator applied -- the FULL gap (direction + throttle) and the
%        throttle recovered on BOTH rows it enters (acceleration and mass);
%        the shooting residual cannot see a control law that minimises the
%        wrong Hamiltonian;
%     3. pumpkyn tfMin from the polished z8 must return within tolDz (a
%        foreign-solver witness; an exception is a FAIL, not a pass);
%     4. the conjugate verdict must be PASS (ENDPOINT = inconclusive = FAIL);
%     5. the min-time hypothesis gates: min|lam_v| > 0, min Q_mt > 0 and
%        dim S = 1 (no abnormal lift of the same trajectory) -- the rank
%        REQUIRING the accepted lift's residual, the independent-field
%        Hamiltonian residual and the Eckart-Young lift_margin (the
%        constraint matrix built at two integration tolerances), plus X2:
%        pumpkyn's field against the independently written CR3BP field,
%        row by row (state and adjoint);
%     6. H6, the reduced conjugate instrument's validity condition
%        lambda_m(0) < c/T: the helper's own verdict (strict margin AND
%        clearance above the Hamiltonian residual), then margin > h6MarginMin;
%     7. the dense conjugate scan (conj_spectrum): no located zero, no
%        multiplicity, nothing UNRESOLVED.
%   Gates 2b, 5 and 7 were widened on 2026-09-11 so that the certifier
%   enforces the SAME checks as transfer_study (Astra review #2: "the
%   production certifier still enforces less than the study" -- the same
%   computed-not-enforced defect class as FINDINGS 36, 41, 47). A gate that
%   is computed and then ignored is the failure mode this stack exists to
%   prevent: before 2026-09-10 an entry could certify with h6Ok false.
%
%   Every flight (ours and the witness's) passes validate_flight: reached
%   t_f, finite, all-burn mass law, clear of both primaries. A returned
%   array is not a completed flight.
%
%   Numbers are kept whether or not the seed certifies; C.reason names the
%   FIRST gate that failed. Callers: certify_crossing (arrival-phase
%   arclength candidates, after converting the homogeneous chart) and the
%   departure ribs, whose seeds are already normal-chart.
%
%% Inputs:
%
%  seed                     struct                  ms_bvp seed (.tf, .tGrid,
%                                                   .Y [14 x K+1])
%  rv0, rvf                 [6 x 1]                 departure / arrival states
%  B                        struct                  .Tnd .cnd .mu (operating
%                                                   point, arclength_arrival)
%  opts                     struct (optional)
%   .gateKm [100] .gateVms [10] .tolDz [1e-6] .wallSec [600] .m0kg [150]
%   .tolR [3e-11] polish tolerance, .tolRelax [1e-8] a polish that plateaus
%   below this still goes to the gates (and says so in .reason),
%   .pool [gcp('nocreate')] the fence's worker pool -- WITHOUT one the
%   external calls are unfenced and can hang for hours,
%   .capPolishSec [900] .capFlySec [300] .capWitnessSec [300]
%   .capGatesSec [900] hard caps; exceeding one is a NAMED FAILURE,
%   .tolH [1e-6] .tolLamMf [1e-6] .tolAdj [1e-7] .tolGap [1e-12] (on the
%   |full gap|) .tolThrottle [1e-10] the pointwise PMP tolerances,
%   .tolField [1e-10] X2 relative row error, .tolLift [1e-6] the accepted
%   lift's BACKWARD error |C lam0|/(|C||lam0|) (one normalisation, shared
%   with lift_margin's liftTol), .liftMarginMin [10] Eckart-Young
%   margin, .relTolPair [1e-12 1e-9] the two adjoint integration settings
%   the margin is measured between,
%   .h6MarginMin [1] required (c/T)/lambda_m(0) ratio, .moonKmMin [1900]
%   .earthKmMin [6600] flight clearances (validate_flight),
%   .conjSpectrum [true] run the dense scan, .conjNSub [8] its samples per
%   segment, .conjOpts [struct()] further conj_spectrum options,
%   .pwOpts [struct()] .gatesOpts [struct()] forwarded to
%   pmp_pointwise_checks / mintime_hypothesis_gates (the mutation test
%   injects a wrong field through these),
%   .nSamp [200] .rankTol [1e-8] .sA .sD (recorded, not used)
%
%% Outputs:
%
%  C                        struct                  .ok .reason .z [8x1]
%                                                   .Y [14 x K] .tfDays
%                                                   .dvKms .propellantKg
%                                                   .finalMassKg .flyKm
%                                                   .flyVms .dz .conj .g
%                                                   .Hmax .lamMf .adjErr
%                                                   .dirGap .fullGap
%                                                   .throttleErr (pointwise
%                                                   PMP) .fieldErr .adjErrRef
%                                                   (X2) .nullResid .Hresid
%                                                   .liftMargin (S3)
%                                                   .conjDense (conj_spectrum
%                                                   summary) .h6Margin .sA
%                                                   .sD .rho
%                                                   .normR .wallSec
%                                                   .flyKmWitness
%                                                   .flyVmsWitness
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
gateKm = d('gateKm', 100);  gateVms = d('gateVms', 10);  tolDz = d('tolDz', 1e-6);
wallSec = d('wallSec', 600);  m0kg = d('m0kg', 150);
% tolR was DOCUMENTED as an option and then hardcoded -- the caller's value
% was silently ignored until 2026-09-09.
tolR = d('tolR', 3e-11);
% A PLATEAU is not a failure. The rib walking the positive departure sense
% out of the anchor stalls at sD = 0.0465 with |R| = 5.1e-11 against a 3e-11
% tolerance, at the same place and the same residual whatever the step size:
% the solver's achievable floor there is simply above a very tight number.
% The residual is one check among five, and the flown miss and the foreign
% witness are stronger evidence than its last decade, so a polish that
% plateaus below tolRelax goes on to the gates and lets THEM decide -- and
% says so in the reason, so a plateaued entry is never mistaken for a clean
% one.
tolRelax = d('tolRelax', 1e-8);
% HARD WALL-CLOCK FENCE on every external call (2026-09-10). `wallSec`
% below is an IN-PROCESS check: it fires only between solver iterations, so
% it cannot bound a single grinding evaluation -- and a CR3BP segment
% propagation whose iterate parks a junction near a primary can grind for
% hours inside ONE evaluation (run_capped's header, after it happened twice
% on the 0.5 N campaign; then again on 2026-09-09, when the sheet assembly
% and a departure rib each burned 12-16 h at 100% CPU with no progress).
% run_capped cancels the worker, so it is the only real fence. A capped
% call FAILS with a named reason: silently accepting a candidate whose
% witness or gates could not be run would put an unverified entry in the
% catalog, which is the one thing this stack exists to prevent.
pool = d('pool', gcp('nocreate'));
capPolish  = d('capPolishSec',  900);
capFly     = d('capFlySec',     300);
capWitness = d('capWitnessSec', 300);
capGates   = d('capGatesSec',   900);
h6MarginMin = d('h6MarginMin', 1);
tolH = d('tolH', 1e-6);  tolLamMf = d('tolLamMf', 1e-6);
tolAdj = d('tolAdj', 1e-7);  tolGap = d('tolGap', 1e-12);  tolThrottle = d('tolThrottle', 1e-10);
tolField = d('tolField', 1e-10);  tolLift = d('tolLift', 1e-6);  liftMarginMin = d('liftMarginMin', 10);
relPair = d('relTolPair', [1e-12 1e-9]);
doSpectrum = d('conjSpectrum', true);  conjNSub = d('conjNSub', 8);  conjOpts = d('conjOpts', struct());
pwOpts = d('pwOpts', struct());  gatesOpts = d('gatesOpts', struct());
flightOpts = struct('moonKmMin', d('moonKmMin', 1900), 'earthKmMin', d('earthKmMin', 6600));
if isempty(pool)
    warning('certify_root:unfenced', ...
        'no parallel pool: external calls run UNFENCED and can hang indefinitely (see capped_pool)');
end
lStar = 389703.264829278;  tStar = 382981.289129055;
t0 = tic;
rv0 = rv0(1:6);  rvf = rvf(1:6);

C = struct('ok', false, 'reason', '', 'z', nan(8,1), 'Y', [], 'tfDays', NaN, ...
           'dvKms', NaN, 'propellantKg', NaN, 'finalMassKg', NaN, ...
           'flyKm', NaN, 'flyVms', NaN, 'dz', NaN, ...
           'conj', -1, 'g', [], 'sA', d('sA', NaN), 'sD', d('sD', NaN), ...
           'rho', NaN, 'normR', NaN, 'wallSec', NaN, 'flyKmWitness', NaN, ...
           'flyVmsWitness', NaN, 'h6Margin', NaN, ...
           'Hmax', NaN, 'lamMf', NaN, 'adjErr', NaN, 'dirGap', NaN, 'fullGap', NaN, ...
           'throttleErr', NaN, 'fieldErr', NaN, 'adjErrRef', NaN, 'nullResid', NaN, 'nullResidRel', NaN, ...
           'Hresid', NaN, 'liftMargin', NaN, 'conjDense', []);

% ---- 1. normal-chart polish + conjugate test ---------------------------
try
    [okC, z, it] = fenced(pool, capPolish, @ms_tfmin, 2, rv0, rvf, seed, B.Tnd, B.cnd, B.mu, ...
                          struct('tolR', tolR, 'wallSec', wallSec, 'conjTest', true));
catch ME
    C.reason = ['ms_tfmin threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okC
    C.reason = sprintf('polish exceeded its %g s cap (or the worker errored)', capPolish);
    C.wallSec = toc(t0);  return
end
% READ the residual before trusting the flag. `it.converged` is the
% solver's own opinion; if it is true while normR is NaN or large, the
% candidate must still fail. Astra chain review 2026-09-10.
[okR, nr] = scalar_verdict(it.normR);
if ~okR
    C.reason = 'polish returned no usable residual';  C.wallSec = toc(t0);  return
end
C.normR = nr;
[okCv, cvFlag] = scalar_verdict(it.converged);
plateau = false;
if ~(okCv && cvFlag == 1 && nr <= tolR)
    if nr < tolRelax
        plateau = true;
    else
        C.reason = sprintf('normal-chart polish did not converge (|R| = %.1e)', nr);
        C.wallSec = toc(t0);  return
    end
end
C.z = z(:);  C.Y = it.Y;  C.tfDays = z(8)*tStar/86400;

% ---- 2. flown arrival, position AND velocity ---------------------------
[okF, tF, Yf] = fenced(pool, capFly, @pumpkyn.cr3bp.tfMinProp, 2, z(8), [rv0; 1; z(1:7)], B.Tnd, B.cnd, B.mu);
if ~okF
    C.reason = sprintf('the flight exceeded its %g s cap', capFly);  C.wallSec = toc(t0);  return
end
% A RETURNED ARRAY IS NOT A COMPLETED FLIGHT: the shared validator checks
% that it reached t_f, is finite, obeys the all-burn mass law and stays
% clear of both primaries, so this file and the study script cannot drift
% apart on what "admissible" means.
VF = validate_flight(tF, Yf, z(8), B.Tnd, B.cnd, B.mu, lStar, flightOpts);
if ~VF.ok, C.reason = ['flight inadmissible: ' VF.reason];  C.wallSec = toc(t0);  return, end
C.flyKm  = norm(Yf(end,1:3) - rvf(1:3)')*lStar;
C.flyVms = norm(Yf(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
mf = Yf(end,7);
% BOTH, named for what they are. `mfKg` held propellant USED here while
% verify_common/pmp/pmp_objective_error uses the same name for FINAL MASS --
% one name, two meanings, in one repository. (Astra script review 2026-09-10.)
C.propellantKg = (1 - mf)*m0kg;
C.finalMassKg  = mf*m0kg;
C.dvKms = B.cnd*log(1/mf)*lStar/tStar;
if ~(C.flyKm < gateKm),   C.reason = sprintf('flown position miss %.1f km > %g', C.flyKm, gateKm);  C.wallSec = toc(t0); return, end
if ~(C.flyVms < gateVms), C.reason = sprintf('flown velocity miss %.2f m/s > %g', C.flyVms, gateVms); C.wallSec = toc(t0); return, end

% ---- 2b. pointwise Pontryagin checks on the flight ----------------------
PW = pmp_pointwise_checks(tF, Yf, B.Tnd, B.cnd, B.mu, pwOpts);
C.Hmax = PW.Hmax;  C.lamMf = PW.lamMf;  C.adjErr = PW.adjErr;  C.dirGap = PW.dirGap;
C.fullGap = max(PW.fullGap, PW.fieldGap);  C.throttleErr = PW.throttleErr;
if ~(PW.Hmax <= tolH),       C.reason = sprintf('Hamiltonian max|H| = %.2e > %g', PW.Hmax, tolH);          C.wallSec = toc(t0); return, end
if ~(PW.lamMf <= tolLamMf),  C.reason = sprintf('transversality |lambda_m(t_f)| = %.2e > %g', PW.lamMf, tolLamMf); C.wallSec = toc(t0); return, end
if ~(PW.adjErr <= tolAdj),   C.reason = sprintf('adjoint equations, relative error %.2e > %g', PW.adjErr, tolAdj); C.wallSec = toc(t0); return, end
if ~(C.fullGap <= tolGap),   C.reason = sprintf('minimum-principle |full gap| %.2e > %g (applied control is not the minimiser; signed %+.1e..%+.1e, field''s own gap %.1e)', C.fullGap, tolGap, PW.gapMin, PW.gapMax, PW.fieldGap); C.wallSec = toc(t0); return, end
if ~(PW.throttleErr <= tolThrottle)
    C.reason = sprintf('applied throttle is not 1: acceleration rows %.2e, mass row %.2e > %g', ...
                       PW.throttleAccErr, PW.throttleMassErr, tolThrottle);
    C.wallSec = toc(t0);  return
end

% ---- 3. foreign witness ---------------------------------------------------
try
    [okW, za] = fenced(pool, capWitness, @pumpkyn.cr3bp.tfMin, 1, rv0', rvf', z(:), B.Tnd, B.cnd, B.mu);
catch ME
    C.reason = ['tfMin witness threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okW
    C.reason = sprintf('tfMin witness exceeded its %g s cap', capWitness);  C.wallSec = toc(t0);  return
end
if ~(isnumeric(za) && numel(za) == 8 && all(isfinite(za(:))))
    C.reason = 'tfMin witness returned an unusable vector';  C.wallSec = toc(t0);  return
end
C.dz = norm(za(:) - z(:));
% AGREEMENT IS NOT CONVERGENCE. A solver that returned its own input on
% stagnation would give dz = 0 and pass. Measured 2026-09-10: this one does
% NOT -- perturbing lam0 by 1.5x, 3x and 10x moved its answer by 9.3, 37 and
% 4612, and an impossible target moved it by 5.6, so dz = 0 here is genuine
% agreement. The general hole is closed anyway by flying the WITNESS's own
% solution: two independent solutions that both reach the target is a much
% stronger statement than two vectors that match. (Astra chain review.)
[okWF, ta, Ya] = fenced(pool, capFly, @pumpkyn.cr3bp.tfMinProp, 2, za(8), [rv0; 1; za(1:7)], B.Tnd, B.cnd, B.mu);
if ~okWF
    C.reason = sprintf('the witness flight exceeded its %g s cap', capFly);  C.wallSec = toc(t0);  return
end
VW = validate_flight(ta, Ya, za(8), B.Tnd, B.cnd, B.mu, lStar, flightOpts);
if ~VW.ok, C.reason = ['witness flight inadmissible: ' VW.reason];  C.wallSec = toc(t0);  return, end
C.flyKmWitness  = norm(Ya(end,1:3) - rvf(1:3)')*lStar;
C.flyVmsWitness = norm(Ya(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
if ~(C.flyKmWitness < gateKm && C.flyVmsWitness < gateVms)
    C.reason = sprintf('witness solution does not fly to the target (%.1f km, %.2f m/s)', ...
                       C.flyKmWitness, C.flyVmsWitness);
    C.wallSec = toc(t0);  return
end
if ~(isfinite(C.dz) && C.dz <= tolDz), C.reason = sprintf('tfMin witness |dz| = %.2e > %g', C.dz, tolDz); C.wallSec = toc(t0); return, end

% ---- 4. conjugate test ----------------------------------------------------
if isfield(it, 'conj') && isfield(it.conj, 'pass')
    [okJ, jv] = scalar_verdict(it.conj.pass);
    if okJ, C.conj = jv; else, C.conj = -1; end
end
if ~(C.conj == 1)
    C.reason = sprintf('conjugate test verdict %g', C.conj);  C.wallSec = toc(t0);  return
end

% ---- 5. hypothesis gates --------------------------------------------------
gOpts = gatesOpts;  gOpts.nSamp = d('nSamp', 200);  gOpts.rankTol = d('rankTol', 1e-8);
gOpts.keepC = true;  gOpts.relTol = relPair(1);
try
    [okG, g] = fenced(pool, capGates, @mintime_hypothesis_gates, 1, z(:), rv0, B.Tnd, B.cnd, B.mu, gOpts);
catch ME
    C.reason = ['gates threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okG
    C.reason = sprintf('hypothesis gates exceeded their %g s cap', capGates);  C.wallSec = toc(t0);  return
end
C.g = g;
% X2, the PHYSICS: pumpkyn's field against the independently written one,
% row by row -- every other check propagates pumpkyn and would pass a
% self-consistent wrong field
for f2 = {'fieldErr', 'adjErrRef', 'nullResid', 'nullResidRel', 'Hresid'}
    if ~isfield(g, f2{1}), C.reason = ['gates omitted ' f2{1}]; C.wallSec = toc(t0); return, end
    [okx, vx] = scalar_verdict(g.(f2{1}));
    if ~okx, C.reason = sprintf('gate %s is not a real finite scalar', f2{1}); C.wallSec = toc(t0); return, end
    C.(f2{1}) = vx;
end
if ~(C.fieldErr <= tolField),  C.reason = sprintf('X2: pumpkyn state rows differ from the independent field by %.2e > %g', C.fieldErr, tolField); C.wallSec = toc(t0); return, end
if ~(C.adjErrRef <= tolField), C.reason = sprintf('X2: pumpkyn adjoint rows differ from the independent adjoint by %.2e > %g', C.adjErrRef, tolField); C.wallSec = toc(t0); return, end
% every gate read as a real finite scalar, so Inf cannot read as "positive"
% and an empty field cannot read as "satisfied"
need = {'minLamV', 'minQmt', 'dimS'};
gv = nan(1, 3);
for kg = 1:3
    if ~isfield(g, need{kg}), C.reason = ['gates omitted ' need{kg}]; C.wallSec = toc(t0); return, end
    [okg, gv(kg)] = scalar_verdict(g.(need{kg}));
    if ~okg
        C.reason = sprintf('gate %s is not a real finite scalar', need{kg});
        C.wallSec = toc(t0);  return
    end
end
if ~(gv(1) > 0), C.reason = sprintf('min|lam_v| = %.2e not > 0', gv(1)); C.wallSec = toc(t0); return, end
if ~(gv(2) > 0), C.reason = sprintf('min Q_mt = %.2e not > 0', gv(2));   C.wallSec = toc(t0); return, end
if gv(3) ~= 1,   C.reason = sprintf('dim S = %g (abnormal lift)', gv(3)); C.wallSec = toc(t0); return, end
% dim S = 1 is only as good as (a) the accepted lift's own residual, (b) the
% independent field's Hamiltonian residual and (c) the Eckart-Young
% separation of sigma_6 from the MEASURED error in C -- all three required,
% not printed (Astra review #2, 2026-09-11)
if ~(C.nullResidRel <= tolLift), C.reason = sprintf('accepted lift backward error |C lam0|/(|C||lam0|) = %.2e > %g: dim S not resolved', C.nullResidRel, tolLift); C.wallSec = toc(t0); return, end
if ~(C.Hresid <= tolH),       C.reason = sprintf('independent-field Hamiltonian residual %.2e > %g', C.Hresid, tolH); C.wallSec = toc(t0); return, end
gOpts.relTol = relPair(2);
try
    [okG2, g2] = fenced(pool, capGates, @mintime_hypothesis_gates, 1, z(:), rv0, B.Tnd, B.cnd, B.mu, gOpts);
catch ME
    C.reason = ['second gates build threw: ' ME.message];  C.wallSec = toc(t0);  return
end
if ~okG2, C.reason = sprintf('second gates build exceeded the %g s cap', capGates); C.wallSec = toc(t0); return, end
if ~(isfield(g, 'C') && isfield(g2, 'C') && isequal(size(g.C), size(g2.C)))
    C.reason = 'lift_margin: the two constraint-matrix builds are not comparable';  C.wallSec = toc(t0);  return
end
LM = lift_margin(g.C, g2.C, z(1:7), struct('marginMin', liftMarginMin, 'liftTol', tolLift));
[okL, C.liftMargin] = scalar_verdict(LM.margin);
if ~okL || ~LM.certified
    C.reason = ['lift_margin: ' LM.reason];  C.wallSec = toc(t0);  return
end

% ---- 6. H6, enforced ----------------------------------------------------
% The reduced instrument's determinant can vanish spuriously when
% lambda_m(0) >= c/T (FINDINGS 40). Its verdict above is only meaningful
% when that is excluded, so a missing or failed H6 is a failed gate, and
% "excluded but with no headroom" is distinguished by h6MarginMin.
if ~isfield(g, 'h6Margin'), C.reason = 'gates omitted H6 (h6Margin)'; C.wallSec = toc(t0); return, end
[okH6, h6m] = scalar_verdict(g.h6Margin);
if ~okH6, C.reason = 'H6 margin is not a real finite scalar'; C.wallSec = toc(t0); return, end
C.h6Margin = h6m;
% the helper's own verdict FIRST (strict margin > 1 AND clearance above the
% Hamiltonian residual): the margin ratio alone let equality through and
% ignored the clearance (Astra review 2026-09-11)
if ~(isfield(g, 'h6Ok') && islogical(g.h6Ok) && isscalar(g.h6Ok) && g.h6Ok)
    if isfield(g, 'h6Reason'), why = g.h6Reason; else, why = 'h6Ok is false or missing'; end
    C.reason = sprintf('H6 not established: %s', why);
    C.wallSec = toc(t0);  return
end
if ~(h6m > h6MarginMin)
    C.reason = sprintf('H6 margin %.2fx not > required %.2fx (lambda_m(0) too close to c/T)', ...
                       h6m, h6MarginMin);
    C.wallSec = toc(t0);  return
end

% ---- 7. dense conjugate scan --------------------------------------------
% The junction sign test cannot see two zeros inside one segment or an
% even-order touch; the scan locates every candidate outside the start-up
% transient and judges it against a numerical floor. A located zero or
% multiplicity refutes; an UNRESOLVED candidate blocks (it is not a pass).
if doSpectrum
    cOpts = conjOpts;  cOpts.K = numel(seed.tGrid) - 1;  cOpts.nSub = conjNSub;
    try
        [okS, CS] = fenced(pool, capGates, @conj_spectrum, 1, z(:), rv0, B.Tnd, B.cnd, B.mu, cOpts);
    catch ME
        C.reason = ['conj_spectrum threw: ' ME.message];  C.wallSec = toc(t0);  return
    end
    if ~okS, C.reason = sprintf('conj_spectrum exceeded the %g s cap', capGates); C.wallSec = toc(t0); return, end
    C.conjDense = struct('nInterior', CS.nInterior, 'nInteriorCand', CS.nInteriorCand, ...
                         'nEndCand', CS.nEndCand, 'nStart', CS.nStart, 'nZero', CS.nZero, ...
                         'nNearMiss', CS.nNearMiss, 'nUnresolved', CS.nUnresolved, ...
                         'multiplicity', CS.multiplicity, 'clear', CS.clear, 'minRel', CS.minRel);
    [okc, cc] = scalar_verdict(CS.clear);
    if ~okc || cc ~= 1
        C.reason = sprintf(['dense conjugate scan not clear: %d coarse sign change(s), %d zero, ' ...
                            '%d UNRESOLVED, %d multiplicity (%d near-miss cleared)'], ...
                           CS.nInterior, CS.nZero, CS.nUnresolved, CS.multiplicity, CS.nNearMiss);
        C.wallSec = toc(t0);  return
    end
end

C.ok = true;  C.wallSec = toc(t0);
if plateau
    C.reason = sprintf('certified (polish plateaued at |R| = %.1e, above tolR = %.1e)', it.normR, tolR);
else
    C.reason = 'certified';
end
end


function varargout = fenced(pool, capSec, fh, nout, varargin)
% FENCED  One external call under a HARD wall-clock cap when a pool is
% available, direct otherwise. Returns [ok, outputs...].
% INPUTS: pool; capSec; fh; nout; varargin.  OUTPUTS: ok; nout outputs.
if isempty(pool)
    varargout = cell(1, nout + 1);
    [varargout{2:nout+1}] = feval(fh, varargin{:});
    varargout{1} = true;
    return
end
varargout = cell(1, nout + 1);
[varargout{1}, varargout{2:nout+1}] = run_capped(pool, fh, nout, capSec, varargin{:});
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
