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
%         .lamMassEndRel    free-final-mass transversality residual (relative)
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
if ~isempty(man.thrRow)
    thr = U(man.thrRow,:);  burn = thr > 0.5;
    lamNB = lamAll;  lamNB([rowsOf('thrLo'), rowsOf('thrHi')]) = 0;
    gNB = gf + s*(A.'*lamNB);           % Lagrangian gradient sans throttle-bound duals
    Sd = zeros(1,N1);
    for k = 1:N1, Sd(k) = gNB(uix(man.thrRow,k)); end
    rep.Sd = Sd;
    rep.signPct = 100*mean((Sd < 0) == burn);      % S<0 <=> full thrust
    checks{end+1} = 'signLaw';
else
    burn = true(1,N1);  rep.Sd = [];  rep.signPct = NaN;   % min-time: all-burn
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
if ~isempty(man.massRow) && man.massFreeAtTf
    lm = rep.lam(man.massRow,:);
    rep.lamMassEndRel = abs(lm(end))/max(abs(lm));
    checks{end+1} = 'transversality';
else
    rep.lamMassEndRel = NaN;
end

% --- (7) singular arcs + regular switching (Sdot != 0) -----------------------
rep.singularArcNodes = 0;  rep.sdotMinRel = NaN;  rep.nSwitches = 0;
if ~isempty(man.thrRow)
    coastS = abs(rep.Sd(~burn));
    scaleS = median(coastS(coastS>0));  if isempty(scaleS)||isnan(scaleS), scaleS = max(abs(rep.Sd)); end
    nearZ = abs(rep.Sd) < max(1e-6*scaleS, 1e-14);
    runL = 0;
    for k = 1:N1
        if nearZ(k), runL = runL+1; else, runL = 0; end
        if runL >= 3, rep.singularArcNodes = rep.singularArcNodes + 1; end
    end
    swI = find(diff(burn) ~= 0);  rep.nSwitches = numel(swI);
    if ~isempty(swI)
        sdot = abs(rep.Sd(min(swI+1,N1)) - rep.Sd(max(swI,1))) / max(scaleS,1e-30);
        rep.sdotMinRel = min(sdot);
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

% --- advisory verdict (REPORT-ONLY burn-in) ----------------------------------
okSign  = isempty(man.thrRow) || rep.signPct >= tolSign;
okTrans = isnan(rep.lamMassEndRel) || rep.lamMassEndRel <= tolTrans;
okSdot  = isnan(rep.sdotMinRel) || rep.sdotMinRel > sdotMin;
rep.pass = rep.kktStatInf <= tolStat && rep.dirTanMax <= tolStat && ...
           okSign && okTrans && rep.singularArcNodes == 0 && okSdot;
end

function v = fcdef(o, f, dflt)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = dflt; end
end
