% TEST_MS_JACOBIAN  Task-4 gate: complex-step Jacobian vs finite differences.
%
% Gate design (v3, review-prescribed): no single FD step h certifies all
% columns here -- stiff perigee columns are roundoff-dominated at small h,
% truncation-dominated at large h, and the h/(h/2) self-consistency estimate
% can deflate when the two FD errors correlate. So Gate 1 sweeps h
% largest-first and passes a column if ANY h satisfies the self-consistency
% formula: each column certifies in its own FD-convergent regime. Gate 2
% adds (a) a CS h-independence check (hCS 1e-16 vs 1e-24 agree to 1e-10;
% measured ~1e-13) and (b) assertion B2: the assembled J block must match
% the independently rebuilt CS column -- catches assembly/indexing bugs
% independent of FD conditioning. Structure gate unchanged.
setup_paths;
rng(7);
epsS = 0.3;
prob = ms_problem(1.03, epsS);
ref  = run_gto_tulip_indirect(false);
lam0 = ref.zSol(1:7);
y0   = [prob.rv0; prob.m0; lam0];
sol  = ode113(@(t,y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsS), [0 prob.tf], y0, prob.odeOpts);
M    = 4;
prob.tJ = arc_boundaries_tau(sol.x, sol.y(1:3,:), M, prob.muStar);
yJ   = deval(sol, prob.tJ);
Z    = ms_pack(lam0, yJ(:, 2:M));

[~, J] = ms_residual(Z, prob);
n = numel(Z);
colSel = unique([1 3 7, 8 15 22, n-13, n, randi(n, 1, 4)]);   % lam0 + joint cols
hSweep = [1e-4 1e-5 1e-6 1e-7];    % largest-first: early break on pass
failMsg = '';

% Gate 1 (any-h): column passes if ANY h satisfies the FD self-consistency
% tolerance max(1e-6*denom, 3*selfErr)
for cIdx = colSel
    scale  = max(1, abs(Z(cIdx)));
    passH  = NaN;
    for hFD = hSweep
        colFD   = fd_col(Z, prob, cIdx, hFD*scale);
        colFDh2 = fd_col(Z, prob, cIdx, hFD*scale/2);
        denom   = max(1, max(abs(colFD)));
        selfErr = max(abs(colFD - colFDh2));
        tolC    = max(1e-6*denom, 3*selfErr);
        errC    = max(abs(full(J(:, cIdx)) - colFD));
        okH     = errC <= tolC;
        fprintf('col %3d  h=%.0e  |J-FD| %.3e   tol %.3e   (selfErr %.3e)  %s\n', ...
                cIdx, hFD, errC, tolC, selfErr, ternStr(okH, 'pass', 'fail'));
        if okH, passH = hFD; break; end
    end
    if isnan(passH)
        failMsg = sprintf('%s col %d: no h in sweep passes FD gate;', failMsg, cIdx);
    else
        fprintf('col %3d  Gate-1 PASS at h=%.0e\n', cIdx, passH);
    end
end

% Gate 2: CS h-independence + assertion B2 (assembled J vs rebuilt column)
for cIdx = colSel
    [kArc, idx] = col_owner(cIdx);
    c16 = cs_col(yJ(:, kArc), idx, prob.tJ(kArc), prob.tJ(kArc+1), prob, 1e-16);
    c24 = cs_col(yJ(:, kArc), idx, prob.tJ(kArc), prob.tJ(kArc+1), prob, 1e-24);
    relCS = max(abs(c16 - c24))/max(1, max(abs(c16)));
    if kArc <= M-1
        b2 = max(abs(full(J(14*(kArc-1)+(1:14), cIdx)) - c16));
    else
        b2 = max(abs(full(J(14*(M-1)+(1:7), cIdx)) - c16([1:6 14])));
    end
    b2rel = b2/max(1, max(abs(c16)));
    fprintf('col %3d  CS h-independence %.3e   B2 assembled-vs-rebuilt %.3e\n', ...
            cIdx, relCS, b2rel);
    if relCS > 1e-10
        failMsg = sprintf('%s col %d: CS h-dependence %.3e > 1e-10;', ...
                          failMsg, cIdx, relCS);
    end
    if b2rel > 1e-9
        failMsg = sprintf('%s col %d: B2 assembled-vs-rebuilt %.3e > 1e-9;', ...
                          failMsg, cIdx, b2rel);
    end
end

% Gate 3 (unchanged): rows of continuity block k depend only on y_k, y_{k+1}
bw = 0;
for k = 1:M-1
    rowsK = 14*(k-1)+(1:14);
    colsAllowed = false(1, n);
    if k == 1, colsAllowed(1:7) = true; else, colsAllowed(7+14*(k-2)+(1:14)) = true; end
    colsAllowed(7+14*(k-1)+(1:14)) = true;
    bw = max(bw, nnz(any(J(rowsK, ~colsAllowed), 1)));
end
fprintf('off-structure cols: %d\n', bw);
if bw ~= 0
    failMsg = sprintf('%s off-structure cols %d ~= 0;', failMsg, bw);
end

if isempty(failMsg)
    fprintf('PASS test_ms_jacobian\n');
else
    error('FAIL test_ms_jacobian:%s', failMsg);   % nonzero exit under -batch
end

% -------------------------------------------------------------------------
function colFD = fd_col(Z, prob, cIdx, h)
% FD_COL  Central-FD column of the MS residual Jacobian.
%
% INPUTS:
%   Z    - unknown vector [(14M-7)x1]
%   prob - problem struct from MS_PROBLEM with tJ set [1x(M+1)]
%   cIdx - column (unknown) index to perturb [scalar]
%   h    - absolute FD step [scalar]
%
% OUTPUTS:
%   colFD - central-difference Jacobian column [(14M-7)x1]
Zp = Z;  Zp(cIdx) = Zp(cIdx) + h;
Zm = Z;  Zm(cIdx) = Zm(cIdx) - h;
colFD = (ms_residual(Zp, prob) - ms_residual(Zm, prob))/(2*h);
end

function [kArc, idx] = col_owner(cIdx)
% COL_OWNER  Map unknown-vector column to (owning arc, state component).
%
% Z = [lam0(7); y_2(14); ...; y_M(14)]: cols 1..7 perturb arc-1 costates
% (components 8..14); col 7+14*(k-2)+d perturbs component d of arc k >= 2.
%
% INPUTS:
%   cIdx - unknown-vector column index [scalar]
%
% OUTPUTS:
%   kArc - owning arc index [scalar]
%   idx  - perturbed component of that arc's initial state, 1..14 [scalar]
if cIdx <= 7
    kArc = 1;  idx = cIdx + 7;
else
    kArc = floor((cIdx - 8)/14) + 2;
    idx  = mod(cIdx - 8, 14) + 1;
end
end

function col = cs_col(y0Arc, idx, t0, t1, prob, hCS)
% CS_COL  One complex-step STM column of a single arc (MS_JACOBIAN_CS scheme).
%
% INPUTS:
%   y0Arc - arc initial augmented state [14x1]
%   idx   - component to perturb, 1..14 [scalar]
%   t0    - arc start time (ND) [scalar]
%   t1    - arc end time (ND) [scalar]
%   prob  - problem struct (Tmax, c, muStar, epsSmooth, odeOpts)
%   hCS   - complex-step size [scalar]
%
% OUTPUTS:
%   col - STM column d yEnd / d y0Arc(idx) [14x1]
scale   = max(1, abs(y0Arc(idx)));
yp      = complex(y0Arc);
yp(idx) = yp(idx) + 1i*hCS*scale;
[~, Yc] = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
                 prob.muStar, prob.epsSmooth), [t0 t1], yp, prob.odeOpts);
col = imag(Yc(end, :).')./(hCS*scale);
end

function s = ternStr(cond, a, b)
% TERNSTR  Ternary string select.
%
% INPUTS:
%   cond - condition [logical scalar]
%   a    - string returned when cond is true
%   b    - string returned when cond is false
%
% OUTPUTS:
%   s - a or b
if cond, s = a; else, s = b; end
end
