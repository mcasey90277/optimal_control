% TEST_SMS_JACOBIAN  Gate B: 16-dim complex-step Jacobian vs finite differences.
%
% Port of the v3 gate design (test_ms_jacobian) to the Sundman-domain
% system: no single FD step h certifies all columns -- stiff columns are
% roundoff-dominated at small h, truncation-dominated at large h, and the
% h/(h/2) self-consistency estimate can deflate when the two FD errors
% correlate. Gate 1 sweeps h largest-first and passes a column if ANY h
% satisfies the self-consistency formula. Gate 2 adds (a) a CS
% h-independence check (hCS 1e-16 vs 1e-24 agree to 1e-10) and (b)
% assertion B2: the assembled J block must match an independently rebuilt
% CS column -- catches assembly/indexing bugs independent of FD
% conditioning. Gate 3: block-bidiagonal structure. Dims 16: arc-1 unknowns
% are the 8 initial costates; terminal rows are Phi_M([1:6 15 8],:).
setup_paths;
rng(7);
epsS = 0.3;
M    = 4;
[Z, prob] = sms_seed_mintime(1.03, M, epsS);
[~, yJ]   = sms_unpack(Z, prob);

[~, J] = sms_residual(Z, prob);
n = numel(Z);
colSel = unique([1 4 8, 9 17 24, n-15, n, randi(n, 1, 4)]);  % lam0 + joint cols
hSweep = [1e-3 1e-4 1e-5 1e-6 1e-7];   % largest-first: early break on pass
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
    c16 = cs_col(yJ(:, kArc), idx, prob.sJ(kArc), prob.sJ(kArc+1), prob, 1e-16);
    c24 = cs_col(yJ(:, kArc), idx, prob.sJ(kArc), prob.sJ(kArc+1), prob, 1e-24);
    relCS = max(abs(c16 - c24))/max(1, max(abs(c16)));
    if kArc <= M-1
        b2 = max(abs(full(J(16*(kArc-1)+(1:16), cIdx)) - c16));
    else
        b2 = max(abs(full(J(16*(M-1)+(1:8), cIdx)) - c16([1:6 15 8])));
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

% Gate 3: rows of continuity block k depend only on Y_k, Y_{k+1}
bw = 0;
for k = 1:M-1
    rowsK = 16*(k-1)+(1:16);
    colsAllowed = false(1, n);
    if k == 1, colsAllowed(1:8) = true; else, colsAllowed(8+16*(k-2)+(1:16)) = true; end
    colsAllowed(8+16*(k-1)+(1:16)) = true;
    bw = max(bw, nnz(any(J(rowsK, ~colsAllowed), 1)));
end
fprintf('off-structure cols: %d\n', bw);
if bw ~= 0
    failMsg = sprintf('%s off-structure cols %d ~= 0;', failMsg, bw);
end

if isempty(failMsg)
    fprintf('PASS test_sms_jacobian\n');
else
    error('FAIL test_sms_jacobian:%s', failMsg);   % nonzero exit under -batch
end

% -------------------------------------------------------------------------
function colFD = fd_col(Z, prob, cIdx, h)
% FD_COL  Central-FD column of the Sundman-MS residual Jacobian.
%
% INPUTS:
%   Z    - unknown vector [(16M-8)x1]
%   prob - problem struct from SMS_PROBLEM with sJ set [1x(M+1)]
%   cIdx - column (unknown) index to perturb [scalar]
%   h    - absolute FD step [scalar]
%
% OUTPUTS:
%   colFD - central-difference Jacobian column [(16M-8)x1]
Zp = Z;  Zp(cIdx) = Zp(cIdx) + h;
Zm = Z;  Zm(cIdx) = Zm(cIdx) - h;
colFD = (sms_residual(Zp, prob) - sms_residual(Zm, prob))/(2*h);
end

function [kArc, idx] = col_owner(cIdx)
% COL_OWNER  Map unknown-vector column to (owning arc, state component).
%
% Z = [lam0(8); Y_2(16); ...; Y_M(16)]: cols 1..8 perturb arc-1 costates
% (components 9..16); col 8+16*(k-2)+d perturbs component d of arc k >= 2.
%
% INPUTS:
%   cIdx - unknown-vector column index [scalar]
%
% OUTPUTS:
%   kArc - owning arc index [scalar]
%   idx  - perturbed component of that arc's initial state, 1..16 [scalar]
if cIdx <= 8
    kArc = 1;  idx = cIdx + 8;
else
    kArc = floor((cIdx - 9)/16) + 2;
    idx  = mod(cIdx - 9, 16) + 1;
end
end

function col = cs_col(y0Arc, idx, s0, s1, prob, hCS)
% CS_COL  One complex-step STM column of a single arc (SMS_JACOBIAN_CS scheme).
%
% INPUTS:
%   y0Arc - arc initial augmented state [16x1]
%   idx   - component to perturb, 1..16 [scalar]
%   s0    - arc start sigma [scalar]
%   s1    - arc end sigma [scalar]
%   prob  - problem struct (Tmax, c, muStar, epsSmooth, pSund, odeOpts)
%   hCS   - complex-step size [scalar]
%
% OUTPUTS:
%   col - STM column d yEnd / d y0Arc(idx) [16x1]
scale   = max(1, abs(y0Arc(idx)));
yp      = complex(y0Arc);
yp(idx) = yp(idx) + 1i*hCS*scale;
[~, Yc] = ode113(@(s, y) sms_eom(s, y, prob.Tmax, prob.c, prob.muStar, ...
                 prob.epsSmooth, prob.pSund), [s0 s1], yp, prob.odeOpts);
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
