function ok = test_conj_resolve()
%% Purpose:
%
%   ADVERSARIAL SYNTHETIC tests of conj_resolve, the candidate-resolution
%   logic behind the dense conjugate scan, driven by matrix-valued
%   functions whose rank loss is known exactly (Astra review #3,
%   2026-09-11: "extract the resolution logic so it can accept a synthetic
%   matrix-valued function", then the list below). Every case asserts the
%   GATE, .clear, not merely that some candidate exists:
%
%     1. a simple transverse zero (sign change) at every grid phase,
%        including exactly on a coarse node and within h/32 of one;
%     2. a corank-two zero WITHOUT a sign change (det ~ (t - t*)^2);
%     3. a quadratic touch of sigma_n without a vanishing column;
%     4. multiple wells in one cluster: two near-misses around a zero;
%     5. an endpoint-only candidate: the zero AT t_f, and just before it;
%     6. a root merged with the start-up transient;
%     7. a vanishing column (no dip in the normalised spectrum at all);
%     8. a positive near-miss well above the floor -> CLEARED;
%     9. a minimum inside the floor band -> UNRESOLVED;
%    10. a scan that never leaves rank deficiency -> not testable, not clear.
%
%   The matrices mix all columns through fixed random orthogonal factors so
%   the column norms stay O(1) except where a vanishing column is intended.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(7);
n = 6;  N = 48;  tf = 4;  tG = (1:N)*tf/N;  h = tf/N;
[Q, ~] = qr(randn(n));  [W, ~] = qr(randn(n));
% base spectrum: sigma = [3 2.5 2 1.5 1.2 s(t)], transient s ~ tanh(t/0.3)
base = @(t) diag([3 2.5 2 1.5 1.2 1]);
startup = @(t) min(1, (t/1.2).^3);                 % structural growth from 0 (the real block is ~t^2..t^3)
mk = @(sfun) @(t) Q * (base(t) .* diag([1 1 1 1 1 sfun(t)])) * W';   % sigma_6(t) = |sfun(t)| * startup-free
% NOTE: base has sigma_6 = 1 before sfun; the last entry becomes sfun(t)

% ---- 1. simple transverse zero at every grid phase ------------------------
tzs = [2.5, 2.5 + h/64, 2.5 + h/32, 2.5 + h/3, 2.5 + h/2, 2.5 + 0.9*h];
for tz = tzs
    R = conj_resolve(mk(@(t) startup(t) .* (t - tz)), tG);
    ok = chk(ok, ~R.clear && R.nZero >= 1 && R.nZeroSign >= 1, ...
             sprintf('transverse zero at t* = 2.5 + %.3fh: NOT clear, %d zero(s) by sign (%s)', (tz-2.5)/h, R.nZeroSign, R.reason));
end

% ---- 2. corank-two zero, no sign change, ON a node and OFF it ------------
for tz = [2.5, 2.5 + h/3, 2.5 + h/32]
    sq2 = @(t) startup(t) .* abs(t - tz);
    M2 = @(t) Q * (diag([3 2.5 2 1.5 sq2(t) sq2(t)])) * W';   % sigma_5 = sigma_6 -> 0 together, det >= 0
    R = conj_resolve(M2, tG);
    ok = chk(ok, ~R.clear && R.nZero >= 1 && R.nZeroFloor >= 1, ...
             sprintf('corank-two zero without a sign change at 2.5 + %.3fh: NOT clear, %d floor-level, corank %d', ...
                     (tz-2.5)/h, R.nZeroFloor, max([R.candidates.nSmall])));
    ok = chk(ok, R.multiplicity >= 1, sprintf('and reported as multiplicity (%d)', R.multiplicity));
end

% ---- 3. quadratic touch (sigma_6 = (t-t*)^2, columns stay O(1)), off-node --
for tz = [2.5, 2.5 + h/3]
    R = conj_resolve(mk(@(t) startup(t) .* (t - tz).^2), tG);
    ok = chk(ok, ~R.clear && R.nZero >= 1, sprintf('quadratic touch at 2.5 + %.3fh: NOT clear (%s)', (tz-2.5)/h, R.reason));
end

% ---- 4. multiple wells in one cluster -------------------------------------
wells = @(t) startup(t) .* (5e-4 + 0.05*abs(t - 2.3)) .* (5e-4 + 0.05*abs(t - 2.7)) .* abs(t - 2.5) * 40;
R = conj_resolve(mk(wells), tG);
ok = chk(ok, ~R.clear && R.nZero >= 1, sprintf('two near-miss wells around a zero in one cluster: NOT clear (%s)', R.reason));
wellsOnly = @(t) startup(t) .* (2e-4 + 0.5*abs(t - 2.3)) .* (2e-4 + 0.5*abs(t - 2.7)) * 25;
R = conj_resolve(mk(wellsOnly), tG);
ok = chk(ok, R.clear && R.nNearMiss >= 1 && R.nZero == 0, ...
         sprintf('the same two wells with the zero removed: CLEAR as near-miss (%s)', R.reason));

% ---- 5. endpoint-only -----------------------------------------------------------
R = conj_resolve(mk(@(t) startup(t) .* (t - tf)), tG);
ok = chk(ok, ~R.clear && R.nEndCand >= 1, sprintf('zero AT t_f: NOT clear, endpoint candidate (%s)', R.reason));
R = conj_resolve(mk(@(t) startup(t) .* (t - (tf - h/5))), tG);
ok = chk(ok, ~R.clear && R.nZero >= 1, sprintf('zero one fifth of a step before t_f: NOT clear (%s)', R.reason));

% ---- 6. a root merged with the start-up transient -------------------------
%    t* = 0.2 lies between coarse samples 2 and 3, inside the region where
%    the transient is still below tolCollapse: at coarse resolution the
%    record is monotone, so this is exactly the absorbed case
R = conj_resolve(mk(@(t) startup(t) .* abs(t - 0.2)), tG);
ok = chk(ok, ~R.clear && R.nZero >= 1, sprintf('a zero inside the start transient (t* = 0.2): NOT clear (%s)', R.reason));
ok = chk(ok, R.tUncovered < 0.2, sprintf('and the uncovered prefix stops before it, at t = %.4f', R.tUncovered));

% ---- 7. a vanishing column (no dip in the normalised spectrum) ----------------
Mv = @(t) Q * diag([(t - 2.5)^2 1 1 1 1 1]);       % column 1 of Q*diag vanishes at t*
R = conj_resolve(Mv, tG);
ok = chk(ok, ~R.clear && (R.nUnresolved >= 1 || R.nZero >= 1), ...
         sprintf('vanishing column: NOT clear (%s); min col ratio %.1e', R.reason, min(R.colRatio)));

% ---- 8. a positive near-miss well above the floor -> cleared -------------------
R = conj_resolve(mk(@(t) startup(t) .* (1e-4 + abs(t - 2.5))), tG);
ok = chk(ok, R.clear && R.nNearMiss == 1 && R.nZero == 0 && R.nUnresolved == 0, ...
         sprintf('near-miss with minimum 1e-4 (1000x the floor): CLEAR, located at t/t_f = %.4f', R.candidates(end).tMinOverTf));
ok = chk(ok, abs(R.candidates(end).tMinOverTf*tf - 2.5) < 1e-6, 'and located to 1e-6');

% ---- 9. a minimum inside the floor band -> unresolved -----------------------
R = conj_resolve(mk(@(t) startup(t) .* (3e-6 + abs(t - 2.5))), tG);
ok = chk(ok, ~R.clear && R.nUnresolved == 1, sprintf('minimum 3e-6 (30x the floor, under the 100x clearance): UNRESOLVED (%s)', R.reason));

% ---- 10. never full rank ------------------------------------------------------------
R = conj_resolve(@(t) Q * diag([3 2.5 2 1.5 1.2 0]) * W', tG);
ok = chk(ok, ~R.testable && ~R.clear, sprintf('rank-deficient throughout: not testable, not clear (%s)', R.reason));

% ---- and a clean scan is clear, with its uncovered start reported ----------------
R = conj_resolve(mk(@(t) startup(t)), tG);
ok = chk(ok, R.clear && R.nStart == 1 && R.tUncovered > 0 && R.tUncovered < 0.6, ...
         sprintf('a clean scan is CLEAR with the start transient reported uncovered to t = %.3f', R.tUncovered));

if ok, fprintf('TEST_CONJ_RESOLVE: ALL PASS\n'); else, fprintf('TEST_CONJ_RESOLVE: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
