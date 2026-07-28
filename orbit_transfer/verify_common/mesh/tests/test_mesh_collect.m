% TEST_MESH_COLLECT  Unit test for per-quantity collection and the study report.
%
% Hand-built ladders with known series -- no solves, no campaign dependencies.
% The properties under test are the ones that decide what the study is allowed
% to CLAIM, not merely what it computes:
%   - a clean O(h^2) physical series earns CONVERGED with p ~ 2 and is labelled
%     as supporting H2;
%   - an integer count is STABLE only if it never moves;
%   - diagnostics never earn CONVERGED, however clean their series;
%   - Richardson-forbidden quantities report no p at all, not a plausible one;
%   - two levels yield INSUFFICIENT rather than a fitted order;
%   - a failed level is DROPPED from the series rather than fitted.
root = fileparts(fileparts(mfilename('fullpath')));  cd(root);  addpath(root);

% ---- (a) clean 2nd-order mass series -> CONVERGED, p ~ 2, H2 -------------
h  = [1 1/2 1/4 1/8];
mf = 1377 + 3*h.^2;
L  = local_ladder([1 2 4 8], mf, [19 19 19 19]);
S  = mesh_collect(L);

assert(abs(S.q.mf.order.p - 2) < 1e-6, 'a: p = 2 expected for m_f, got %.6f', S.q.mf.order.p);
assert(abs(S.q.mf.order.rich - 1377) < 1e-6, 'a: Richardson should recover 1377, got %.9f', S.q.mf.order.rich);
assert(strcmp(S.q.mf.class, 'physical'), 'a: m_f must be classed physical');
assert(S.q.mf.richardsonOK, 'a: Richardson is permitted on m_f');
T = mesh_report(S, 'TEST_a', '');
vm = T.verdict{strcmp(T.quantity, 'final mass')};
assert(strcmp(vm, 'CONVERGED'), 'a: expected CONVERGED for m_f, got %s', vm);
hm = T.hypothesis{strcmp(T.quantity, 'final mass')};
assert(strcmp(hm, 'H2'), 'a: p ~ 2 must be labelled H2, got %s', hm);

% ---- (b) a 1st-order series is labelled H1 ------------------------------
L  = local_ladder([1 2 4 8], 1377 + 3*h, [19 19 19 19]);
S  = mesh_collect(L);
assert(abs(S.q.mf.order.p - 1) < 1e-6, 'b: p = 1 expected, got %.6f', S.q.mf.order.p);
T  = mesh_report(S, 'TEST_b', '');
hm = T.hypothesis{strcmp(T.quantity, 'final mass')};
assert(strcmp(hm, 'H1'), 'b: p ~ 1 must be labelled H1, got %s', hm);

% ---- (c) switch count: STABLE only when it never moves ------------------
S = mesh_collect(local_ladder([1 2 4 8], mf, [19 19 19 19]));
assert(S.switchStability.stable, 'c: constant counts must read STABLE');
S = mesh_collect(local_ladder([1 2 4 8], mf, [19 20 20 20]));
assert(~S.switchStability.stable, 'c: a count that moves must NOT read STABLE');
T = mesh_report(S, 'TEST_c', '');
vs = T.verdict{strcmp(T.quantity, 'switch count')};
assert(strcmp(vs, 'NOT-STABLE'), 'c: expected NOT-STABLE, got %s', vs);

% ---- (d) diagnostics NEVER earn CONVERGED -------------------------------
% Give kktStatInf a textbook 2nd-order series. A physical quantity with this
% series is CONVERGED; a diagnostic with the same series must not be.
L = local_ladder([1 2 4 8], mf, [19 19 19 19]);
for k = 1:4, L(k).rep.kktStatInf = 1e-12 + 3e-9*h(k)^2; end
S = mesh_collect(L);
assert(strcmp(S.q.kktStatInf.class, 'diagnostic'), 'd: kktStatInf must be classed diagnostic');
assert(isnan(S.q.kktStatInf.order.p), 'd: no p may be reported for a diagnostic');
T = mesh_report(S, 'TEST_d', '');
vk = T.verdict{strcmp(T.quantity, 'KKT stationarity')};
assert(~strcmp(vk, 'CONVERGED'), 'd: a diagnostic must never read CONVERGED, got %s', vk);

% ---- (e) an O(h) diagnostic is called out as a representation artifact ---
L = local_ladder([1 2 4 8], mf, [19 19 19 19]);
for k = 1:4, L(k).rep.lamMassEndMapped = 1.2e-3 * h(k); end   % halves per doubling
S = mesh_collect(L);
T = mesh_report(S, 'TEST_e', '');
vl = T.verdict{strcmp(T.quantity, 'mapped transversality')};
assert(strcmp(vl, 'O(h)-CONSISTENT'), 'e: expected O(h)-CONSISTENT, got %s', vl);

% ---- (f) two levels -> INSUFFICIENT, never a fitted order ---------------
S = mesh_collect(local_ladder([1 2], [1377.1 1377.2], [19 20]));
assert(isnan(S.q.mf.order.p), 'f: two levels cannot support an order');
T = mesh_report(S, 'TEST_f', '');
vm = T.verdict{strcmp(T.quantity, 'final mass')};
assert(strcmp(vm, 'INSUFFICIENT'), 'f: expected INSUFFICIENT, got %s', vm);

% ---- (g) a FAILED level is dropped, not fitted --------------------------
L = local_ladder([1 2 4 8], mf, [19 19 19 19]);
L(4).ok = false;  L(4).why = 'synthetic failure';
w = warning('off', 'mesh_collect:droppedLevels');
S = mesh_collect(L);
warning(w);
assert(S.nLevels == 3, 'g: the failed level must be dropped, got %d levels', S.nLevels);
assert(isequal(S.factors, [1 2 4]), 'g: factors should be [1 2 4], got [%s]', num2str(S.factors));

% ---- (h) the per-switch window is used, and is reported ----------------
% Two levels whose switches sit where the local step differs by 20x. A global
% window built from the median step would mis-handle one end; the per-switch
% window must pair both.
S = mesh_collect(local_ladder([1 2], mf(1:2), [2 2]));
assert(all(isfinite(S.switchStability.tolMedian(1))), 'h: a window must be recorded');
assert(S.switchStability.matched(1) == 2, ...
    'h: both switches should pair, got %d', S.switchStability.matched(1));

fprintf('test_mesh_collect: ALL PASS\n');

% ===========================================================================
function L = local_ladder(factors, mfVals, nSw)
% LOCAL_LADDER  Build a synthetic ladder with prescribed series.
%
% Each level carries a small but INTERNALLY CONSISTENT solution: a physical
% time row that is deliberately NONUNIFORM (so the per-switch window is
% exercised), a throttle with nSw(k) transitions, and a switching function
% that changes sign across each of them so the sub-grid extractor engages.
%
% INPUTS:  factors [1xL]; mfVals [1xL] final mass per level; nSw [1xL] switches
% OUTPUTS: L - ladder struct array [1xL]
L = struct('factor',{},'N',{},'sigma',{},'out',{},'rep',{},'wall',{},'ok',{},'why',{});
for k = 1:numel(factors)
    N   = 40 * factors(k);
    sg  = linspace(0, 1, N+1).';
    % nonuniform physical time: quadratic in sigma -> local step varies ~20x
    t   = 30 * (0.05*sg + 0.95*sg.^2).';
    thr = zeros(1, N+1);
    Sh  = ones(1, N+1);
    % place nSw(k) transitions at evenly spaced interior node pairs
    if nSw(k) > 0
        idx = round(linspace(0.15*N, 0.85*N, nSw(k)));
        state = false;
        for q = 1:numel(idx)
            state = ~state;
            thr(idx(q)+1:end) = double(state);
            Sh(idx(q)+1:end)  = -Sh(idx(q));      % sign change at each switch
        end
    end
    X = ones(7, N+1);
    X(6,:) = linspace(1, mfVals(k)/1500, N+1);
    X(7,:) = t;
    U = [zeros(3, N+1); thr];
    U(1,:) = 1;                                    % unit beta
    out = struct('X', X, 'U', U, 'dL', 46, 'm_f_kg', mfVals(k), ...
        'dV_kms', 1.6765, 'termErr', 1e-12, 'maxDefect', 1e-15, ...
        'switches', nSw(k), 'success', true, 'ipoptStatus', 'Solve_Succeeded');
    rep = struct('kktStatInf', 1e-14, 'lamMassEndMapped', 1e-18, ...
        'lamTimeCoV', 1e-9, 'sdotMinRel', 0.5, 'signPct', 100, ...
        'dirTanMax', 0, 'Sdeweighted', Sh, 'nSwitches', nSw(k));
    L(end+1) = struct('factor', factors(k), 'N', N, 'sigma', sg, 'out', out, ...
        'rep', rep, 'wall', 1, 'ok', true, 'why', ''); %#ok<AGROW>
end
end
