function ok = test_conj_fixedtf()
%% Purpose:
%
%   Validates the FIXED-FINAL-TIME conjugate-point spec for
%   ms_conjugate_test (freeTime = false, quotientDir = [], rows = the
%   components vanishing under the terminal conditions, cols = all initial
%   costates) against the textbook analytic case:
%
%     minimize Int_0^T (u^2 - x^2)/2 dt,  xdot = u,  x(0) = x(T) = 0
%
%   PMP: H = u^2/2 - x^2/2 + lam*u, u* = -lam; the variational flow is a
%   rotation, Phi_xl(t) = -sin(t) -- the FIRST CONJUGATE POINT IS AT
%   t = pi EXACTLY. So:
%     T = 2.0  (< pi): no interior sign change  -> pass
%     T = 4.0  (> pi): one interior sign change -> fail, crossing near pi
%     T = 6.6  (> 2pi): two interior sign changes
%
%   Segment STMs are built by chaining exact rotations (what ms_bvp's
%   keepSTMs would return for this field), junctions on a uniform grid.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All cases passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/02/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
spec = struct('stateRows', 1, 'costateCols', 2, ...
              'quotientDir', [], 'freeTime', false);
K = 48;

cases = {2.0, 0, false;  4.0, 1, false;  6.6, 2, false};
for kc = 1:size(cases, 1)
    T = cases{kc,1};
    info = lqInfo(T, K);
    out = ms_conjugate_test(info, spec);
    okC = out.nCrossings == cases{kc,2} && out.atFinal == cases{kc,3};
    ok = chk(ok, okC, sprintf('LQ T=%.1f: %d interior crossings (want %d), atFinal=%d', ...
             T, out.nCrossings, cases{kc,2}, out.atFinal));
    if T > pi && out.nCrossings >= 1
        % locate the first crossing: must bracket pi at junction resolution
        s = sign(out.detScaled);
        kf = find(diff(s) ~= 0, 1);
        tlo = out.t(kf);  thi = out.t(kf+1);
        ok = chk(ok, tlo < pi && pi < thi, ...
                 sprintf('LQ T=%.1f: first crossing brackets pi (%.3f..%.3f)', ...
                         T, tlo, thi));
    end
end

%% Review fixes 2026-09-05 (Astra/sol findings, minfuel_code_review):
% (a) the LAST segment must be monitored: put the pi crossing INSIDE the
%     final segment (K = 4, T = 3.5: t_K = 2.625 < pi < t_f = 3.5). The
%     old instrument stopped at t_K and passed this extremal.
info = lqInfo(3.5, 4);
out  = ms_conjugate_test(info, spec);
ok = chk(ok, out.nCrossings == 1 && ~out.pass, ...
         sprintf('LQ T=3.5 K=4: crossing in the FINAL segment counted (n=%d, pass=%d)', ...
                 out.nCrossings, out.pass));
ok = chk(ok, numel(out.t) == 4 && abs(out.t(end) - 3.5) < 1e-12, ...
         sprintf('LQ T=3.5 K=4: samples run through t_f (last sample %.3f)', out.t(end)));

% (b) an INITIAL COAST makes the state block structurally zero (no costate
%     dependence until the first burn). Those samples are not focal points:
%     skip them, report where full rank is first attained, and PASS when the
%     rest of the arc has no sign change. Fixture: 2 coast segments (state
%     block exactly zero) then the LQ rotation flow with T = 2 (< pi).
info = lqInfo(2.0, 8);
Z = eye(2);                                       % frozen subsystem: Phi_xl == 0, still a fundamental matrix
info.PHI(1:2) = {Z};                              % Phi_xl == 0 through junctions 2,3
out = ms_conjugate_test(info, spec);
ok = chk(ok, out.pass && out.nCrossings == 0, ...
         sprintf('initial coast: structural zeros skipped, no refutation (n=%d)', out.nCrossings));
ok = chk(ok, isfield(out, 'firstFullRank') && out.firstFullRank == 3, ...
         sprintf('initial coast: firstFullRank = %d (want 3)', ...
                 ternary(isfield(out,'firstFullRank'), out.firstFullRank, -1)));

% (c) provenance: the spec actually used is echoed in the output so a verdict
%     record can prove which block it monitored (P0.3 traceability).
ok = chk(ok, isfield(out, 'stateRows') && isequal(out.stateRows, spec.stateRows) && ...
             isfield(out, 'costateCols') && isequal(out.costateCols, spec.costateCols), ...
         'spec echo: out.stateRows / out.costateCols present and equal to the spec');

%% Astra review #2 (2026-09-06):
% (d) NOTHING TESTABLE must not be a PASS: state row never sees the costate.
info = lqInfo(2.0, 6);  info.PHI(:) = {eye(2)};
out = ms_conjugate_test(info, spec);
okD = isfield(out,'tested') && isfield(out,'verdict') && ~out.tested && ~out.pass && strcmp(out.verdict, 'UNDETERMINED');
ok = chk(ok, okD, sprintf('no testable sample -> UNDETERMINED, pass=false (tested=%d pass=%d)', ...
                 ternary(isfield(out,'tested'), double(out.tested), -1), out.pass));

% (e) raw-determinant UNDERFLOW must not create fictitious zero roots: TWO
%     independent LQ copies (2x2 state/costate block), each segment STM
%     scaled so the block's ENTRIES stay representable (~1e-160) while its
%     determinant (~1e-330) underflows to exactly 0 in double. T = 2 < pi:
%     no root expected; the old raw-det path reported a zero root here.
K = 4;  tg = linspace(0, 2, K+1);  h = tg(2);
R = [cos(h), -sin(h); sin(h), cos(h)];  sc = 10^(-165/K);
P4 = sc*blkdiag(R, R);                                % states 1,3; costates 2,4
info = struct('PHI', {repmat({P4}, 1, K)}, 'Y', zeros(4, K), 'tGrid', tg);
out = ms_conjugate_test(info, struct('stateRows', [1 3], 'costateCols', [2 4], ...
                                     'quotientDir', [], 'freeTime', false));
ok = chk(ok, out.pass && out.nCrossings == 0 && all(isfinite(out.detScaled)) && out.tested, ...
         sprintf('underflowing raw det (2x2 block, entries 1e-160): no fictitious roots (n=%d, pass=%d)', ...
                 out.nCrossings, out.pass));

% (f) an EXACT-ZERO sample adjacent to a sign flip is ONE root, not two:
%     LQ T = 4 with the junction nearest pi forced to det == 0.
%     Fixture: three quarter-turn rotations (exact matrices). Phi_xl at the
%     samples: -1 (t = pi/2), 0 EXACTLY (t = pi, [0 -1;1 0]^2 = -I), +1.
%     Samples [-, 0, +]: one root at pi, not "one flip plus one zero".
tg = [0 pi/2 pi 3*pi/2];
info = struct('PHI', {{[0 -1; 1 0], [0 -1; 1 0], [0 -1; 1 0]}}, 'Y', zeros(2,3), 'tGrid', tg);
out = ms_conjugate_test(info, spec);
ok = chk(ok, out.nCrossings == 1, ...
         sprintf('[-, 0, +] samples (exact zero at pi, flip across it) count ONE root (n=%d)', out.nCrossings));

% (g) a root ONLY on the bracket ending at t_f is ENDPOINT, not a refutation:
%     LQ with T = pi + tiny: the crossing sits at pi, inside the last bracket.
K = 48;  T = pi + 1e-3;
info = lqInfo(T, K);
out = ms_conjugate_test(info, spec);
ok = chk(ok, out.atFinal && out.nInterior == 0 && strcmp(out.verdict, 'ENDPOINT') && ~out.pass, ...
         sprintf('root in the last bracket only -> ENDPOINT (atFinal=%d nInterior=%d verdict=%s)', ...
                 out.atFinal, ternary(isfield(out,'nInterior'), out.nInterior, -1), ...
                 ternary(isfield(out,'verdict'), out.verdict, '?')));
% and sampledThrough is reported:
ok = chk(ok, isfield(out, 'sampledThrough') && abs(out.sampledThrough - T) < 1e-12, ...
         'sampledThrough == t_f when Yend is not needed (fixed tf)');

if ok, fprintf('TEST_CONJ_FIXEDTF: ALL PASS\n');
else,  fprintf('TEST_CONJ_FIXEDTF: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function v = ternary(c, a, b)
% TERNARY  a if c else b.  INPUTS: c logical; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

% ------------------------------------------------------------------------
function info = lqInfo(T, K)
% LQINFO  ms_bvp-shaped info (PHI segment STMs, Y junctions, tGrid) for
% the LQ accessory system ydot = [ -y(2); y(1) ] (rotation flow).
% INPUTS: T horizon; K segments.  OUTPUTS: info struct.
tG = linspace(0, T, K+1);
h  = tG(2) - tG(1);
Phi = [cos(h), -sin(h); sin(h), cos(h)];       % exact segment STM
info = struct('PHI', {repmat({Phi}, 1, K)}, ...
              'Y', zeros(2, K), 'tGrid', tG);
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
