function ok = test_mintime_gates()
%% Purpose:
%
%   Tests mintime_hypothesis_gates -- the per-entry checks of the
%   Bonnard-Caillau-Trelat sufficiency hypotheses the min-time catalog
%   test had been assuming (doc/mintime_second_order_audit.tex, section 4):
%
%     H2  strong Legendre: min_t |lam_v(t)| > 0 along the flown arc
%     H3  all-burn is the PMP control: min_t Q_mt(t) > 0,
%         Q_mt = |lam_v|/m + lam_m/c
%     H1' no abnormal lift: the linear space S of costate lifts of the
%         flown trajectory (fixed-control adjoint, lam_v parallel to the
%         flown thrust direction at every sample, lam_m(tf) = 0) is
%         ONE-dimensional. Self-consistency: the accepted lam0 must lie
%         in S (constraint residual ~ 0), and lam.f == -1 along the arc
%         (H == 0 with lambda_0 = 1).
%
%   Fixtures: the three engine golden cells (dro/halo/dpo, certified 1 N
%   entries) -- all three must pass all gates.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
G = load(fullfile(here, 'golden_cells_data.mat'));
for kc = 1:numel(G.cells)
    c = G.cells(kc);
    g = mintime_hypothesis_gates(c.z8, c.rv0, c.Tnd, c.cnd, c.muStar);
    ok = chk(ok, g.minLamV > 1e-3, sprintf('%s H2: min|lam_v| = %.3e > 0 (at t/tf = %.2f)', c.name, g.minLamV, g.tMinLamV/c.z8(8)));
    ok = chk(ok, g.minQmt > 0,     sprintf('%s H3: min Q_mt = %.3e > 0 (at t/tf = %.2f)', c.name, g.minQmt, g.tMinQmt/c.z8(8)));
    ok = chk(ok, g.dimS == 1,      sprintf('%s H1'': dim S = %d (sv ratio s7/s6 = %.1e)', c.name, g.dimS, g.svRatio));
    ok = chk(ok, g.nullResid < 1e-6, sprintf('%s probe self-consistency: |C lam0|/|lam0| = %.1e', c.name, g.nullResid));
    ok = chk(ok, g.Hresid < 1e-6,  sprintf('%s normal lift: max|lam.f + 1| = %.1e', c.name, g.Hresid));
end
if ok, fprintf('TEST_MINTIME_GATES: ALL PASS\n');
else,  fprintf('TEST_MINTIME_GATES: FAILURE (see lines above)\n');
end
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
