function ok = test_gates_h6_wiring()
%% Purpose:
%
%   Tests that H6 is actually WIRED INTO mintime_hypothesis_gates, not merely
%   available beside it.
%
%   h6_margin exists and is tested on its own; this checks the thing that
%   matters operationally -- that every entry certified from today on CARRIES
%   the H6 numbers, so the reduced problem's spurious-zero mechanism is
%   excluded per entry rather than argued in a document. An instrument nobody
%   calls is the same failure mode as a gate computed and not enforced
%   (FINDINGS 36), one step earlier.
%
%   It also checks the SECOND half of the wiring: that the constraint matrix
%   is exposed so lift_margin can be built at two settings, which is what
%   turns the rank threshold into a measured margin.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

[B, anc] = arclength_arrival('setup');
C = certify_crossing(anc.p, anc.sA, B, anc);
assert(C.ok, 'fixture: the anchor must certify (%s)', C.reason);
z8 = C.z;  rv0 = B.rv0(1:6);

g = mintime_hypothesis_gates(z8, rv0, B.Tnd, B.cnd, B.mu, struct());

% (1) H6 must be present and must agree with the standalone unit
ok = chk(ok, isfield(g, 'h6Margin') && isfield(g, 'h6Ok'), 'the gates carry H6');
H = h6_margin(z8, B.Tnd, B.cnd);
if isfield(g, 'h6Margin')
    ok = chk(ok, abs(g.h6Margin - H.margin) < 1e-12 && g.h6Ok == H.ok, ...
             sprintf('and it agrees with h6_margin: %.4fx', g.h6Margin));
    ok = chk(ok, g.h6Ok && g.h6Margin > 1, sprintf('the anchor passes H6 (%.1fx)', g.h6Margin));
end

% (2) the constraint matrix must be obtainable, at a requested tolerance,
%     so the rank statement can become a MEASURED margin
g2 = mintime_hypothesis_gates(z8, rv0, B.Tnd, B.cnd, B.mu, struct('keepC', true));
ok = chk(ok, isfield(g2, 'C') && size(g2.C, 2) == 7, ...
         'the constraint matrix C is exposed on request');
g3 = mintime_hypothesis_gates(z8, rv0, B.Tnd, B.cnd, B.mu, ...
                              struct('keepC', true, 'relTol', 1e-7));
if isfield(g2, 'C') && isfield(g3, 'C') && isequal(size(g2.C), size(g3.C))
    M = lift_margin(g2.C, g3.C, z8(1:7), struct());
    ok = chk(ok, isfinite(M.margin), ...
             sprintf('lift_margin runs on the two builds: sigma_6 %.2e, error %.2e, margin %.1fx', ...
                     M.sigma6, M.errEst, M.margin));
    ok = chk(ok, M.dimS >= 1, sprintf('and reports dim S = %d (%s)', M.dimS, M.reason));
end

if ok, fprintf('TEST_GATES_H6_WIRING: ALL PASS\n'); else, fprintf('TEST_GATES_H6_WIRING: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
