function ok = test_lift_margin()
%% Purpose:
%
%   Tests lift_margin -- the replacement for the heuristic rank threshold
%   behind dim S.
%
%   The current rule counts singular values below a tolerance built from
%   three ad hoc numbers, and its header claims the cap prevents a noisy
%   lift from reporting a spurious nullity. It prevents no such thing
%   (Astra, 2026-09-10). What CAN be claimed is a theorem plus a measurement:
%
%     dim S >= 1   is CONSTRUCTIVE -- we exhibit the lift.
%     dim S <= 1   follows from Eckart-Young: the distance from C to the
%                  nearest rank-deficient matrix is sigma_6, so if
%                  sigma_6 > ||dC|| then rank(C) = 6 EXACTLY.
%
%   So the honest output is a MARGIN, sigma_6 / ||dC||, with ||dC|| measured
%   by rebuilding C at a second numerical setting -- not a pass against a
%   number nobody can defend.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(7);

% a matrix of known rank 6 in R^7, with a known null direction
[U, ~] = qr(randn(20, 7), 0);
sv = [3 1 0.4 0.2 0.05 0.01 0].';
[V, ~] = qr(randn(7));
C = U*diag(sv)*V.';
lam = V(:, 7);                                  % the exact null vector

% (1) a clean pair: the two builds differ by almost nothing -> huge margin
M = lift_margin(C, C + 1e-14*randn(20,7), lam, struct());
ok = chk(ok, M.dimS == 1, sprintf('rank 6 recognised: dim S = %d', M.dimS));
ok = chk(ok, M.margin > 1e3, sprintf('clean pair gives a large margin: %.2e', M.margin));
ok = chk(ok, M.certified, sprintf('and is certified (%s)', M.reason));
ok = chk(ok, abs(M.sigma6 - 0.01) < 1e-9, sprintf('sigma_6 read correctly: %.4g', M.sigma6));

% (2) noise COMPARABLE to sigma_6 must destroy the certificate, not the answer
M2 = lift_margin(C, C + 0.02*randn(20,7)/sqrt(20), lam, struct());
ok = chk(ok, ~M2.certified && M2.margin < 10, ...
         sprintf('noise near sigma_6 removes the certificate: margin %.2f (%s)', M2.margin, M2.reason));

% (3) the constructive half must be reported separately and honestly
ok = chk(ok, M.nullResid < 1e-12, sprintf('exhibited lift residual %.1e', M.nullResid));
Mz = lift_margin(C, C + 1e-14*randn(20,7), zeros(7,1), struct());
ok = chk(ok, ~Mz.certified && contains(lower(Mz.reason), 'zero'), ...
         sprintf('a zero vector is not a lift: %s', Mz.reason));
Mn = lift_margin(C, C + 1e-14*randn(20,7), [1; NaN; zeros(5,1)], struct());
ok = chk(ok, ~Mn.certified && contains(lower(Mn.reason), 'finite'), 'a non-finite lift is refused');
Mb = lift_margin(C, C, randn(7,1), struct());
ok = chk(ok, ~Mb.certified && contains(lower(Mb.reason), 'lift'), ...
         sprintf('a vector that is NOT a lift is refused: %s', Mb.reason));

% (4) genuine rank 5 (a two-dimensional null space) must not read as dim S = 1
sv2 = [3 1 0.4 0.2 0.05 0 0].';
C2 = U*diag(sv2)*V.';
M3 = lift_margin(C2, C2 + 1e-14*randn(20,7), V(:,6), struct());
ok = chk(ok, M3.dimS == 2, sprintf('a 2-D null space reads dim S = %d', M3.dimS));
ok = chk(ok, ~M3.certified, 'and is not certified as dim S = 1');

if ok, fprintf('TEST_LIFT_MARGIN: ALL PASS\n'); else, fprintf('TEST_LIFT_MARGIN: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
