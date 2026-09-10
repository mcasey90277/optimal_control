function ok = test_h6_margin()
%% Purpose:
%
%   Tests h6_margin -- the gate that excludes the REDUCED problem's
%   spurious-zero mechanism.
%
%   In the 6-state reduction the reduced Hamiltonian is not conserved and
%   p'J = 0, so det[JP, F] = 0 can mean h(t) = 0 rather than a rank drop: a
%   zero that is NOT a loss of optimality, and one that none of H1-H4
%   excluded (Astra proof review, 2026-09-10). In closed form
%
%       lambda_rv . f_rv = -1 + (T/c) lambda_m
%
%   which vanishes exactly at lambda_m = c/T, and lambda_m decreases
%   monotonically to lambda_m(t_f) = 0. So the whole objection reduces to
%
%       H6:  lambda_m(0) < c/T.
%
%   By the identity proved the same day, integral(T Q_mt dt) = lambda_m(0),
%   this also says the total strict-bang margin must stay below c/T.
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

Tnd = 0.1756418;  cnd = 8.6737;          % the 70 mN operating point
thresh = cnd/Tnd;

% (1) a real catalog value: comfortable
H = h6_margin([zeros(6,1); 5.500404; 4.0151], Tnd, cnd);
ok = chk(ok, abs(H.threshold - thresh) < 1e-9, sprintf('threshold c/T = %.4f', H.threshold));
ok = chk(ok, H.ok, sprintf('the anchor passes (%s)', H.reason));
ok = chk(ok, abs(H.margin - thresh/5.500404) < 1e-9, sprintf('margin %.2fx', H.margin));
ok = chk(ok, abs(H.hMin - (-1 + (Tnd/cnd)*5.500404)) < 1e-12, ...
         sprintf('the worst |lambda_rv.f_rv| is reported: %.4f', H.hMin));

% (2) exactly at the threshold: must FAIL, not round in our favour
H2 = h6_margin([zeros(6,1); thresh; 4.0], Tnd, cnd);
ok = chk(ok, ~H2.ok && abs(H2.margin - 1) < 1e-9, ...
         sprintf('lambda_m(0) = c/T exactly is refused (margin %.4f)', H2.margin));

% (3) beyond it: the spurious mechanism CAN fire
H3 = h6_margin([zeros(6,1); 2*thresh; 4.0], Tnd, cnd);
ok = chk(ok, ~H3.ok && H3.margin < 1, sprintf('beyond the threshold is refused: %s', H3.reason));

% (4) a negative mass costate is not admissible here -- lambda_m decreases to
%     zero, so lambda_m(0) < 0 means the transversality condition is violated
H4 = h6_margin([zeros(6,1); -1; 4.0], Tnd, cnd);
ok = chk(ok, ~H4.ok && contains(lower(H4.reason), 'negative'), ...
         sprintf('a negative lambda_m(0) is refused: %s', H4.reason));

if ok, fprintf('TEST_H6_MARGIN: ALL PASS\n'); else, fprintf('TEST_H6_MARGIN: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
