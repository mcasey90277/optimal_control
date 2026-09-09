function ok = test_conjugate_pole_predict()
%% Purpose:
%
%   Tests conjugate_pole_predict -- the EARLY WARNING that a continuation is
%   walking into a conjugate point, fitted from the conditioning it already
%   computes at every step.
%
%   The physics it rests on (FINDINGS 39): a conjugate point is a
%   non-trivial solution of the linearised BVP, so the multiple-shooting
%   Jacobian goes singular there and cond(J) has a POLE. Log-linear
%   extrapolation therefore under-shoots badly; a pole fit does not.
%
%   The fixture is the measured series from the DRO -> tulip departure walk
%   at 70 mN, arrival phase 0.0754, whose conjugate point was located at
%   sD = 0.04665 +- 0.00005 by direct bisection of the conjugate verdict.
%   The test hands the predictor only the points up to sD = 0.0450 -- that
%   is 0.0017 of a period SHORT of the answer, and two decades of cond(J)
%   before the verdict flips -- and requires the prediction to land.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

% measured 2026-09-09 (asym3 probe): departure phase, cond(J), conj verdict
s    = [0.0000 0.0100 0.0200 0.0300 0.0380 0.0420 0.0440 0.0450 0.0460 0.0465 0.0466 0.0467];
cj   = [4.401e9 4.668e9 5.393e9 9.197e9 2.785e10 8.601e10 2.350e11 5.272e11 ...
        2.143e12 9.050e12 1.458e13 2.734e13];
sTrue = 0.04665;                       % located by bisecting the verdict

% (1) predict from data that stops well short of the answer
m = s <= 0.0450;
W = conjugate_pole_predict(s(m), cj(m));
ok = chk(ok, W.ok, sprintf('a pole is detected from data ending at 0.0450 (%s)', W.reason));
ok = chk(ok, abs(W.sPole - sTrue) < 1e-3, ...
         sprintf('predicted %.5f vs measured %.5f (error %.1e of a period)', ...
                 W.sPole, sTrue, abs(W.sPole - sTrue)));
ok = chk(ok, W.sPole > s(find(m, 1, 'last')), 'the pole is predicted AHEAD of the data');

% (2) a log-linear extrapolation is much worse -- this is why the fit is a pole
yy = log10(cj(m));  xx = s(m);
p = polyfit(xx(end-2:end), yy(end-2:end), 1);
sLin = (log10(2.734e13) - p(2))/p(1);
ok = chk(ok, abs(sLin - sTrue) > 3*abs(W.sPole - sTrue), ...
         sprintf('pole fit beats log-linear: %.5f vs %.5f (truth %.5f)', W.sPole, sLin, sTrue));

% (3) a benign series must NOT raise a warning
sB = linspace(0, 0.05, 10);
W2 = conjugate_pole_predict(sB, 4e9*(1 + 0.3*sB/0.05));
ok = chk(ok, ~W2.ok, sprintf('a flat series raises nothing (%s)', W2.reason));

% (4) too few points is refused, not guessed
W3 = conjugate_pole_predict(s(1:2), cj(1:2));
ok = chk(ok, ~W3.ok && contains(lower(W3.reason), 'point'), ...
         sprintf('too few points refused: %s', W3.reason));

if ok, fprintf('TEST_CONJUGATE_POLE_PREDICT: ALL PASS\n'); else, fprintf('TEST_CONJUGATE_POLE_PREDICT: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
