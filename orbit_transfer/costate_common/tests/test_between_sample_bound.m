function ok = test_between_sample_bound()
% TEST_BETWEEN_SAMPLE_BOUND  A lower bound on a function's minimum OVER AN
% INTERVAL from its samples and a bound on its slope, against a function
% whose true minimum is known and falls BETWEEN the samples:
%
%       v(t) = 1.2 + cos(3 t),  min = 0.2 at t = pi/3,  |v'| <= 3 .
%
% The sampled minimum overstates the true one; the bound must never do so,
% and must close on it as the sampling is refined. Then the two hypothesis
% gates: mintime_hypothesis_gates reports the bounds beside the sampled
% minima (H2: |lam_v|, H3: Q_mt) on the library's anchor.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % costate_common
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
v = @(t) 1.2 + cos(3*t);  trueMin = 0.2;

% ---- coarse sampling: the sampled minimum misses, the bound does not ------
t = linspace(0, 2, 12).';                       % pi/3 = 1.047 falls between samples
[lo, k] = between_sample_bound(t, v(t), 3*ones(size(t)));
ok = chk(ok, min(v(t)) > trueMin + 1e-3, sprintf('the sampled minimum %.4f overstates the true one %.4f', min(v(t)), trueMin));
ok = chk(ok, lo <= trueMin, sprintf('the bound %.4f does not', lo));
ok = chk(ok, t(k) <= pi/3 && pi/3 <= t(k+1), 'and it names the interval that holds the true minimum');
% ---- it is no looser than slope * step / 2 below the samples ---------------
ok = chk(ok, lo >= min(v(t)) - 3*max(diff(t))/2 - 1e-12, 'never looser than L dt / 2 below the sampled minimum');
% ---- refinement closes the gap ---------------------------------------------
t = linspace(0, 2, 4001).';
lo = between_sample_bound(t, v(t), 3*ones(size(t)));
ok = chk(ok, lo <= trueMin && trueMin - lo < 3*max(diff(t))/2 + 1e-9, sprintf('fine sampling: bound %.6f within L dt/2 of %.1f', lo, trueMin));
% ---- a slope bound that varies along the samples: the larger end is used ---
t = [0; 1];  [lo1, ~] = between_sample_bound(t, [1; 1], [0; 2]);
ok = chk(ok, abs(lo1 - (1 - 2*1/2)) < 1e-12, 'per interval the LARGER of its two slope bounds applies');
ok = chk(ok, throws(@() between_sample_bound([0; 1], [1; 1; 1], [1; 1])) && throws(@() between_sample_bound([0; 0], [1; 1], [1; 1])), ...
         'mismatched sizes and non-increasing times are refused');

% ---- the gates report the bounds (the library's anchor) ---------------------
lib = dro_tulip_library();  E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup', struct('physicsOnly', true));
rv0 = B.stateD(E.sD);
g = mintime_hypothesis_gates(E.z(:), rv0(1:6), B.Tnd, B.cnd, B.mu, struct());
ok = chk(ok, all(isfield(g, {'minLamVBound', 'minQmtBound', 'dtMax'})), 'mintime_hypothesis_gates reports .minLamVBound .minQmtBound .dtMax');
if all(isfield(g, {'minLamVBound', 'minQmtBound'}))
    ok = chk(ok, g.minLamVBound <= g.minLamV && g.minLamVBound > 0.95*g.minLamV, sprintf('H2: bound %.4f just under the sampled %.4f', g.minLamVBound, g.minLamV));
    ok = chk(ok, g.minQmtBound <= g.minQmt && g.minQmtBound > 0.95*g.minQmt, sprintf('H3: bound %.4f just under the sampled %.4f', g.minQmtBound, g.minQmt));
end
if ok, fprintf('test_between_sample_bound: ALL PASS\n'); else, fprintf('test_between_sample_bound: FAIL\n'); end
end

function tf = throws(f)
% THROWS  Does calling f throw?  INPUTS: f (handle).  OUTPUTS: tf.
tf = false;
try, f(); catch, tf = true; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
