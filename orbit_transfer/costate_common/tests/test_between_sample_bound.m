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
ok = chk(ok, all(isfield(g, {'minLamVEstimate', 'minQmtEstimate', 'dtMax'})) && ~any(isfield(g, {'minLamVBound', 'minQmtBound'})), ...
         'the gates report .minLamVEstimate .minQmtEstimate .dtMax (an ESTIMATE is not called a bound)');
if all(isfield(g, {'minLamVEstimate', 'minQmtEstimate'}))
    ok = chk(ok, g.minLamVEstimate <= g.minLamV && g.minLamVEstimate > 0.95*g.minLamV, sprintf('H2: estimate %.4f just under the sampled %.4f', g.minLamVEstimate, g.minLamV));
    ok = chk(ok, g.minQmtEstimate <= g.minQmt && g.minQmtEstimate > 0.95*g.minQmt, sprintf('H3: estimate %.4f just under the sampled %.4f', g.minQmtEstimate, g.minQmt));
end
% ---- the slope rule, against the flight itself (an oracle, not the gate's own algebra) ----
% The Coriolis term of lam_v' is skew, so it cannot change |lam_v|: |d|lam_v|/dt| <= |lam_r|.
% Checked by differencing |lam_v| along an independently flown arc.
[tt, Y] = pumpkyn.cr3bp.tfMinProp(E.z(8), [rv0(1:6); 1; E.z(1:7)], B.Tnd, B.cnd, B.mu);
[tt, iu] = unique(tt);  Y = Y(iu, :);
rho = sqrt(sum(Y(:, 11:13).^2, 2));  lamR = sqrt(sum(Y(:, 8:10).^2, 2));
kk = (2:numel(tt)-1).';
dRho = (rho(kk+1) - rho(kk-1)) ./ (tt(kk+1) - tt(kk-1));
worst = max(abs(dRho) ./ max(lamR(kk-1), max(lamR(kk), lamR(kk+1))));
ok = chk(ok, worst <= 1 + 1e-3, sprintf('on the flown anchor |d|lam_v|/dt| <= |lam_r| (largest ratio %.4f)', worst));
if isfield(g, 'minLamVEstimate')
    gap = g.minLamV - g.minLamVEstimate;
    ok = chk(ok, gap <= 0.5*max(lamR)*g.dtMax + 1e-12, sprintf('and the gate uses that rule: gap %.2e <= max|lam_r| dtMax/2 = %.2e', gap, 0.5*max(lamR)*g.dtMax));
end
% ---- inputs a bound cannot be made from --------------------------------------
t2 = [0; 1; 2];
ok = chk(ok, throws(@() between_sample_bound(t2, [1; NaN; 1], [1; 1; 1])) && throws(@() between_sample_bound(t2, [1; 1; 1], [1; Inf; 1])) ...
          && throws(@() between_sample_bound([0; NaN; 2], [1; 1; 1], [1; 1; 1])), 'non-finite samples, slopes or times are refused');
ok = chk(ok, throws(@() between_sample_bound(t2, [1; 1; 1], [1; -1; 1])), 'a negative slope bound is refused');
ok = chk(ok, throws(@() between_sample_bound(t2, [1; 1i; 1], [1; 1; 1])), 'complex samples are refused');
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
