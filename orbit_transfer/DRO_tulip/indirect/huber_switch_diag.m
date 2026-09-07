function D = huber_switch_diag(Y1, tf, Tmax, c, muStar, smooth)
% HUBER_SWITCH_DIAG  Switch-structure diagnostics of a converged (or stalled)
% smoothed-fuel extremal: how many times Q crosses 1, how transversal the
% least transversal crossing is, and whether a local extremum of Q hovers
% just short of 1 -- an IMMINENT GRAZING BIFURCATION (a burn arc about to be
% born or annihilated as the smoothing parameter moves).
%
% Built for the Huber-wall diagnosis of 2026-09-06 (FINDINGS 24): both
% Huber walls showed exactly these signatures (a near-tangent crossing with
% |dQ/dt| = 0.045, and a Q maximum of 0.9907 without crossing), the clean
% cells none. Family-agnostic: works for 'eps' too.
%
% INPUTS:
%   Y1     - full PMP state at t = 0 [14x1]: [r; v; m; lam_r; lam_v; lam_m]
%   tf     - time of flight, ND [scalar]
%   Tmax, c, muStar - as cr3bp_minfuel_pmp
%   smooth - struct('family', 'eps'|'huber', 'p', value)
%
% OUTPUTS:
%   D - struct: .T [nT x 1] times, .Q [nT x 1] switch quantity along the arc,
%       .Qdot [nT x 1] EXACT dQ/dt (cr3bp_minfuel_qdot), .tCross [1 x nC]
%       interpolated crossing times of Q = 1, .dQdtCross [1 x nC] exact
%       Qdot there (transversality; small = near-grazing), .minAbsDQdt,
%       .tTouch tangential touches (Q == 1 without a sign change),
%       .graze struct array of grazing-RISK extrema (.t, .Q, .kind
%       'max-below-1' | 'min-above-1', .Qdot) -- a screen, not a verdict,
%       .nCross, .nTouch, .note
%
% REFERENCES:
%   [1] DRO_tulip/FINDINGS.md section 24 (the measurement this packages).
%   [2] costate_common/cr3bp_minfuel_prop.m (saltation denominator n'F- is
%       dQ/dt at the crossing; grazing = that denominator -> 0).

[~, ~, T, Yt] = cr3bp_minfuel_prop(tf, Y1(:), false, Tmax, c, muStar, smooth);
% drop (near-)duplicate event samples: tolerance-aware, keep the first
keep = [true; diff(T) > 1e-12 * max(tf, 1)];
T = T(keep);  Yt = Yt(keep, :);
Q  = Tmax * (sqrt(sum(Yt(:, 11:13).^2, 2)) ./ Yt(:, 7) + Yt(:, 14) / c);
Qd = cr3bp_minfuel_qdot(Yt', Tmax, c)';           % EXACT dQ/dt at every sample

% Crossings vs touches (Astra review #2): a sample with Q == 1 exactly is a
% crossing only if the signs on its two sides differ; same sign = a TOUCH.
g = Q - 1;  sg = sign(g);
kx = [];  touch = [];
n = numel(sg);  k = 1;
while k < n
    if sg(k) ~= 0 && sg(k+1) ~= 0
        if sg(k) ~= sg(k+1), kx(end+1) = k; end                 %#ok<AGROW> crossing in (k, k+1)
        k = k + 1;
    elseif sg(k+1) == 0
        k2 = k + 1;  while k2 < n && sg(k2+1) == 0, k2 = k2 + 1; end   % zero run k+1..k2
        prevS = sg(k);  nextS = 0;  if k2 < n, nextS = sg(k2+1); end
        if prevS ~= 0 && nextS ~= 0 && prevS ~= nextS, kx(end+1) = k+1;   %#ok<AGROW> crossing AT the sample
        elseif prevS ~= 0 && nextS ~= 0,              touch(end+1) = k+1; %#ok<AGROW> tangency
        end
        k = k2 + 1;
    else
        k = k + 1;
    end
end
% crossing time by linear interpolation of g, transversality = exact Qdot there
tCross = zeros(1, numel(kx));  dQdtCross = zeros(1, numel(kx));
for j = 1:numel(kx)
    a = kx(j);
    if g(a) == 0 || a == n
        tCross(j) = T(a);  dQdtCross(j) = Qd(a);
    else
        w = g(a) / (g(a) - g(a+1));
        tCross(j) = T(a) + w*(T(a+1) - T(a));
        dQdtCross(j) = (1-w)*Qd(a) + w*Qd(a+1);
    end
end
minAbsDQdt = NaN;  if ~isempty(kx), minAbsDQdt = min(abs(dQdtCross)); end

% Grazing-RISK screen (not a bifurcation verdict): local extrema of Q
% located as sign changes of the EXACT Qdot, not at a crossing, with Q
% within a configurable window of 1.
win = 0.1;
ext = find(diff(sign(Qd)) ~= 0) + 1;
ext = ext(arrayfun(@(e) all(abs(T(e) - tCross) > 1e-3 * tf), ext));
isMax = Qd(max(ext-1, 1)) > 0;
gz = ext((isMax & Q(ext) < 1 & Q(ext) > 1 - win) | (~isMax & Q(ext) > 1 & Q(ext) < 1 + win));
graze = struct('t', {}, 'Q', {}, 'kind', {}, 'Qdot', {});
for e = gz(:)'
    k = 'min-above-1';  if Q(e) < 1, k = 'max-below-1'; end
    graze(end+1) = struct('t', T(e), 'Q', Q(e), 'kind', k, 'Qdot', Qd(e)); %#ok<AGROW>
end

D = struct('T', T, 'Q', Q, 'Qdot', Qd, 'tCross', tCross, 'dQdtCross', dQdtCross, ...
           'minAbsDQdt', minAbsDQdt, 'tTouch', T(touch)', 'graze', graze, ...
           'nCross', numel(kx), 'nTouch', numel(touch), ...
           'note', 'graze = risk screen (extremum of Q within a window of 1), not a proven bifurcation');
end
