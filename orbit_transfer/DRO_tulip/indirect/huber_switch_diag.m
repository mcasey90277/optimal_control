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
%       .tCross [1 x nC] crossing times of Q = 1, .dQdtCross [1 x nC] dQ/dt
%       there (transversality; small = near-grazing), .minAbsDQdt,
%       .graze struct array of imminent grazes (.t, .Q, .kind
%       'max-below-1' | 'min-above-1'), .nCross
%
% REFERENCES:
%   [1] DRO_tulip/FINDINGS.md section 24 (the measurement this packages).
%   [2] costate_common/cr3bp_minfuel_prop.m (saltation denominator n'F- is
%       dQ/dt at the crossing; grazing = that denominator -> 0).

[~, ~, T, Yt] = cr3bp_minfuel_prop(tf, Y1(:), false, Tmax, c, muStar, smooth);
[T, iu] = unique(T, 'stable');  Yt = Yt(iu, :);     % drop duplicated event samples
Q = Tmax * (sqrt(sum(Yt(:, 11:13).^2, 2)) ./ Yt(:, 7) + Yt(:, 14) / c);

s = sign(Q - 1);  s(s == 0) = 1;
kx = find(diff(s) ~= 0);
dQ = gradient(Q, T);
tCross = T(kx)';  dQdtCross = dQ(kx)';
minAbsDQdt = NaN;  if ~isempty(kx), minAbsDQdt = min(abs(dQdtCross)); end

% local extrema not coincident with a crossing, within 10% of the jump
d2  = diff(sign(diff(Q)));  ext = find(d2 ~= 0) + 1;
ext = ext(arrayfun(@(e) all(abs(T(e) - T(kx)) > 1e-3 * tf), ext));
isMax = d2(ext - 1) < 0;
gz = ext((isMax & Q(ext) < 1 & Q(ext) > 0.9) | (~isMax & Q(ext) > 1 & Q(ext) < 1.1));
graze = struct('t', {}, 'Q', {}, 'kind', {});
for e = gz(:)'
    k = 'min-above-1';  if Q(e) < 1, k = 'max-below-1'; end
    graze(end+1) = struct('t', T(e), 'Q', Q(e), 'kind', k); %#ok<AGROW>
end

D = struct('T', T, 'Q', Q, 'tCross', tCross, 'dQdtCross', dQdtCross, ...
           'minAbsDQdt', minAbsDQdt, 'graze', graze, 'nCross', numel(kx));
end
