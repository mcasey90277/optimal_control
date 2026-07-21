% P0A_GRAZE_MARGIN  Preflight: switching-slope margin at every PSR bang-bang
% solution on disk (PLAN_PRONG_Z.md P0a).
%
% For each PSR_data/psr_data_tf*_minEps0.mat, take the dual-mapped switching
% function S(sigma) (mesh accuracy ~1%), locate its sign crossings, and
% estimate the local slope |dS/dtau| (Sundman time) and |dS/dt| (physical
% time) at each crossing. The MINIMUM over crossings is the graze margin that
% the ZTL saltation matrix Psi = I + (f+ - f-)(dS/dy)/Sdot divides by --
% Zhang's known weakness. This pre-registers the Z4/Z5 risk with a number and
% sets the graze-guard floor.
%
% Slope estimator: least-squares line through the +/-W nodes bracketing each
% crossing (W=3), robust to single-node dual noise.
%
% Output: table per file + summary; saved to results/p0a_graze_margin.mat.

here = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', 'PSR_data');
resDir  = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

d = dir(fullfile(dataDir, 'psr_data_tf*_minEps0.mat'));
assert(~isempty(d), 'no bang-bang PSR data products found in %s', dataDir);

W = 3;                                   % half-width of the slope fit stencil
rows = struct('file',{},'factor',{},'nCross',{}, ...
              'minAbsdSdtau',{},'medAbsdSdtau',{},'tauAtMin',{}, ...
              'minAbsdSdt',{},'medAbsdSdt',{},'ratioMinMed_tau',{});

fprintf('=== P0a: graze margin from PSR bang-bang products ===\n');
for kf = 1:numel(d)
    D   = load(fullfile(dataDir, d(kf).name));
    S   = D.costate.S(:);
    tau = D.mesh.tau(:);
    t   = D.mesh.t(:);
    n   = numel(S);

    cr = find(S(1:end-1).*S(2:end) < 0);          % bracketing indices
    dSdtau = zeros(numel(cr),1);  dSdt = zeros(numel(cr),1);
    tauC   = zeros(numel(cr),1);
    for kc = 1:numel(cr)
        idx = max(1, cr(kc)-W+1) : min(n, cr(kc)+W);
        ptau = polyfit(tau(idx) - tau(cr(kc)), S(idx), 1);
        pt   = polyfit(t(idx)   - t(cr(kc)),   S(idx), 1);
        dSdtau(kc) = ptau(1);
        dSdt(kc)   = pt(1);
        % crossing location from the local line
        tauC(kc) = tau(cr(kc)) - ptau(2)/ptau(1);
    end

    [mn, im] = min(abs(dSdtau));
    rows(end+1) = struct('file', d(kf).name, 'factor', D.factor, ...
        'nCross', numel(cr), ...
        'minAbsdSdtau', mn, 'medAbsdSdtau', median(abs(dSdtau)), ...
        'tauAtMin', tauC(im), ...
        'minAbsdSdt', min(abs(dSdt)), 'medAbsdSdt', median(abs(dSdt)), ...
        'ratioMinMed_tau', mn/median(abs(dSdtau))); %#ok<SAGROW>

    fprintf(['%-38s f=%.3f  cross=%2d  |dS/dtau|: min=%.3e med=%.3e ' ...
             '(min/med=%.2f, at tau=%.2f/%.2f)  |dS/dt|: min=%.3e med=%.3e\n'], ...
        d(kf).name, D.factor, numel(cr), mn, median(abs(dSdtau)), ...
        mn/median(abs(dSdtau)), tauC(im), tau(end), ...
        min(abs(dSdt)), median(abs(dSdt)));
end

save(fullfile(resDir, 'p0a_graze_margin.mat'), 'rows');
fprintf('saved %s\n', fullfile(resDir, 'p0a_graze_margin.mat'));

% Verdict guidance: saltation is safe when min|Sdot| is not << the typical
% scale. Flag any file whose min/med ratio < 0.05 as graze-risky.
risky = rows([rows.ratioMinMed_tau] < 0.05);
if isempty(risky)
    fprintf('VERDICT: no graze-risky crossings (all min/med >= 0.05).\n');
else
    fprintf('VERDICT: %d file(s) with graze-risky crossings (min/med < 0.05):\n', numel(risky));
    fprintf('  %s\n', risky.file);
end
