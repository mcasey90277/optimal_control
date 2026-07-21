function dataFile = psr_export_data(solFile, dataDir, opts)
% PSR_EXPORT_DATA  Generate costates + save the PSR data products to PSR_data.
%
% Takes a converged direct/PSR solution (seed layout .mat) and writes ONE
% self-contained data-product file for independent analysis and for seeding
% the IFS method. Two layers in the same file:
%
%  (a) SEED-COMPATIBLE layer (top-level: out, sigma, tauf0, rv0, rvf, factor)
%      -- byte-for-byte the layout ifs_seed / sms_seed_duals /
%      verify_direct_pmp / prep_refine_seed already consume, so this file IS
%      a valid IFS seed with zero conversion.
%
%  (b) UNPACKED data products (structs, all documented below):
%      mesh    - sigma [1xnN] normalized mesh; tau [1xnN] Sundman variable
%                (= sigma*tauf0); tauf0; t [1xnN] physical time at nodes;
%                tDays; pSund; nN
%      traj    - r [3xnN], v [3xnN], m [1xnN] (mass fraction), X [8xnN] the
%                full state [r;v;m;t]
%      ctrl    - alpha [3xnN] thrust direction, s [1xnN] throttle,
%                nSwitchS (dual-S crossing count -- the CERTIFIED counting),
%                tauSwitchS [1xk] sub-node S=0 crossing times,
%                nSwitchThrottle / tauSwitchThrottle (raw s=0.5 crossings --
%                can overcount on shallow dips, kept for cross-check)
%      costate - lam [8xnN] = [lamR;lamV;lamM;lamT] node costates from the
%                NLP KKT defect duals via the ADJUDICATED mode-'d' midpoint
%                map, scaled by the fitted beta (PSR does not natively
%                produce costates; this is the recovery machinery). S [1xnN]
%                switching function 1-||lamV||c/m-lamM; beta; spreadPct
%                (beta-fit consistency, %); accuracy note: these are O(h)
%                mesh-accuracy costates (~1%), NOT integration-accuracy --
%                see ifs/RESULTS_RUNG01_RUNG2.md Rung A for the proof that
%                no better single-trajectory recovery exists at 40 revs.
%      pmp     - first-order condition errors computable without a full
%                verify run: lamM_end (free-final-mass transversality,
%                want 0), termPosErr/termVelErr (rendezvous residual, ND),
%                termTimeErr (t(end)-tf, ND), SsignAgree (% of nodes where
%                sign(S) matches the throttle), primerAlignDeg (NLP-reported
%                mean primer alignment, deg, NaN if absent)
%      scal    - factor, tf, tf_days, dV (km/s), prop_kg, mf, maxDefect,
%                switches (NLP edge-count), edge
%      const   - full cr3bp_lt_params struct (muStar, lStar, tStar, Tmax, c,
%                m0kg, ...) + tfMin, so the file needs nothing else
%      provenance - date, source file, git hash, dual-map settings, pipeline
%
% The verify stage of run_psr APPENDS a `verify` struct (certificate summary)
% to the same file when it runs.
%
% INPUTS:
%   solFile - solution .mat in seed layout: out (X [8xnN], U [4xnN],
%             lamDef [8xN]), sigma, tauf0, rv0, rvf, factor
%   dataDir - destination directory (created if missing) [char]
%   opts    - (optional) struct: M dual-map arcs [default 40], quiet
%             suppress prints [default false], epsMin the homotopy endpoint
%             epsilon of the run [default 0] -- recorded in the filename
%
% OUTPUTS:
%   dataFile - written path, named
%              psr_data_tf<factor>_sw<k>_<insertionLabel>_minEps<e>.mat, e.g.
%              psr_data_tf1p150_sw25_tulipCampaign_minEps0.mat (bang-bang) or
%              psr_data_tf1p200_sw0_tulipCampaign_minEps0p01.mat (smooth
%              eps=0.01); factor and epsMin encode '.'->'p', k = certified
%              dual-S switch count, insertionLabel = insMeta.label (the
%              declared endpoint criterion, reconstructed via insertion_states
%              since this exporter is tulip-only) [char]. The tag is inserted
%              BEFORE _minEps so ztl/p0a_graze_margin.m's glob
%              ('psr_data_tf*_minEps0.mat') keeps matching.
%
% REFERENCES:
%   [1] ms_band/sms_seed_duals.m (mode-'d' dual->costate map, adjudicated
%       2026-07-10) + ms_band/beta_from_duals.m (scale fit).
%   [2] ifs/ifs_seed.m (the IFS consumer of this layout).
%   [3] PSR/run_psr.m section 4 (pipeline caller).

if nargin < 3, opts = struct(); end
if ~isfield(opts,'M'),      opts.M = 40;       end
if ~isfield(opts,'quiet'),  opts.quiet = false; end
if ~isfield(opts,'epsMin'), opts.epsMin = 0;   end
if ~exist(dataDir, 'dir'), mkdir(dataDir); end

% ---- declared insertion point (provenance + filename tag) -------------------
% This exporter is tulip-only (no target parameter), so the criterion is the
% pipeline's declared default ('campaign' -- see PSR/run_psr.m, psr_run_one.m).
% Reconstructed here (rather than threaded through opts) so this file's only
% callers (run_psr.m, psr_run_one.m) need no changes; the drift-guard assert
% below still catches a criterion mismatch loudly instead of silently
% mislabeling the export.
here = fileparts(mfilename('fullpath'));
if isempty(which('insertion_states'))
    addpath(fullfile(here, '..', 'sundman_minfuel'));
end
[rv0Decl, rvfDecl, insMeta] = insertion_states('tulip', 'campaign');

% ---- load the solution (seed layout) ----------------------------------------
Ssol   = load(solFile);
out    = Ssol.out;    sigma = Ssol.sigma;  tauf0 = Ssol.tauf0;
rv0    = Ssol.rv0;    rvf   = Ssol.rvf;    factor = Ssol.factor;
rvfC   = rvf(:);      % column view for residuals (legacy files store 1x6 rows;
                      % a row would broadcast X(1:3,end)-rvf(1:3) into a 3x3)

% drift guard: the solution being exported must match the reconstructed
% insertion point (catches a future criterion change loudly instead of
% silently mislabeling the export with the wrong 'insertion' tag).
assert(norm(rvf(:).' - rvfDecl) < 1e-10 && norm(rv0(:).' - rv0Decl) < 1e-10, ...
    'insertion:drift', ['solution endpoints differ from the reconstructed %s ' ...
    'insertion (rvf %.2e, rv0 %.2e) -- psr_export_data assumes the tulip ' ...
    'campaign default; update it if the pipeline default has changed'], ...
    insMeta.label, norm(rvf(:).'-rvfDecl), norm(rv0(:).'-rv0Decl));
X = out.X;  U = out.U;  nN = size(X, 2);
p   = cr3bp_lt_params(0.025, 15, 2100);
cfg = minfuel_config();
tf  = X(8, end);

% ---- costate generation (the machinery PSR itself does not run) -------------
% sms_seed_duals mode 'd': KKT defect duals -> node costates via the
% midpoint-principled map, with beta (positive scale + sign) fitted by the
% switching-law consistency fit (beta_from_duals). info.Y16 rows 9:16 are the
% scaled node costates on the solution's own node grid.
[~, ~, info] = sms_seed_duals(solFile, opts.M, 1e-4, 'd');
lam  = info.Y16(9:16, :);
tauN = info.tauN;                              % Sundman variable at nodes
if info.spreadPct > 5
    warning('psr_export_data:betaSpread', ...
        ['beta-fit spread %.1f%% > 5%% -- the dual->costate scale is poorly ' ...
         'determined for this solution (cf. the 1.85x failure); treat the ' ...
         'costate layer with suspicion'], info.spreadPct);
end

% switching function + certified (dual-S) switch times, sub-node by secant
c  = p.c;
S  = 1 - sqrt(sum(lam(4:6,:).^2, 1)).*c./X(7,:) - lam(7,:);
cr = find(diff(sign(S)) ~= 0);
tauSwitchS = zeros(1, numel(cr));
for q = 1:numel(cr)
    kk = cr(q);
    tauSwitchS(q) = tauN(kk) + (0-S(kk))*(tauN(kk+1)-tauN(kk))/(S(kk+1)-S(kk));
end
% raw throttle crossings (overcounts on shallow dips; cross-check only)
sT = U(4, :);
crT = find(diff(sign(sT - 0.5)) ~= 0);
tauSwitchThrottle = zeros(1, numel(crT));
for q = 1:numel(crT)
    kk = crT(q);
    tauSwitchThrottle(q) = tauN(kk) + (0.5-sT(kk))*(tauN(kk+1)-tauN(kk))/(sT(kk+1)-sT(kk));
end

% ---- data-product structs ----------------------------------------------------
mesh = struct('sigma', sigma(:).', 'tau', tauN, 'tauf0', tauf0, ...
              't', X(8,:), 'tDays', X(8,:)*p.tStar/86400, ...
              'pSund', cfg.pSund, 'nN', nN);
traj = struct('r', X(1:3,:), 'v', X(4:6,:), 'm', X(7,:), 'X', X);
ctrl = struct('alpha', U(1:3,:), 's', sT, ...
              'nSwitchS', numel(tauSwitchS), 'tauSwitchS', tauSwitchS, ...
              'nSwitchThrottle', numel(tauSwitchThrottle), ...
              'tauSwitchThrottle', tauSwitchThrottle);
costate = struct('lam', lam, 'S', S, 'beta', info.beta, ...
                 'spreadPct', info.spreadPct, 'mode', 'd', 'M', opts.M, ...
                 'accuracy', 'O(h) mesh-accuracy (~1%); see ifs/RESULTS_RUNG01_RUNG2.md Rung A');
pmp  = struct( ...
    'lamM_end',    lam(7, end), ...                       % transversality, want 0
    'termPosErr',  norm(X(1:3,end) - rvfC(1:3)), ...      % rendezvous residual (ND)
    'termVelErr',  norm(X(4:6,end) - rvfC(4:6)), ...
    'termTimeErr', X(8,end) - factor*cfg.tfMin, ...       % fixed-t_f residual (ND)
    'SsignAgree',  100*mean((S < 0) == (sT > 0.5)), ...   % sign law vs throttle
    'primerAlignDeg', get_field(out, 'primerAlignDeg', NaN));
scal = struct('factor', factor, 'tf', tf, 'tf_days', tf*p.tStar/86400, ...
              'dV', p.c*log(1/X(7,end))*p.lStar/p.tStar, ...
              'prop_kg', p.m0kg*(1 - X(7,end)), 'mf', X(7,end), ...
              'maxDefect', get_field(out, 'maxDefect', NaN), ...
              'switches', get_field(out, 'switches', NaN), ...
              'edge', get_field(out, 'edge', NaN));
const = p;  const.tfMin = cfg.tfMin;
provenance = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
    'source', char(solFile), 'gitHash', git_hash(fileparts(mfilename('fullpath'))), ...
    'dualMap', sprintf('sms_seed_duals mode d, M=%d, beta=%.6g', opts.M, info.beta), ...
    'pipeline', 'PSR/run_psr.m');

% ---- write, seed-compatible layer + products together -------------------------
fTag = strrep(sprintf('%.3f', factor), '.', 'p');
eTag = strrep(sprintf('%g', opts.epsMin), '.', 'p');   % 0 -> '0', 0.001 -> '0p001'
insertion = insMeta.label; %#ok<NASGU>
% tag inserted BEFORE _minEps (not appended at the end) so the existing
% ztl/p0a_graze_margin.m glob ('psr_data_tf*_minEps0.mat') still matches.
dataFile = fullfile(dataDir, ...
    sprintf('psr_data_tf%s_sw%d_%s_minEps%s.mat', fTag, ctrl.nSwitchS, insMeta.label, eTag));
provenance.epsMin = opts.epsMin;
save(dataFile, 'out', 'sigma', 'tauf0', 'rv0', 'rvf', 'factor', 'insertion', ...
     'mesh', 'traj', 'ctrl', 'costate', 'pmp', 'scal', 'const', 'provenance');

if ~opts.quiet
    fprintf(['psr_export_data: wrote %s\n' ...
             '  k=%d (dual-S; throttle raw %d)  beta=%.4g (spread %.2f%%)\n' ...
             '  transversality lamM(sigf)=%.2e  term rv err=(%.2e, %.2e)  S-sign %.1f%%\n'], ...
            dataFile, ctrl.nSwitchS, ctrl.nSwitchThrottle, info.beta, ...
            info.spreadPct, pmp.lamM_end, pmp.termPosErr, pmp.termVelErr, pmp.SsignAgree);
end
end

% -----------------------------------------------------------------------------
function v = get_field(s, f, dflt)
% Field with default (raw casadi structs carry different subsets).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

% -----------------------------------------------------------------------------
function h = git_hash(here)
% Short git hash for provenance (empty on failure -- never blocks the export).
[rc, s] = system(sprintf('cd "%s" && git rev-parse --short HEAD 2>/dev/null', here));
if rc == 0, h = strtrim(s); else, h = ''; end
end
