function cfg = minfuel_config()
% MINFUEL_CONFIG  Single source of truth for the min-fuel GTO->tulip campaign.
%
% Replaces the scattered magic numbers -- most importantly the minimum
% transfer time, which was previously recovered by the fragile hack
% `load(...,'tf'); tfMin = tf/1.15` in three different drivers.
%
% INPUTS:  none
% OUTPUTS:
%   cfg - struct:
%     .tfMin        - minimum transfer time [ND], from the certified indirect
%                     min-time solution (27.8845 d); anchor of all t_f factors
%     .pSund        - Sundman power in dt/dtau = r1^pSund [scalar]
%     .thrustN      - thrust [N];  .m0kg - wet mass [kg];  .ispS - Isp [s]
%     .schedSharpen - energy->fuel homotopy schedule for sharpening a TIGHT
%                     energy solution at fixed t_f (ends at exactly 0) [1xK]
%     .schedNeighbor- light re-sharpen schedule for continuing a bang-bang
%                     solution from a NEIGHBORING t_f [1xK]
%     .maxIter      - default IPOPT iteration cap [scalar]
%     .dirs         - result directories (.energy .minfuel .fronts .plots .logs)
%     .fname        - @(kind,factor) canonical result filename, e.g.
%                     fname('minfuel',1.20) -> 'minfuel_f1200.mat' (collision-
%                     free at 0.001 factor granularity, unlike the old %.2f)
%     .fparse       - @(name) factor recovered from a canonical filename
%
% REFERENCES:
%   [1] LOW_THRUST_MINFUEL_CAMPAIGN.md (campaign record; schedule provenance).
%   [2] HONEST_EVALUATION_DV_TF_FRONT.md sec "Honest problems" items 3, 5.

here = fileparts(mfilename('fullpath'));

% --- physics anchors --------------------------------------------------------
% tfMin: certified indirect min-time answer (lowThrust_GTO_tulip solver),
% cross-checked = X(8,end)/1.15 of sundman_minfuel_certified.mat to 1e-9.
cfg.tfMin   = 6.2906939607;      % [ND]  (= 27.8845 days)
cfg.pSund   = 1.5;               % Sundman power, dt/dtau = r1^1.5
cfg.thrustN = 0.025;  cfg.m0kg = 15;  cfg.ispS = 2100;

% --- homotopy schedules -----------------------------------------------------
% Sharpen an energy (eps=1) solution at fixed t_f down to pure fuel. NOTE the
% final exact 0: the legacy solve_tf_minfuel stopped at 0.001, leaving a tiny
% objective bias (flagged in HONEST_EVALUATION item on schedule inconsistency).
cfg.schedSharpen  = [0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0];
% Light re-sharpen when continuing an already-bang-bang neighbor to a new t_f.
cfg.schedNeighbor = [0.05 0.02 0.008 0.003 0.001 0];
cfg.maxIter = 1500;

% --- results layout ---------------------------------------------------------
r = fullfile(here, 'results');
cfg.dirs = struct('root',r, 'energy',fullfile(r,'energy'), ...
                  'minfuel',fullfile(r,'minfuel'), 'fronts',fullfile(r,'fronts'), ...
                  'plots',fullfile(r,'plots'), 'logs',fullfile(r,'logs'));

% --- canonical filenames ----------------------------------------------------
% factor encoded as milli-units: 1.20x -> f1200 (no %.2f collisions).
cfg.fname  = @(kind,factor) sprintf('%s_f%04d.mat', kind, round(1000*factor));
cfg.fparse = @(name) local_fparse(name);
end

% ---------------------------------------------------------------------------
function factor = local_fparse(name)
% Recover the t_f factor from a canonical result filename ('*_f####.mat').
tok = regexp(name, '_f(\d{4})\.mat$', 'tokens', 'once');
if isempty(tok), factor = NaN; else, factor = str2double(tok{1})/1000; end
end
