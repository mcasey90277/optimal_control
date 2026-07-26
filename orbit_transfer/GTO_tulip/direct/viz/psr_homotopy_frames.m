function framesFile = psr_homotopy_frames(factor, opts)
% PSR_HOMOTOPY_FRAMES  Capture the throttle + metrics at every epsilon of the
% energy->fuel homotopy, for the min-fuel movie.
%
% Replays MINFUEL_AT_TF's 'energy'-seed recipe (tight re-clean at eps=1, then a
% densified Bertrand-Epenoy energy->fuel sharpen down to eps=0) BUT records the
% converged throttle s(tau), the carried time state, and the summary metrics at
% EACH epsilon step -- so a movie can show the control sharpening from the smooth
% min-energy ramp (eps=1) into the bang-bang min-fuel solution (eps=0). The
% Sundman node set sigma is FIXED across every epsilon (the solver takes the same
% sigma each call), so the throttle profiles are directly comparable frame to
% frame. Each solve is warm-started from the previous epsilon, exactly as the
% production driver, so the captured path is the genuine homotopy path.
%
% INPUTS:
%   factor - t_f / t_f_min [scalar] (1.80 -> the f1800 case)
%   opts   - (optional) struct:
%            epsList  - capture schedule, 1 -> 0 [1xK]  (default: densified)
%            maxIter  - IPOPT cap [scalar]              (default cfg.maxIter)
%            outFile  - frames .mat path [char]         (default results/...)
%
% OUTPUTS:
%   framesFile - path to the saved cache .mat, containing:
%     frames  - 1xK struct array (only converged steps kept): fields
%               eps, s [1x(N+1)], sigma [1x(N+1)], tDays [1x(N+1)], dV, prop_kg,
%               switches, edge, mf, defect, primerDeg, ipoptStatus
%     meta    - factor, tf_days, tfMin, rv0, rvf, insertionLabel, date, epsList
%
% REFERENCES:
%   [1] PSR/lib/minfuel_at_tf.m (the recipe this mirrors, capturing every step).
%   [2] PSR/lib/casadi_minfuel_sundman.m (the per-epsilon solver).
%   [3] LOW_THRUST_MINFUEL_CAMPAIGN.md ("SOLVED": Sundman + energy->fuel homotopy).

here = fileparts(mfilename('fullpath'));
addpath(here);  setup_paths();
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

if nargin < 2, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
% Densified capture schedule: finer than cfg.schedSharpen up top (where the
% control changes the most) so the movie leans less on interpolation, and finer
% near 0 (coarse steps there break the basin -- campaign lesson). Ends at 0.
epsList = d('epsList', [1.0 0.85 0.7 0.6 0.5 0.42 0.35 0.28 0.22 0.18 0.14 ...
                        0.11 0.09 0.07 0.055 0.045 0.035 0.028 0.022 0.017 ...
                        0.013 0.010 0.007 0.005 0.003 0.002 0.001 0]);
maxIter = d('maxIter', cfg.maxIter);
tag     = sprintf('f%04d', round(1000*factor));
outFile = d('outFile', fullfile(here, 'results', ['psr_homotopy_frames_' tag '.mat']));

% --- seed: the min-energy backbone solution at this factor -------------------
seedFile = fullfile(cfg.dirs.energy, cfg.fname('energy', factor));
assert(isfile(seedFile), 'psr_homotopy_frames:noSeed', ...
    'no energy backbone seed %s -- run the energy sweep first', seedFile);
E = load(seedFile);                         % sigma, rv0, rvf, tauf0, X, U
sigma = E.sigma;  rv0 = E.rv0;  rvf = E.rvf;  tauf0 = E.tauf0;
Xk = E.X;  Uk = E.U;
tf = factor * cfg.tfMin;

fprintf('\n=== PSR_HOMOTOPY_FRAMES factor %.3f  t_f=%.4f ND (%.2f d) ===\n', ...
        factor, tf, tf*p.tStar/86400);
fprintf('    seed: %s   |   %d epsilon steps 1 -> 0\n', seedFile, numel(epsList));

frames = struct('eps',{},'s',{},'sigma',{},'tDays',{},'dV',{},'prop_kg',{}, ...
                'switches',{},'edge',{},'mf',{},'defect',{},'primerDeg',{}, ...
                'ipoptStatus',{});
for ke = 1:numel(epsList)
    e = epsList(ke);
    % warmTight=true throughout: fixed t_f, re-solving from a tight warm start
    % (no continuation MOVE), same as minfuel_at_tf for an energy seed.
    o = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                               Xk, Uk, tauf0, cfg.pSund, maxIter, e, true);
    ok = o.success && o.maxDefect < 1e-6;
    fprintf('  eps=%6.4f  ok=%d  defect=%.2g  sw=%3d  edge=%5.1f%%  dV=%.4f\n', ...
            e, ok, o.maxDefect, o.switches, 100*o.edge, ...
            p.c*log(1/o.mf)*p.lStar/p.tStar);
    if ok
        fr = struct();
        fr.eps       = e;
        fr.s         = o.U(4,:);
        fr.sigma     = sigma(:).';
        fr.tDays     = o.X(8,:) * p.tStar/86400;
        fr.dV        = p.c*log(1/o.mf) * p.lStar/p.tStar;
        fr.prop_kg   = p.m0kg*(1 - o.mf);
        fr.switches  = o.switches;
        fr.edge      = o.edge;
        fr.mf        = o.mf;
        fr.defect    = o.maxDefect;
        fr.primerDeg = o.primerAlignDeg;
        fr.ipoptStatus = o.ipoptStatus;
        frames(end+1) = fr; %#ok<AGROW>
        Xk = o.X;  Uk = o.U;             % advance the warm start ONLY on success
    end
end

assert(~isempty(frames), 'psr_homotopy_frames:noConverged', ...
    'no epsilon step converged tight -- nothing to animate');

meta = struct('factor', factor, 'tf_days', tf*p.tStar/86400, 'tfMin', cfg.tfMin, ...
    'rv0', rv0, 'rvf', rvf, 'seedFile', seedFile, 'epsList', epsList, ...
    'nCaptured', numel(frames), ...
    'date', char(datetime('now','Format','yyyy-MM-dd HH:mm'))); %#ok<NASGU>

if ~exist(fullfile(here,'results'),'dir'), mkdir(fullfile(here,'results')); end
save(outFile, 'frames', 'meta');
fprintf('=== captured %d/%d epsilon frames -> %s ===\n', ...
        numel(frames), numel(epsList), outFile);
framesFile = outFile;
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
