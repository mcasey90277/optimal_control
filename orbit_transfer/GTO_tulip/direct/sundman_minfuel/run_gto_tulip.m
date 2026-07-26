% RUN_GTO_TULIP  Front door: min-fuel GTO -> south-pole tulip transfer (CR3BP).
%
% One self-contained entry script (house run_gergaud / run_cr3bp_geo pattern):
% set the parameters in section 1, run, and get a solved transfer + a saved
% data product + plots, all under a run name you choose.
%
% The chain is the certified campaign pipeline, and this script does not
% re-implement it -- it drives minfuel_at_tf, THE canonical per-t_f driver, so
% every run inherits that function's seed-fingerprint checks, schedule policy
% from minfuel_config, certification gate and provenance stamping:
%
%   energy backbone at this factor  ->  tight re-clean at eps=1
%                                   ->  energy->fuel homotopy down to epsMin
%
%   epsMin = 1   -> minimum ENERGY solution (smooth throttle, no sharpening)
%   epsMin = 0   -> minimum FUEL solution (bang-bang; the campaign objective)
%   0<epsMin<1   -> eps-optimal (partially sharpened) solution
%
% CERTIFICATION IS NOT ASSUMED. minfuel_at_tf certifies a run only if the
% homotopy reaches the REQUESTED endpoint eps with a tight solve, and it
% refuses to write an uncertified result (a loose iterate must never become a
% neighbour seed). This script reports that verdict prominently rather than
% burying it.
%
% THE TWO KNOBS THAT LOOK ADJUSTABLE BUT ARE NOT, AND WHY:
%
%   t_f factor. Reliable band is ~1.12-1.95, which is exactly where energy
%   backbones exist (results/energy/energy_f1120..f1950). The hard band
%   1.01-1.11 approaching min-time is an OPEN campaign problem, not a setting:
%   the energy backbone itself will not converge there. Ask for a factor with
%   no backbone and this script says so, and names the generator.
%
%   thrust. Fixed at the campaign's 25 mN. Off-nominal thrust is the OPEN
%   thrust-ladder problem: the 20 mN pilot (pilot_rung_20mN.m) was an honest
%   FAILURE against a fixed-tau_f topology wall -- a chained 25 mN winding
%   cannot grow to the rev count a lower-thrust rung needs, so the eps=1 energy
%   step will not close. This script refuses off-nominal thrust with that
%   explanation instead of running a solve that is expected to fail. See
%   process/LADDER_PREP_PILOT_FINDINGS.md and GTO_tulip/TODO.md.
%
% This is why the companion sweep script (run_tulip_front.sh) sweeps t_f and
% NOT thrust: for this campaign t_f is the dimension that works, and the
% Delta-V/t_f front is the campaign's mapped result. The earth and CR3BP
% campaigns sweep thrust because for them thrust is the dimension that works.
%
% OUTPUTS (files, under results/minfuel/):
%   <runName>.mat      - out (solver struct + .certified .dV .prop_kg .meta
%                        provenance), sigma, tauf0, rv0, rvf, factor, fp
%   <runName>.png      - throttle history + rotating-frame trajectory
%
% NO MOVIE MODE, deliberately: this campaign's renderers are either a script
% with hard-coded absolute paths (direct/movie/gen_movie_data.m) or PSR-owned
% (PSR/psr_movie.m), and reaching into PSR would breach the isolation its
% lib/README.md defines. Render separately from those tools.
%
% REFERENCES:
%   [1] minfuel_at_tf.m (the canonical per-t_f driver this script drives).
%   [2] process/LOW_THRUST_MINFUEL_CAMPAIGN.md (campaign record).
%   [3] process/LADDER_PREP_PILOT_FINDINGS.md (the thrust topology wall).
%   [4] ../../../earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m (pattern).

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------
factor    = 1.15;        % t_f / t_f^min. Backbones exist 1.12-1.95; 1.15 is
                         %   the certified reference (25 switches, 2.264 kg)
epsMin    = 0;           % homotopy endpoint: 1 = min-ENERGY, 0 = min-FUEL
seedKind  = 'energy';    % 'energy' (backbone at this factor) | 'neighbor'
seedFactor= NaN;         % required for 'neighbor': factor of the seed solution
thrustN   = 0.025;       % max thrust [N]. Off-nominal is REFUSED -- see header
runName   = '';          % '' -> canonical minfuel_f####_<branch>.mat name
maxIter   = 1500;        % IPOPT cap per homotopy step
doPlot    = true;        % throttle + trajectory figure
saveResult= true;        % write the .mat (uncertified runs are never written)
rerun     = false;       % true -> ignore an existing result and solve cold

% Programmatic override hook (sweep driver): if the caller's workspace defines
% TULIP_OVERRIDES (struct), its fields replace the defaults above. Absent ->
% byte-identical interactive behaviour. Note this script does NOT clear the
% workspace, precisely so this hook can work -- see run_tulip_front.sh.
if exist('TULIP_OVERRIDES','var') && isstruct(TULIP_OVERRIDES)
    ovf = fieldnames(TULIP_OVERRIDES);
    for ovk = 1:numel(ovf), eval([ovf{ovk} ' = TULIP_OVERRIDES.(ovf{ovk});']); end
    clear ovf ovk
end

%% ------------------------------------------------------------------------
%% 2. SETUP + PRECONDITIONS  (fail loudly and specifically, before solving)
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
cfg  = minfuel_config();
p    = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);   % ND constants + tStar

% Thrust: refuse off-nominal rather than run a solve known to hit the wall.
if abs(thrustN - cfg.thrustN) > 1e-12
    error('run_gto_tulip:offNominalThrust', ...
        ['thrustN = %.4g N requested, but this campaign is certified only at ' ...
         '%.4g N.\n' ...
         'Off-nominal thrust is the OPEN thrust-ladder problem, not a setting: ' ...
         'the 20 mN pilot (pilot_rung_20mN.m) failed against a fixed-tau_f ' ...
         'topology wall -- a chained %.4g N winding cannot grow to the rev count ' ...
         'a lower-thrust rung needs, so the eps=1 energy step does not close.\n' ...
         'To work ON that problem, run pilot_rung_20mN.m directly and read ' ...
         'process/LADDER_PREP_PILOT_FINDINGS.md. This front door will not ' ...
         'pretend the rung solves.'], thrustN, cfg.thrustN, cfg.thrustN);
end

% epsMin: build the schedule by truncating the campaign's sharpening schedule.
% Truncation (rather than a fresh schedule) keeps every run on the same
% certified step sequence, so a partial run is a prefix of the full one.
validateattributes(epsMin, {'numeric'}, {'scalar','>=',0,'<=',1}, mfilename, 'epsMin');
if epsMin >= 1
    sched = 1;                                  % energy only: one tight eps=1 solve
else
    sched = cfg.schedSharpen(cfg.schedSharpen > epsMin);
    sched = [sched, epsMin];                    % always END at the request
end

% Seed availability: for an energy seed, say plainly if no backbone exists.
if strcmpi(seedKind, 'energy')
    ef = fullfile(cfg.dirs.energy, cfg.fname('energy', factor));
    if ~isfile(ef)
        have = dir(fullfile(cfg.dirs.energy, 'energy_f*.mat'));
        haveF = sort(cellfun(@(n) cfg.fparse(n), {have.name}));
        error('run_gto_tulip:noBackbone', ...
            ['no energy backbone at factor %.3f (looked for %s).\n' ...
             'Backbones on disk: %s\n' ...
             'The band below ~1.12 is an OPEN problem -- the energy backbone ' ...
             'itself will not converge approaching min-time. Inside the band, ' ...
             'generate one with gen_tulip_energy_2p, or use ' ...
             'seedKind = ''neighbor'' with seedFactor set to a solved factor.'], ...
            factor, ef, strjoin(cellstr(compose('%.3f', haveF(:))).', ', '));
    end
end

% Resolve the artifact name once, so resume and save agree.
branch = 'en';  if strcmpi(seedKind, 'neighbor'), branch = 'nb'; end
if isempty(runName)
    runName = strrep(cfg.fname('minfuel', factor), '.mat', ['_' branch]);
end
if ~exist(cfg.dirs.minfuel, 'dir'), mkdir(cfg.dirs.minfuel); end
outFile = fullfile(cfg.dirs.minfuel, [runName '.mat']);
pngFile = fullfile(cfg.dirs.minfuel, [runName '.png']);

fprintf(['\n=== RUN_GTO_TULIP =======================================\n' ...
         '  factor   : %.3f  (t_f = %.4f ND = %.2f days)\n' ...
         '  epsMin   : %.4g   -> %s\n' ...
         '  seed     : %s\n' ...
         '  run name : %s\n' ...
         '=========================================================\n'], ...
        factor, factor*cfg.tfMin, factor*cfg.tfMin*p.tStar/86400, epsMin, ...
        epsLabel(epsMin), seedKind, runName);

%% ------------------------------------------------------------------------
%% 3. SOLVE  (or resume)
%% ------------------------------------------------------------------------
resumed = false;
if ~rerun && isfile(outFile)
    S = load(outFile);
    out = S.out;  sigma = S.sigma;  rv0 = S.rv0;  rvf = S.rvf;
    resumed = true;
    fprintf('RESUMED from %s (set rerun = true to solve cold)\n', outFile);
else
    args = {'seed', seedKind, 'sched', sched, 'maxIter', maxIter, ...
            'branch', branch, 'outFile', outFile, 'save', saveResult};
    if strcmpi(seedKind, 'neighbor'), args = [args, {'seedFactor', seedFactor}]; end
    out = minfuel_at_tf(factor, args{:});
    % minfuel_at_tf returns the solution but keeps its inputs internal; reload
    % the endpoints/mesh it saved so the plots below describe the SAME artifact.
    if out.certified && saveResult && isfile(outFile)
        S = load(outFile);  sigma = S.sigma;  rv0 = S.rv0;  rvf = S.rvf;
    else
        E = load(fullfile(cfg.dirs.energy, cfg.fname('energy', factor)));
        sigma = E.sigma;  rv0 = E.rv0;  rvf = E.rvf;
    end
end

%% ------------------------------------------------------------------------
%% 4. REPORT
%% ------------------------------------------------------------------------
fprintf('\n=========================== RESULT ============================\n');
fprintf('  transfer time    : %.4f ND  (%.2f days,  %.3f x min-time)\n', ...
        out.tf, out.tf_days, out.factor);
fprintf('  delta-V          : %.4f km/s\n', out.dV);
fprintf('  propellant       : %.4f kg   (final mass fraction %.6f)\n', ...
        out.prop_kg, out.mf);
fprintf('  throttle switches: %d   (bang-bang fraction %.1f%%)\n', ...
        out.switches, 100*out.edge);
fprintf('  max defect       : %.2e   primer alignment: %.3f deg\n', ...
        out.maxDefect, out.primerAlignDeg);
fprintf('  eps reached      : %.4g  (requested %.4g)\n', out.epsReached, epsMin);
fprintf('  CERTIFIED        : %d%s\n', out.certified, ...
        repmat('   <-- tight solve at the requested eps', 1, out.certified));
if resumed
    fprintf('  source           : resumed artifact %s\n', outFile);
elseif out.certified && saveResult
    fprintf('  saved to         : %s\n', outFile);
end
fprintf('===============================================================\n');

% --- cross-check against a known certified solution at the SAME t_f ---------
% MEASURED 2026-07-26, and the reason this check exists: at factor 1.150 the
% energy-backbone route (what minfuel_at_tf drives, hence what this front door
% drives) and the run_certified_minfuel chain both converge machine-tight to 25
% switches at the SAME t_f -- and land on DIFFERENT local optima:
%
%   run_certified_minfuel chain : m_f 0.849066, dV 3.3696 km/s, 2.2640 kg
%   energy-backbone route       : m_f 0.847086, dV 3.4176 km/s, 2.2937 kg
%                                 ------------------------------------------
%                                 +1.43% dV, +0.0297 kg  (WORSE)
%
% Both are valid extremals; min-fuel here has closely spaced optima and the SEED
% ROUTE decides which one you land on. Without this check a front-door run at the
% flagship factor would quietly report numbers 1.4% off the campaign's headline
% result and look completely healthy. Higher m_f is better (minimum fuel means
% maximum final mass), so a negative delta below means this run found the worse
% basin -- not a failure, but not the number to quote.
if out.certified && isfile(fullfile(here, 'sundman_minfuel_certified.mat'))
    K = load(fullfile(here, 'sundman_minfuel_certified.mat'));
    if isfield(K, 'out') && isfield(K.out, 'X') && ...
            abs(K.out.X(8,end) - out.X(8,end)) < 1e-9      % same t_f
        dmf = out.mf - K.out.mf;
        fprintf(['\n  cross-check vs sundman_minfuel_certified.mat (same t_f):\n' ...
                 '    m_f  this run %.6f  vs certified %.6f   (%+.6f)\n' ...
                 '    switches   %3d           vs %3d\n'], ...
                out.mf, K.out.mf, dmf, out.switches, K.out.switches);
        if dmf < -1e-6
            fprintf(['    -> this run found a WORSE (lower-mass) optimum at the same t_f.\n' ...
                     '       Both are certified extremals; the seed route decides which.\n' ...
                     '       Quote run_certified_minfuel''s numbers as the campaign result.\n']);
        elseif dmf > 1e-6
            fprintf(['    -> this run found a BETTER (higher-mass) optimum than the banked\n' ...
                     '       certified solution. Worth investigating before quoting either.\n']);
        else
            fprintf('    -> same optimum.\n');
        end
    end
end

if ~out.certified
    warning('run_gto_tulip:uncertified', ...
        ['this run is NOT certified: the homotopy reached eps=%.4g, not the ' ...
         'requested %.4g. It was deliberately NOT saved. Do not use it as a ' ...
         'neighbour seed and do not quote its numbers as a campaign result.'], ...
        out.epsReached, epsMin);
end

%% ------------------------------------------------------------------------
%% 5. PLOTS
%% ------------------------------------------------------------------------
if doPlot
    tDays = out.X(8,:) * p.tStar / 86400;      % physical time per node [days]
    thr   = out.U(4,:);                         % throttle history

    fig = figure('Name', sprintf('GTO->tulip  f=%.3f  eps=%.4g', factor, epsMin), ...
                 'Color', 'w', 'Position', [100 100 1100 450]);

    subplot(1,2,1);
    plot(tDays, thr, '-', 'LineWidth', 1.2);  grid on;
    xlabel('time [days]');  ylabel('throttle s');  ylim([-0.05 1.05]);
    title(sprintf('throttle: %d switches, %.4f km/s', out.switches, out.dV));

    subplot(1,2,2);
    plot(out.X(1,:), out.X(2,:), '-', 'LineWidth', 0.8);  hold on;  grid on;  axis equal;
    plot(-p.muStar,  0, 'b.', 'MarkerSize', 22);                                  % Earth
    plot(1-p.muStar, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [.5 .5 .5]);    % Moon
    plot(rv0(1), rv0(2), 'gs', 'MarkerSize', 8,  'MarkerFaceColor', 'g');         % GTO start
    plot(rvf(1), rvf(2), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');         % tulip target
    xlabel('x (rotating, ND)');  ylabel('y (rotating, ND)');
    legend({'trajectory','Earth','Moon','start','target'}, 'Location', 'best');
    title('GTO \rightarrow south-pole tulip');

    if saveResult
        exportgraphics(fig, pngFile, 'Resolution', 150);
        fprintf('  figure saved to  : %s\n', pngFile);
    end
end

% ---------------------------------------------------------------------------
function s = epsLabel(e)
% EPSLABEL  Human-readable name for a homotopy endpoint.
% INPUTS:  e - epsMin [scalar in [0,1]]
% OUTPUTS: s - description [char]
if e >= 1,      s = 'minimum ENERGY (smooth throttle)';
elseif e <= 0,  s = 'minimum FUEL (bang-bang)';
else,           s = 'eps-optimal (partially sharpened)';
end
end
