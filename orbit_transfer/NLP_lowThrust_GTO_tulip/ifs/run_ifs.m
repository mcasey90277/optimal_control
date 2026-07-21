% RUN_IFS  Entry-level driver for a full IFS (Indirect Finishing Solve) run.
%
% IFS is "point 3" of the direct<->indirect roadmap: take a good DIRECT / PSR
% bang-bang min-fuel solution and hand it to an ACTUAL INDIRECT solver that
% holds the switch STRUCTURE fixed and places every switch exactly at S(tau)=0,
% producing exact sub-mesh switch times, exact costates, and a continuous-time
% first-order PMP certificate. It is the sibling of PSR (the direct-side
% mesh-refinement stage); IFS goes the rest of the way to a costate-carrying
% indirect solution. One run of this script does:
%
%   1. PARAMETERS   - set t_f, the seed direct solution, and solver knobs here.
%   2. BUILD SEED   - build the IFS unknown vector Z (initial costates + node
%                     states + switch times) from a direct/PSR .mat: costates
%                     from the KKT-dual map ('dual') or an adjoint sweep/smoother
%                     ('adjoint'); switch times from the dual-S crossings; arc
%                     throttles from the direct sign pattern.
%   3. INDIRECT     - drive the square multiple-shooting PMP residual to zero
%      SOLVE          with a scaled, rank-revealing truncated-SVD Gauss-Newton
%                     step (ifs_solve2), complex-step Jacobian, no epsilon layer.
%   4. CERTIFY      - first-order PMP certificate (ifs_certify): S=0 at each
%                     switch, the bang-bang sign law on every arc interior,
%                     terminal residual, rendezvous transversality lamM(tau_f)=0.
%   5. RECONSTRUCT  - integrate the IFS arcs (ifs_reconstruct) into a seed-layout
%      + DATA EXPORT  trajectory + costates and save ALL products (mesh, state,
%                     control, costates, switching function, switch times,
%                     certificate, provenance) to ../IFS_data/.
%   6. MOVIE        - control movie (reuses PSR/psr_movie): rotating-frame
%                     transfer colored burn/coast, primer thrust arrows,
%                     throttle strip, running Delta-V.
%
% *** HONEST STATUS (2026-07-12): IFS is OPEN. *** The machinery is validated
% (unit tests green) and the min-time k=0 anchor converges, but the full
% multi-switch solve does NOT converge from a cold direct/PSR seed: ifs_solve2
% DESCENDS the residual (e.g. 1.96 -> ~0.4) then FLOORS on the small convergence
% basin of 40-rev shooting -- the same conditioning wall that this whole project
% keeps hitting. So stage 3 will report a residual floor, not tolR, and stage 4
% yields at best a PARTIAL certificate. This driver runs the attempt end to end
% and reports honestly. See README.md, RESULTS_RUNG01_RUNG2.md, PLAN_RUNG_B.md.
%
% REFERENCES:
%   [1] README.md, RESULTS.md, RESULTS_RUNG01_RUNG2.md (this folder: the arc)
%   [2] PLAN_OF_ATTACK_2.md / PLAN_RUNG_B.md (next levers)
%   [3] Zhang, Topputo, Bernelli-Zazzera, Zhao, JGCD 38(8), 2015 (indirect min-fuel)

%% ------------------------------------------------------------------------
%% 0. Paths and output folder
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
addpath(fullfile(here, '..', 'PSR'));          % psr_movie for stage 6
addpath(fullfile(here, '..', 'PSR', 'lib'));   % cr3bp_lt_params etc.
resDir  = fullfile(here, 'results');           % pipeline intermediates
dataDir = fullfile(here, '..', 'IFS_data');    % exported data products (stage 5)
if ~exist(resDir, 'dir'), mkdir(resDir); end

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------

% ---- transfer time --------------------------------------------------------
% t_f = factor * tfMin, tfMin = 6.2906939607 ND. 1.12x is the canonical gate
% (k=10 certified switches, the smallest many-switch system on disk).
factor = 1.25;                   % t_f / t_f_min

% ---- SEED (the DIRECT / PSR bang-bang solution to finish indirectly) -------
% Any direct-solution .mat with out.X/out.U, out.lamDef (KKT duals), sigma,
% tauf0, rv0, rvf works -- legacy campaign solutions, or a PSR result
% (../PSR_data/psr_data_*.mat or ../PSR/results/psr_refined_*.mat). The IFS
% seed is built FROM this bang-bang solution (its switch structure + duals).
seedFile = fullfile(here, '..', 'sundman_minfuel', 'results', 'minfuel', ...
                    'legacy_ms_f1120.mat');    % the 1.12x gate

% ---- seed method ----------------------------------------------------------
% 'dual'    (default) costates from the adjudicated mode-'d' KKT dual->costate
%           map (ifs_seed). The cold seed the campaign characterized.
% 'adjoint' costates from a backward adjoint sweep / smoother (ifs_seed_adjoint,
%           Rung A). Built and FALSIFIED as a cold seed (the 40-rev adjoint flow
%           amplifies ~1e12); kept for the record. See RESULTS_RUNG01_RUNG2.md.
seedMethod = 'dual';             % 'dual' | 'adjoint'

% ---- switch-time parameterization -----------------------------------------
tauParam = 'sigmoid';            % 'sigmoid' (monotone stick-breaking) | 'direct'

% ---- indirect-solve knobs (ifs_solve2) ------------------------------------
solveOpts = struct( ...
    'tolR',     1e-8, ...        % success threshold on ||R||_2 (OPEN: not reached)
    'maxIter',  200, ...         % Gauss-Newton iterations
    'relTrunc', 1e-2, ...        % starting SVD truncation ratio
    'verbose',  true);

% ---- movie ----------------------------------------------------------------
movieMode = 'movie';             % 'preview' (3 stills, fast) | 'movie' | 'none'

% ---- stage rerun control --------------------------------------------------
rerunSolve = false;              % skip stage 3 if its output exists

% ---- derived file names (canonical; do not edit) --------------------------
tag       = sprintf('f%04d_%s', round(1000*factor), seedMethod);
solveFile = fullfile(resDir, ['ifs_solve_' tag '.mat']);
reconFile = fullfile(resDir, ['ifs_recon_' tag '.mat']);   % seed layout (movie)

fprintf('\n=== IFS PIPELINE: factor=%.3f, seed=%s, method=%s ===\n', ...
        factor, seedFile, seedMethod);

%% ------------------------------------------------------------------------
%% 2. BUILD IFS SEED  (costates + node states + switch times from the direct sol)
%% ------------------------------------------------------------------------
fprintf('\n[stage 2] BUILD IFS SEED (%s)...\n', seedMethod);
switch seedMethod
    case 'dual'
        [Z0, prob, meta] = ifs_seed(seedFile, struct('mode','full','tauParam',tauParam));
    case 'adjoint'
        [Z0, prob, meta] = ifs_seed_adjoint(seedFile, struct('method','smooth','tauParam',tauParam));
    otherwise
        error('run_ifs:seedMethod', 'unknown seedMethod %s', seedMethod);
end
prob.factor = factor;                                 % for reconstruction/export
fprintf('[stage 2] k=%d switches  seed ||R||=%.4e\n', prob.k, meta.seedResNorm);

%% ------------------------------------------------------------------------
%% 3. INDIRECT SOLVE  (truncated-SVD Gauss-Newton on the PMP residual)
%% ------------------------------------------------------------------------
% ifs_solve2 descends the residual with a scaled, rank-revealing truncated-SVD
% step (drops near-null directions below relTrunc*sigma_max, so the step moves
% along the well-determined directions instead of blowing up along the weakly
% determined lambda_r0 direction), a Levenberg fallback, an alpha-floor line
% search, and adaptive truncation continuation. NO epsilon / smoothing layer.
% OPEN: from a cold seed this DESCENDS then FLOORS (the 40-rev shooting basin).
if isfile(solveFile) && ~rerunSolve
    fprintf('\n[stage 3] solve exists (%s) -- skipping. Set rerunSolve=true to redo.\n', solveFile);
    L = load(solveFile);  out = L.out;
else
    fprintf('\n[stage 3] INDIRECT SOLVE (ifs_solve2, maxIter=%d)...\n', solveOpts.maxIter);
    out = ifs_solve2(Z0, prob, solveOpts);
    save(solveFile, 'out', 'prob', 'meta');
end
fprintf(['[stage 3] seed ||R||=%.3e -> ||R||=%.3e  (tolR=%.0e)  success=%d  ' ...
         'iters=%d  flag=%d\n'], out.seedResNorm, out.resNorm, solveOpts.tolR, ...
         out.success, out.iterations, out.flag);
if ~out.success
    fprintf(['[stage 3] NOTE: did not reach tolR -- IFS is OPEN; ||R|| floored on the ' ...
             '40-rev cold-seed shooting basin (expected). Proceeding with the best iterate.\n']);
end

%% ------------------------------------------------------------------------
%% 4. CERTIFY  (first-order PMP certificate -- extremality)
%% ------------------------------------------------------------------------
% ifs_certify checks the continuous-time first-order PMP conditions of the best
% iterate: S=0 at each switch node, the bang-bang sign law on every arc interior
% (S<0 burn / S>0 coast, bounded from 0 => no singular arc), the terminal
% rendezvous residual, and the free-final-mass transversality lamM(tau_f)=0. It
% also reports structure diagnostics (vanishing arc => spurious switch; in-arc
% sign violation => missing switch). With ||R|| not at tolR the certificate is
% PARTIAL (residual-limited) -- reported honestly, not issued as PASS.
fprintf('\n[stage 4] FIRST-ORDER PMP CERTIFICATE...\n');
cert = [];
try
    cert = ifs_certify(out.Z, prob, meta);
    fprintf('[stage 4] cert.ok=%d  max|S(switch)|=%.2e  signViol=%.2e  |lamM(tf)|=%.2e  termRes=%.2e\n', ...
        cert.ok, max(abs(cert.Sswitch)), cert.signViol, abs(cert.lamMend), cert.termResNorm);
catch ME
    fprintf('[stage 4] ifs_certify errored (%s) -- continuing without a certificate.\n', ME.message);
end

%% ------------------------------------------------------------------------
%% 5. RECONSTRUCT + DATA EXPORT  (trajectory + costates -> ../IFS_data/)
%% ------------------------------------------------------------------------
% The IFS unknown Z is costates + nodes + switch times, not a sampled path.
% ifs_reconstruct propagates the k+1 arcs into a seed-layout trajectory (state,
% primer-vector control, costates) so the result can be analyzed and animated
% with the same tooling as a direct/PSR solution. Saved to ../IFS_data/ with the
% seed layout at top level (out, sigma, tauf0, rv0, rvf, factor) PLUS the
% costates, switching function, solve residual, and certificate.
% NB: this runs on the BEST ITERATE regardless of convergence -- even when the
% solve floored (IFS OPEN), we still write the .mat and animate the control we
% got (the reconstructed trajectory will miss the terminal rendezvous by the
% residual, which is the honest picture of a non-converged indirect finish).
fprintf('\n[stage 5] RECONSTRUCT + DATA EXPORT...\n');
recon = ifs_reconstruct(out.Z, prob, 4000);
% movie-consumable seed-layout file (psr_movie reads out.X/out.U, sigma, rvf)
mov = struct('out', recon.out, 'sigma', recon.sigma, 'tauf0', recon.tauf0, ...
             'rv0', recon.rv0, 'rvf', recon.rvf, 'factor', recon.factor);
save(reconFile, '-struct', 'mov');
% full data product
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
provenance = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
    'seedFile', char(seedFile), 'seedMethod', seedMethod, 'pipeline', 'ifs/run_ifs.m', ...
    'resNorm', out.resNorm, 'success', out.success);
dataFile = fullfile(dataDir, sprintf('ifs_data_%s.mat', tag));
ifsRecon = recon;  ifsSolve = out;  ifsCert = cert;  %#ok<NASGU>
save(dataFile, 'ifsRecon', 'ifsSolve', 'ifsCert', 'prob', 'meta', 'provenance');
fprintf('[stage 5] switches(recon)=%d  reconFile=%s\n[stage 5] data product: %s\n', ...
        recon.out.switches, reconFile, dataFile);

%% ------------------------------------------------------------------------
%% 6. CONTROL MOVIE  (transfer + control law, synced -- reuses PSR/psr_movie)
%% ------------------------------------------------------------------------
if ~strcmp(movieMode, 'none')
    fprintf('\n[stage 6] CONTROL MOVIE (%s)...\n', movieMode);
    titleStr = sprintf('IFS min-fuel GTO\\rightarrowtulip, t_f = %.2fx min-time (%d-switch, indirect)', ...
                       factor, recon.out.switches);
    psr_movie(reconFile, fullfile(resDir, ['ifs_movie_' tag]), titleStr, movieMode);
end

fprintf('\n=== IFS PIPELINE DONE (factor %.3f). Intermediates: %s  Data products: %s ===\n', ...
        factor, resDir, dataFile);
