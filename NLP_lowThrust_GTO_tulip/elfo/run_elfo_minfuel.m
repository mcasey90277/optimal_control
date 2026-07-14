% RUN_ELFO_MINFUEL  Entry driver: generate a min-fuel GTO->ELFO transfer.
%
% The ELFO analog of PSR/run_psr.m. Edit the PARAMETERS section and run. Unlike
% the tulip PSR pipeline (fixed-t_f casadi_minfuel_sundman), the ELFO transfer
% runs on the FREE-t_f two-primary solver casadi_energy_freetf -- the lunar
% capture leg needs the two-primary Sundman clock kappa=(r1^-q+(r2/D)^-q)^(-p/q)
% that a single-primary (Earth-only) clock cannot resolve. One run does:
%
%   1. PARAMETERS - pick target tag, t_f (factor), the energy seed, epsMin.
%   2. SOLVE      - Bertrand-Epenoy energy->fuel homotopy eps 1 -> epsMin on the
%                   two-primary solver at PINNED t_f (gen_elfo_minfuel core).
%                   epsMin=0 -> bang-bang FUEL; epsMin>0 -> smooth regularized.
%   3. VERIFY     - independent endpoint + solver-free defect check.
%
% OUTPUT NAME carries the TARGET: results/minfuel_<target>_tf<fTag>_sw<k>_
% minEps<eTag>.mat  e.g.  minfuel_ELFO_tf1p200_sw34_minEps0.mat  (parallel to the
% tulip PSR_data name psr_data_tf1p200_sw34_minEps0.mat, with 'ELFO'/'Tulip'
% inserted). Set target='Tulip' to tag a tulip run the same way.
%
% PREREQUISITE: an ELFO min-ENERGY seed at this t_f must exist -- the base seed
% energy_elfo_freetf.mat (tf=7.5488 ND = 1.20x, from gen_elfo_energy_gravhom), or
% a tf-grid seed energy_elfo_tf####.mat (from gen_elfo_energy_tfsweep). Remember
% the energy band is WIDER than the eps=0-convergent band: some t_f reach fuel,
% some stall -- that map is the point of a t_f sweep (minfuel-tf-grid-strategy).
%
% REFERENCES:
%   [1] PSR/run_psr.m (the tulip entry driver this mirrors)
%   [2] gen_elfo_minfuel.m (the homotopy core); casadi_energy_freetf.m (solver)
%   [3] PSR/ELFO_RETARGET.md (the GTO->ELFO build record)

%% ------------------------------------------------------------------------
%% 0. Paths
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();  addpath(fullfile(here,'..','PSR'));
resDir = fullfile(here,'results');
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------
target   = 'ELFO';        % name tag written into the output filename
                          %   'ELFO' (this pipeline) | 'Tulip' (to tag a tulip run)
factor   = 1.20;          % t_f / tfMin_tulip (tfMin=6.2906939607 ND). Selects the
                          % ELFO energy seed at t_f = factor*tfMin. 1.20 = base seed.
seedSpec = 'auto';        % 'auto' finds the ELFO energy seed for this factor, or
                          % give an explicit /path/to/energy_seed.mat
epsMin   = 0;             % homotopy endpoint: 0 = bang-bang FUEL (default);
                          % >0 = smooth eps-optimal (e.g. 0.2) -- a regularized ramp
doExport = true;          % stage 4: write data products to ../PSR_data (target-tagged)
doVerify = true;          % stage 5: independent endpoint + defect check
movieMode= 'movie';       % stage 6: 'movie' (MP4+GIF) | 'preview' (3 stills) | 'none'
rerun    = false;         % true -> ignore any checkpoint and re-solve from the seed

% solver knobs (campaign defaults; edit only for experiments)
step0    = 0.20;          % initial epsilon step
stepMin  = 0.005;         % give up below this step (declares the sharpening wall)
maxIter  = 2000;          % IPOPT cap (tight solves)
looseIter= 500;           % IPOPT cap (loose probe)

%% ------------------------------------------------------------------------
%% Seed selection + output name
%% ------------------------------------------------------------------------
tf = factor * cfg.tfMin;
if strcmpi(seedSpec, 'auto')
    seedFile = '';
    cand = fullfile(resDir, sprintf('energy_elfo_tf%04d.mat', round(1000*tf)));  % tf-grid seed
    base = fullfile(resDir, 'energy_elfo_freetf.mat');                            % base seed (1.20x)
    if isfile(cand)
        seedFile = cand;
    elseif isfile(base)
        B = load(base, 'X');  if abs(B.X(8,end) - tf) < 0.02, seedFile = base; end
    end
    assert(~isempty(seedFile), ['no ELFO energy seed for factor %.3f (t_f=%.4f ND). ' ...
        'Build one: gen_elfo_energy_tfsweep (grid) or gen_elfo_energy_gravhom (base at 1.20x).'], factor, tf);
else
    seedFile = seedSpec;  assert(isfile(seedFile), 'seedSpec not found: %s', seedFile);
end

fTag = strrep(sprintf('%.3f', factor), '.', 'p');
eTag = strrep(sprintf('%g', epsMin), '.', 'p');
fprintf('\n=== RUN_ELFO_MINFUEL: target=%s  factor=%.3f (t_f=%.4f ND = %.2f d)  epsMin=%g ===\n', ...
        target, factor, tf, tf*p.tStar/86400, epsMin);
fprintf('    energy seed: %s\n', seedFile);

%% ------------------------------------------------------------------------
%% 2. SOLVE  (energy->fuel homotopy on the free-t_f two-primary solver)
%% ------------------------------------------------------------------------
if rerun
    ck = fullfile(resDir, sprintf('minfuel_%s_tf%s_minEps%s_ckpt.mat', target, fTag, eTag));
    if isfile(ck), delete(ck); end
end
fprintf('\n[stage 2] energy->fuel homotopy (eps 1 -> %g)...\n', epsMin);
outFile = gen_elfo_minfuel(struct('seedFile',seedFile,'target',target,'epsMin',epsMin, ...
    'step0',step0,'stepMin',stepMin,'maxIter',maxIter,'looseIter',looseIter,'resume',~rerun));

%% ------------------------------------------------------------------------
%% 3. PSR REFINEMENT  (indirect-steered mesh refinement) -- DEFERRED for ELFO
%% ------------------------------------------------------------------------
% run_psr stage 3 (refine_loop + pmp_refine_indicator) recovers costates and
% localizes switches with the tulip SINGLE-primary FIXED-t_f machinery; it does
% not apply to the two-primary free-t_f ELFO model without porting. The ELFO
% switch times therefore stand at the solver's converged mesh (N=4001). Porting
% switch-localization refinement to the two-primary model is future work.

%% ------------------------------------------------------------------------
%% 4. DATA EXPORT  (data products -> ../PSR_data/, target-tagged)
%% ------------------------------------------------------------------------
dataDir = fullfile(here, '..', 'PSR_data');
if doExport
    fprintf('\n[stage 4] DATA EXPORT (-> %s)...\n', dataDir);
    dataFile = elfo_export_data(outFile, dataDir);   % costates from the two-primary KKT duals
end

%% ------------------------------------------------------------------------
%% 5. VERIFY  (independent endpoint + solver-free defect check)
%% ------------------------------------------------------------------------
% First-order PMP diagnostics (primer alignment, mass-costate transversality)
% are carried in the stage-4 data product. The full 16-dim per-arc PMP
% PROPAGATION certificate (verify_direct_pmp) is tulip single-primary machinery
% and, like refinement, is deferred for the two-primary ELFO model.
if doVerify
    fprintf('\n[stage 5] VERIFY (endpoints + solver-free defect)...\n');
    SEEDFILE = outFile;   %#ok<NASGU>  verify_elfo_seed reads SEEDFILE
    verify_elfo_seed;
end

%% ------------------------------------------------------------------------
%% 6. CONTROL MOVIE  (transfer + control law, synced; ELFO orbit backdrop)
%% ------------------------------------------------------------------------
if ~strcmpi(movieMode, 'none')
    fprintf('\n[stage 6] CONTROL MOVIE (%s)...\n', movieMode);
    [~, ~, elfoTrace] = gto_elfo_endpoints(p, struct('point','apolune'));  % one-period ELFO backdrop
    Smf = load(outFile,'out');  ss = Smf.out.U(4,:);  nsw = sum(abs(diff(ss>0.5)));
    if epsMin == 0
        titleStr = sprintf('min-fuel GTO\\rightarrow%s, t_f=%.2fx min-time (%d-switch bang-bang)', target, factor, nsw);
    else
        titleStr = sprintf('GTO\\rightarrow%s, t_f=%.2fx min-time (smooth \\epsilon=%.3g)', target, factor, epsMin);
    end
    movieStem = fullfile(resDir, sprintf('movie_%s_tf%s_minEps%s', target, fTag, eTag));
    psr_movie(outFile, movieStem, titleStr, movieMode, elfoTrace(:,1:3));
end

fprintf('\n=== RUN_ELFO_MINFUEL DONE. Solution: %s ===\n', outFile);
if doExport, fprintf('    data product: %s\n', dataFile); end
