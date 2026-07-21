function outFile = gen_energy_seed(targetFactor, opts)
% GEN_ENERGY_SEED  Generate a min-ENERGY backbone (eps=1) for a given t_f factor.
%
% run_psr's default 'energy' seed reads a backbone file
% <energy>/energy_f####.mat at the requested factor. Those exist only for the
% factors someone pre-walked (1.12-1.15, 1.20-1.85); ask for a new t_f (e.g.
% factor=2.0) and minfuel_at_tf errors "no energy backbone for factor ...".
% This script makes one.
%
% Method (the energy_step / backbone_walk recipe, in a single callable): the
% convex eps=1 energy problem is the only primitive that continues crash-free
% across t_f, but only in SMALL steps -- so this WALKS from the nearest existing
% backbone to targetFactor in steps of opts.step, and at each rung does
%   (a) a LOOSE-warm-start continuation to the new t_f (a genuine move), then
%   (b) a TIGHT re-clean at that SAME t_f (keeps the KKT duals consistent so the
%       next step starts clean -- chaining loose steps without this blows inf_du
%       up to ~1e14 and diverges).
% Every rung is saved as its own energy_f####.mat, so a walk to 2.0 also fills
% in 1.90, 1.95, ... on the way.
%
% INPUTS:
%   targetFactor - t_f / t_f_min to generate a backbone for [scalar]
%   opts - (optional) struct:
%          step       - continuation step in factor [0.05]
%          seedFactor - starting backbone factor [auto: nearest existing]
%          maxIter    - IPOPT cap per solve [1500]
%          force      - regenerate rungs even if their file exists [false]
%
% OUTPUTS:
%   outFile - path to the generated energy_f####.mat at targetFactor (the file
%             run_psr's 'energy' seed then finds). Errors if the walk fails.
%
% EMPIRICAL USABLE BAND (2026-07-12): the walk converges only for factors in
% ~[1.12, 1.95]. Below 1.12x it DIVERGES on the first step (the transfer
% approaches min-time = all-burn, the throttle saturates, the convexity that
% makes energy easy degrades -- a near-min-time conditioning wall that hits the
% SMOOTH problem too, not just bang-bang). Above 1.95x it also diverges (2.00x:
% inf_du ~1e11, MUMPS OOM). Outside that band, expect failure without a finer
% perigee mesh / tighter warm start / different continuation. See
% ../../process/LOW_THRUST_MINFUEL_CAMPAIGN.md "ENERGY-BACKBONE floor/ceiling".
%
% RESUMABLE (important): every rung is saved as it succeeds and existing rungs
% are skipped, so if a step dies -- especially the sporadic UNCATCHABLE CasADi/
% IPOPT MEX FATAL crash (~1 in 10 solves) that kills the whole MATLAB process --
% just RE-RUN gen_energy_seed(targetFactor); it picks up from the last saved
% rung. (backbone_walk.sh isolates each step in its own process to survive that
% crash automatically; this single-process version relies on the re-run instead.)
%
% REFERENCES:
%   [1] ../sundman_minfuel/energy_step.m (the per-step recipe this mirrors)
%   [2] ../sundman_minfuel/orchestrate/backbone_walk.sh (the shell walker)
%   [3] ../../process/LOW_THRUST_MINFUEL_CAMPAIGN.md ("Down-sweep CRACKED": energy backbone)

if nargin < 2, opts = struct(); end
if ~isfield(opts,'step'),    opts.step = 0.05;     end
if ~isfield(opts,'maxIter'), opts.maxIter = 1500;  end
if ~isfield(opts,'force'),   opts.force = false;   end

here = fileparts(mfilename('fullpath'));
addpath(here);  setup_paths();
addpath(fullfile(here, '..', 'sundman_minfuel'));   % insertion_states (single-source; PSR vendors the rest)
cfg  = minfuel_config();
p    = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
edir = cfg.dirs.energy;
if ~exist(edir, 'dir'), mkdir(edir); end
efile = @(f) fullfile(edir, cfg.fname('energy', f));

% ---- INSERTION POINT (edit here to retarget) --------------------------------
insertion = 'campaign';        % tulip: 'campaign'|'maxydot'|'apoapsis'  (elfo: 'nearest'|'apolune'|'perilune')
% insertion = 'maxydot';       % uncomment to use the max-ydot point (needs a matching energy seed)
% insertion = 'apoapsis';      % uncomment to use the slowest/apoapsis point (needs a matching seed)
[rv0Decl, rvfDecl, insMeta] = insertion_states('tulip', insertion);

% already have it?
outFile = efile(targetFactor);
if isfile(outFile) && ~opts.force
    Lo = load(outFile, 'rv0', 'rvf');
    assert(norm(Lo.rvf(:).' - rvfDecl) < 1e-10 && norm(Lo.rv0(:).' - rv0Decl) < 1e-10, ...
        'insertion:drift', ['cached target backbone %s is for a different insertion than the ' ...
        'declared %s -- delete it or change the criterion'], outFile, insMeta.label);
    fprintf('gen_energy_seed: backbone already exists: %s\n', outFile);
    return
end

% ---- find the nearest existing backbone to seed from -----------------------
if isfield(opts,'seedFactor') && ~isempty(opts.seedFactor)
    seedFactor = opts.seedFactor;
    assert(isfile(efile(seedFactor)), 'seedFactor %.3f has no backbone (%s)', ...
           seedFactor, efile(seedFactor));
else
    d = dir(fullfile(edir, 'energy_f*.mat'));
    have = zeros(1, numel(d));
    for k = 1:numel(d), have(k) = cfg.fparse(d(k).name); end
    have = have(~isnan(have));
    assert(~isempty(have), 'gen_energy_seed:noBackbones', ...
        ['no existing energy backbone to seed from in %s -- bootstrap one first ' ...
         '(the campaign''s first energy solution)'], edir);
    [~, ix] = min(abs(have - targetFactor));
    seedFactor = have(ix);
end
fprintf('gen_energy_seed: target %.3f, nearest backbone %.3f, step %.3f\n', ...
        targetFactor, seedFactor, opts.step);

% drift guard: the base backbone this walk starts from must match the
% declared insertion point (rv0/rvf are propagated unchanged through every
% rung -- see the save() below -- so checking the base backbone is sufficient).
base = load(efile(seedFactor), 'rv0', 'rvf');
assert(norm(base.rvf(:).' - rvfDecl) < 1e-10 && norm(base.rv0(:).' - rv0Decl) < 1e-10, ...
    'insertion:drift', ['seed endpoints differ from the declared %s insertion ' ...
    '(rvf %.2e, rv0 %.2e) -- regenerate the seed for this criterion'], ...
    insMeta.label, norm(base.rvf(:).'-rvfDecl), norm(base.rv0(:).'-rv0Decl));

% ---- build the continuation ladder seedFactor -> targetFactor --------------
sgn   = sign(targetFactor - seedFactor);
rungs = [];
f = seedFactor;
while abs(f - targetFactor) > 1e-9
    f = round(f + sgn*min(opts.step, abs(targetFactor - f)), 6);
    rungs(end+1) = f; %#ok<AGROW>
end
if isempty(rungs)
    error('gen_energy_seed:sameFactor', 'target %.3f equals the seed factor', targetFactor);
end
fprintf('gen_energy_seed: ladder %s\n', mat2str(rungs));

% ---- walk it ---------------------------------------------------------------
prevFactor = seedFactor;
for f = rungs
    rungFile = efile(f);
    if isfile(rungFile) && ~opts.force
        fprintf('  [%.3f] exists -- skip\n', f);  prevFactor = f;  continue
    end
    E   = load(efile(prevFactor));
    tf  = f * cfg.tfMin;   tfp = prevFactor * cfg.tfMin;
    Xk  = E.X;  Xk(8,:) = Xk(8,:) * (tf / tfp);          % rescale time state
    % (a) loose continuation to the new t_f
    o = casadi_minfuel_sundman(E.sigma, tf, E.rv0, E.rvf, p.Tmax, p.c, p.muStar, ...
            Xk, E.U, E.tauf0, cfg.pSund, opts.maxIter, 1, false);
    fprintf('  [%.3f<-%.3f] loose ok=%d defect=%.2g', f, prevFactor, o.success, o.maxDefect);
    if ~(o.success && o.maxDefect < 1e-6)
        fprintf('  -- FAILED, walk stops\n');
        error('gen_energy_seed:stepFailed', ...
            ['loose continuation %.3f -> %.3f failed (ok=%d defect=%.2g); try a smaller ' ...
             'opts.step'], prevFactor, f, o.success, o.maxDefect);
    end
    % (b) tight re-clean at the same t_f (keeps duals consistent for the next step)
    oT = casadi_minfuel_sundman(E.sigma, tf, E.rv0, E.rvf, p.Tmax, p.c, p.muStar, ...
             o.X, o.U, E.tauf0, cfg.pSund, opts.maxIter, 1, true);
    fprintf('  reclean ok=%d defect=%.2g edge=%.1f%%\n', oT.success, oT.maxDefect, 100*oT.edge);
    if ~(oT.success && oT.maxDefect < 1e-6)
        error('gen_energy_seed:recleanFailed', ...
            'tight re-clean at %.3f failed (ok=%d defect=%.2g)', f, oT.success, oT.maxDefect);
    end
    X = oT.X;  U = oT.U;  factor = f;  sigma = E.sigma;  tauf0 = E.tauf0; %#ok<NASGU>
    rv0 = E.rv0;  rvf = E.rvf;  insertion = insMeta.label; %#ok<NASGU>
    % NOTE: filename (energy_f####.mat) is NOT tagged with the insertion label --
    % it is the shared cfg.fname('energy',...) convention read by minfuel_at_tf.m,
    % elfo/smoke_fixedtf.m, elfo/smoke_energy_freetf.m, gen_elfo_energy_gravhom.m,
    % and others outside this task's touched-file set. Retagging would require
    % updating all of those readers too (see task-4-report.md concerns); the
    % 'insertion' field still records the criterion for provenance.
    save(rungFile, 'X', 'U', 'factor', 'tauf0', 'sigma', 'rv0', 'rvf', 'insertion');
    fprintf('  saved %s\n', rungFile);
    prevFactor = f;
end

outFile = efile(targetFactor);
fprintf('gen_energy_seed DONE: %s\n', outFile);
end
