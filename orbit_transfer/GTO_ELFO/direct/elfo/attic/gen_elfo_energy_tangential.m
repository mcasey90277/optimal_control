function outFile = gen_elfo_energy_tangential(opts)
% GEN_ELFO_ENERGY_TANGENTIAL  From-SCRATCH min-ENERGY (eps=1) ELFO backbone with
% NO tulip seed: propagate max-thrust velocity-aligned ("tangential") steering
% from the GTO, map that arc into Sundman coordinates targeting the ELFO
% apolune, and attempt the energy solve directly.  [INDEPENDENT CHECK route.]
%
% This is the target-agnostic analog of how the tulip backbone was first
% bootstrapped (build_guess.m), but WITHOUT the tulip-specific converged
% indirect costate: the tangential guess reaches the GTO neighborhood only, so
% the NLP must close a large multi-rev rendezvous gap. build_guess's own note
% warns this "may converge slowly, to a different local minimum, or not at all."
% If it converges it gives a tulip-free ELFO seed; if not, the homotopy
% (gen_elfo_energy_backbone) carries the campaign.
%
% INPUTS:
%   opts - (optional) struct: .factor [1.50] .maxIter [1500] .elfo (gto_elfo opts)
%
% OUTPUTS:
%   outFile - results/energy_elfo_tangential.mat if it converges (empty if not)
%
% REFERENCES:
%   [1] attic/build_guess.m ('tangential' mode), attic/lt_dynamics.m
%   [2] sundman_seed_map.m (physical arc -> Sundman mesh, endpoints pinned)

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
factor  = gd('factor', 1.50);
maxIter = gd('maxIter', 1500);
elfoOpts = gd('elfo', struct());

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
% sundman_seed_map + lt_dynamics are not vendored into PSR/lib. Append with
% '-end' so PSR/lib's 14-arg casadi_minfuel_sundman keeps priority over the
% stale 13-arg copies in sundman_minfuel/attic.
addpath(fullfile(here,'..','..','..','..','GTO_tulip','direct','sundman_minfuel'),'-end');
addpath(fullfile(here,'..','attic'),'-end');
resDir = fullfile(here,'results');  if ~exist(resDir,'dir'), mkdir(resDir); end
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
tf  = factor * cfg.tfMin;

[rv0, rvf_elfo] = gto_elfo_endpoints(p, elfoOpts);
fprintf('=== GEN_ELFO_ENERGY_TANGENTIAL (from scratch): factor=%.2f tf=%.4f ND ===\n', factor, tf);

% --- propagate max-thrust tangential steering from the GTO ------------------
odeo  = odeset('RelTol',1e-9,'AbsTol',1e-11);
steer = @(t,x) lt_dynamics(x, x(4:6)./sqrt(sum(x(4:6).^2)), p.Tmax, p.c, p.muStar);
[tau, y] = ode113(steer, [0 tf], [rv0(:); 1], odeo);
[tau, keep] = unique(tau, 'stable');  y = y(keep, :);
fprintf('  tangential arc: %d nodes, final dist Moon = %.0f km, mass frac = %.3f\n', ...
        numel(tau), norm(y(end,1:3)-[1-p.muStar 0 0])*p.lStar, y(end,7));

% --- map to Sundman mesh, target the ELFO apolune (endpoints pinned) --------
Xseed = y(:,1:7).';                          % [7xM] [r;v;m]
Vg    = Xseed(4:6,:);
alpha = Vg ./ max(sqrt(sum(Vg.^2,1)), 1e-9); % unit tangential direction
Useed = [alpha; ones(1,size(Xseed,2))];      % full throttle guess
[sigma, X0, U0, tauf0] = sundman_seed_map(Xseed, Useed, tf, tau, cfg.pSund, p.muStar, rv0, rvf_elfo);
fprintf('  seed mapped: %d Sundman nodes, tauf0=%.4f\n', numel(sigma), tauf0);

% --- attempt the energy solve (loose, then tight) ---------------------------
oL = casadi_minfuel_sundman(sigma,tf,rv0,rvf_elfo,p.Tmax,p.c,p.muStar, ...
                            X0,U0,tauf0,cfg.pSund,maxIter,1,false);
fprintf('  energy LOOSE : ok=%d defect=%.2g edge=%.1f%% status=%s\n', ...
        oL.success, oL.maxDefect, 100*oL.edge, oL.ipoptStatus);
best = oL;
if oL.success
    oT = casadi_minfuel_sundman(sigma,tf,rv0,rvf_elfo,p.Tmax,p.c,p.muStar, ...
                                oL.X,oL.U,tauf0,cfg.pSund,maxIter,1,true);
    fprintf('  energy TIGHT : ok=%d defect=%.2g edge=%.1f%% status=%s\n', ...
            oT.success, oT.maxDefect, 100*oT.edge, oT.ipoptStatus);
    if oT.success, best = oT; end
end

outFile = '';
if best.success && best.maxDefect < 1e-6
    X = best.X;  U = best.U;  rvf = rvf_elfo(:).';  pSund = cfg.pSund; %#ok<NASGU>
    outFile = fullfile(resDir, 'energy_elfo_tangential.mat');
    save(outFile, 'X','U','factor','tauf0','sigma','rv0','rvf','pSund');
    fprintf('GEN_ELFO_ENERGY_TANGENTIAL CONVERGED: %s\n', outFile);
else
    fprintf('GEN_ELFO_ENERGY_TANGENTIAL DID NOT CONVERGE (defect=%.2g) -- from-scratch route fails here; homotopy carries.\n', best.maxDefect);
end
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
