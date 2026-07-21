function outFile = gen_tulip_mintime(opts)
% GEN_TULIP_MINTIME  Solve the GTO->tulip minimum-time (hard all-burn, s==1)
% transfer directly, to certify tfMin_tulip and cross-check the pumpkyn/MS
% indirect reference (cfg.tfMin). The tulip analog of elfo/gen_elfo_mintime: it
% reuses the SAME target-agnostic solver casadi_mintime_freetf (in elfo/), with
% the SINGLE-primary clock (moonZone=0 -> kappa=r1^pSund) matching the tulip
% energy backbone. Loads the lowest tulip energy backbone, overrides throttle to
% s==1, minimizes t(tau_f), saves + reports the acceptance diagnostics.
%
% Warm start: sundman_minfuel/results/energy/energy_f1120.mat (tf 7.046 ND =
% 1.12x, single-primary, targets the far-Moon tulip point rvf the min-fuel front
% uses). X is 8-row [r;v;m;t] (solver appends cScale); U is 4-row [alpha;s]
% (solver drops the throttle row, enforcing s==1).
%
% INPUTS:
%   opts - (optional) struct: .seedFile[results/energy/energy_f1120.mat]
%          .maxIter[3000] .warmTight[false]
%
% OUTPUTS:
%   outFile - results/mintime_tulip_<insertionLabel>.mat (e.g. mintime_tulip_
%             tulipCampaign.mat): X[9xnN], U[4xnN] (alpha + s==1 row for
%             drop-in verify_elfo_seed compat), sigma, rv0, rvf, tauf0, tf(=tfMin),
%             mf, cScale, maxDefect, minR1, pSund, qSund, moonZone(=0), insertion
%             (= insMeta.label, the declared endpoint criterion; provenance only)
%
% REFERENCES:
%   [1] elfo/casadi_mintime_freetf.m (the solver, target-agnostic);
%   [2] elfo/gen_elfo_mintime.m (the ELFO sibling driver);
%   [3] Route B design spec docs/superpowers/specs/2026-07-15-elfo-mintime-route-b-design.md.

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here, '..', 'elfo'));      % casadi_mintime_freetf + verify_elfo_seed
resDir = fullfile(here,'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

% ---- INSERTION POINT (edit here to retarget) ---------------------------------
insertion = 'campaign';        % tulip: 'campaign'|'maxydot'|'apoapsis'  (elfo: 'nearest'|'apolune'|'perilune')
% insertion = 'maxydot';       % uncomment to use the max-ydot point (needs a matching energy seed)
% insertion = 'apoapsis';      % uncomment to use the slowest/apoapsis point (needs a matching seed)
[rv0, rvf, insMeta] = insertion_states('tulip', insertion);

% Default: the TWO-primary tulip energy seed (gen_tulip_energy_2p) -- re-meshed
% for the two-primary clock so the exact Hessian survives the near-Moon terminal.
% (The single-primary backbone energy_f####.mat crashes the exact-Hessian solve.)
% Filename carries the insertion label (gen_tulip_energy_2p writes it tagged the
% same way) -- both sides are in-scope drivers so they stay consistent.
seedFile = gd('seedFile', fullfile(resDir, sprintf('energy_tulip_2p_%s.mat', insMeta.label)));
S = load(seedFile);
% drift guard: the seed must be for the declared insertion point
assert(norm(S.rvf(:).' - rvf) < 1e-10 && norm(S.rv0(:).' - rv0) < 1e-10, ...
    'insertion:drift', ['seed endpoints differ from the declared %s insertion ' ...
    '(rvf %.2e, rv0 %.2e) -- regenerate the seed for this criterion'], ...
    insMeta.label, norm(S.rvf(:).'-rvf), norm(S.rv0(:).'-rv0));
fprintf('=== GEN_TULIP_MINTIME: seed %s (tf=%.4f ND, %.2f d) ===\n', ...
        seedFile, S.X(8,end), S.X(8,end)*p.tStar/86400);

% warm start: states from the seed, throttle overridden to s==1 (drop U row 4).
% TWO-primary clock (moonZone>0) tames the near-Moon terminal Hessian -- the
% single-primary clock leaves lunar gravity untamed near the tulip terminal
% (dMoon 28k km), which overflows IPOPT's exact Hessian -> MUMPS bus error.
X0 = S.X;  U0 = S.U(1:3,:);
o = struct('pSund',1.5,'qSund',4,'moonZone',gd('moonZone',0.15), ...
           'cBox',[0.10 8],'tfCapMult',4,'maxIter',gd('maxIter',3000), ...
           'warmTight',gd('warmTight',false));
out = casadi_mintime_freetf(S.sigma, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                            X0, U0, S.tauf0, o);

% acceptance diagnostics
rperi = norm(S.rv0(1:3) - [-p.muStar 0 0]);   % GTO perigee radius (ND)
rferr = norm(out.X(1:6,end) - S.rvf(:));
fprintf('  ipopt: %s   success=%d\n', out.ipoptStatus, out.success);
fprintf('  tfMin_tulip = %.6f ND (%.4f d)   mf=%.4f (prop %.1f%%)   cScale=%.4f\n', ...
        out.tf, out.tf*p.tStar/86400, out.mf, 100*(1-out.mf), out.cScale);
fprintf('  maxDefect=%.2e  maxUnit=%.2e  rendezvous=%.2e\n', out.maxDefect, out.maxUnit, rferr);
fprintf('  minR1=%.4f (GTO perigee=%.4f)  tMonotone=%d  primerAlign=%.3f deg\n', ...
        out.minR1, rperi, out.tMonotone, out.primerAlignDeg);

% cross-check vs the pumpkyn/MS indirect reference (cfg.tfMin)
dRef = out.tf - cfg.tfMin;
fprintf('\n  CROSS-CHECK vs indirect reference cfg.tfMin = %.6f ND:\n', cfg.tfMin);
fprintf('    direct tfMin_tulip = %.6f ND   delta = %+.3e ND (%+.4f%%)\n', ...
        out.tf, dRef, 100*dRef/cfg.tfMin);
if abs(dRef) < 1e-3
    fprintf('    -> MATCH: direct solve CERTIFIES the indirect reference.\n');
else
    fprintf('    -> DIFFERS: 6.2907 likely a different tulip target; this is the\n');
    fprintf('       true direct min-time to the front target (rvf). See notes.\n');
end

% save with a 4-row U (s==1 row) for drop-in compat with verify_elfo_seed / movie
X = out.X;  U = [out.U; ones(1,size(out.U,2))];  sigma = S.sigma; %#ok<NASGU>
tauf0 = S.tauf0; %#ok<NASGU>  % rv0/rvf: keep the DECLARED values from insertion_states (line 42) --
                              % guard-equal to the seed's (<=1e-10) but these are the ones the
                              % 'insertion' label above actually describes.
pSund = 1.5;  qSund = 4;  moonZone = o.moonZone; %#ok<NASGU>
tf = out.tf;  mf = out.mf;  cScale = out.cScale; %#ok<NASGU>
maxDefect = out.maxDefect;  minR1 = out.minR1; %#ok<NASGU>
insertion = insMeta.label; %#ok<NASGU>  provenance: the declared insertion criterion
fp = cr3bp_fingerprint(p, struct('tf', tf, 'insertion', insertion)); %#ok<NASGU>
outFile = fullfile(resDir, sprintf('mintime_tulip_%s.mat', insMeta.label));
save(outFile,'X','U','sigma','rv0','rvf','tauf0','tf','mf','cScale', ...
     'maxDefect','minR1','pSund','qSund','moonZone','insertion','fp');
fprintf('  saved %s\n', outFile);
fprintf('GEN_TULIP_MINTIME DONE\n');
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
