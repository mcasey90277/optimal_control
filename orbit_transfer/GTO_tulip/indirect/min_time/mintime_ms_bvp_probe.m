function out = mintime_ms_bvp_probe()
% MINTIME_MS_BVP_PROBE  Route the GTO->tulip flagship min-time problem
%   through the SHARED costate-pipeline engine (costate_common/ms_bvp via
%   ms_tfmin) and gate the result like a catalog entry.
%
%   Stage A of the GTO->tulip costate-catalog program (STATUS_AND_ROADMAP
%   item 5): prove the generic multiple-shooting engine closes in the
%   25 mN / ~40-rev regime, cross-check against the module's own root
%   (mintime_ms_solve, ||R|| = 4e-9, 2026-07-13), and produce the first
%   tfMin-accepted + conjugate-tested entry in this regime.
%
%   Seed: fly the stored z8 with pumpkyn tfMinProp and cut into K
%   junctions (the conj_catalog_pass pattern). NOTE the harvest-from-
%   direct leg is NOT run: gen_tulip_mintime's artifact stores no defect
%   multipliers (finding, 2026-08-25) -- dual capture is a Stage-B item.
%
% INPUTS:
%   none (reads results/mintime_tulip_ms.mat beside this module)
%
% OUTPUTS:
%   out - struct: .z8 (accepted [lam0(7); tf]), .gates (per-gate pass),
%         .info (ms_bvp info incl. .conj), .pass (all gates)
%
% REFERENCES:
%   [1] costate_common/ms_bvp.m, ms_conjugate_test.m (the shared engine)
%   [2] DRO_tulip/indirect/ms_tfmin.m (the min-time wrapper)
%   [3] mintime_ms_solve.m (the module's own MS root, the cross-check)

here = fileparts(mfilename('fullpath'));
ot   = fileparts(fileparts(fileparts(here)));            % orbit_transfer
addpath(fullfile(ot, 'costate_common'));
addpath(fullfile(ot, 'DRO_tulip', 'indirect'));          % ms_tfmin (promotion pending)

[rv0, rvf, P] = mintime_params();
R = load(fullfile(here, 'results', 'mintime_tulip_ms.mat'));
z8ref = R.zt(:);                                          % [lam0(7); tf]
K = 60;                                                   % match the module root

fprintf('[probe] flagship min-time: tf_ref = %.6f ND, K = %d\n', z8ref(8), K);

%% Seed: fly the known root, cut into K junctions (shared builder):
seed = seed_from_z8(z8ref, rv0(1:6), K, P.Tmax25, P.c, P.muStar);

%% Solve through the shared engine, conjugate test on:
t0 = tic;
[z8, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, P.Tmax25, P.c, P.muStar, ...
                      struct('conjTest', true, 'wallSec', 1800));
wall = toc(t0);

%% Gates (catalog-entry standard):
g = struct();
g.converged = info.converged;
g.normR     = info.normR;                                 % <= 1e-10 target
g.dzRef     = norm(z8 - z8ref);                           % same extremal?
% flown arrival:
[~, yf] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(:); 1; z8(1:7)], ...
                                  P.Tmax25, P.c, P.muStar);
g.flownKm = norm(yf(end,1:3) - rvf(1:3)) * P.lStar;
% tfMin acceptance (accept unchanged; ss floor 1e-6):
zAcc = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z8, P.Tmax25, P.c, P.muStar);
g.accDz = norm(zAcc - z8);
g.conjPass  = info.conj.pass;
g.conjNCross = info.conj.nCrossings;
g.conjAtFinal = info.conj.atFinal;

pass = g.converged && g.normR <= 1e-9 && g.flownKm < 100 && g.accDz < 1e-6;

fprintf(['[probe] %.1f s | conv=%d ||R||=%.2e iters=%d | |z-zref|=%.2e | ', ...
         'flown %.3f km | accept |dz|=%.2e | conj pass=%d (nCross=%d atFinal=%d)\n'], ...
        wall, g.converged, g.normR, info.iters, g.dzRef, g.flownKm, ...
        g.accDz, g.conjPass, g.conjNCross, g.conjAtFinal);
fprintf('[probe] OVERALL: %s\n', ternary(pass, 'PASS', 'FAIL'));

out = struct('z8', z8, 'gates', g, 'info', info, 'pass', pass, 'wall', wall);
save(fullfile(here, 'results', 'mintime_ms_bvp_probe.mat'), '-struct', 'out');
end

% ------------------------------------------------------------------------
function s = ternary(c, a, b)
% TERNARY  a if c else b.  INPUTS: c;a;b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
