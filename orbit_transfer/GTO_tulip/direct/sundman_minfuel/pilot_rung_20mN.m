function best = pilot_rung_20mN()
% PILOT_RUNG_20MN  Ladder-prep validation: one warm-chained 20 mN fuel rung.
%
% Chains from the certified nominal 1.15x solution (same t_f, new thrust) via
% chain_rung_seed_tulip, re-cleans the energy problem at 20 mN, then sharpens
% eps 1->0 through sundman_homotopy's hardened gates. PASS = certified=1 +
% clean boundSat + fp recorded. Artifacts under _T20mN tags; certified caches
% untouched. (2026-07-21 ladder-prep T6; spec sec 6.)
%
% OUTPUTS: best - sundman_homotopy best struct (.certified .epsReached ...)
% REFERENCES: [1] spec 2026-07-21-ladder-prep-design.md sec 6.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','..','..','cr3bp_common'));  setup_cr3bp_common();
cfg = minfuel_config(struct('thrustN', 0.020));
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
C   = load('sundman_minfuel_certified.mat');
tf  = C.out.X(8,end);
[sg, X0, U0, tauf0, fp] = chain_rung_seed_tulip(C, tf, p, struct('pilot','20mN')); %#ok<ASGLU>
tag = thrust_tag(cfg.thrustN);
saveFile = fullfile(here, 'results', sprintf('pilot_minfuel%s.mat', tag));
sched = [1 0.6 0.35 0.2 0.12 0.07 0.04 0.02 0.01 0.005 0.002 0.001 0];
[best, tbl] = sundman_homotopy(p, C.rv0, C.rvf, sg, X0, U0, tauf0, cfg.pSund, ...
                               sched, 3000, saveFile); %#ok<ASGLU>
sat = 'n/a'; if isfield(best,'boundSat'), sat = best.boundSat.worst; end
fprintf('\nPILOT 20mN TULIP: certified=%d epsReached=%.4g defect=%.2g sw=%d boundSatWorst=%s\n', ...
        best.certified, best.epsReached, best.maxDefect, best.switches, sat);
end
