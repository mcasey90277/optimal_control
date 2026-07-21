% SMOKE_FIXEDTF  Confirm the pinned-t_f energy solve is well-posed (no drift):
% the leg-0 conversion (backbone -> free-t_f representation, mu=1, single-primary)
% should now converge cleanly instead of wandering off the warm start.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
E   = load(fullfile(cfg.dirs.energy, cfg.fname('energy', 1.20)));
tf0 = E.X(8,end);
o = struct('moonZone',0,'muGain',1,'tfTarget',tf0,'maxIter',150,'warmTight',false);
r = casadi_energy_freetf(E.sigma, E.rv0, E.rvf, p.Tmax, p.c, p.muStar, E.X, E.U, E.tauf0, o);
fprintf('\n[pinned-tf conversion] ok=%d status=%s\n', r.success, r.ipoptStatus);
fprintf('    maxDefect=%.2e  tf=%.4f (target %.4f)  cScale=%.4f  edge=%.1f%%\n', ...
        r.maxDefect, r.tf, tf0, r.cScale, 100*r.edge);
