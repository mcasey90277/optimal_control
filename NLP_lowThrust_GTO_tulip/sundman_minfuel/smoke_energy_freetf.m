% SMOKE_ENERGY_FREETF  Machinery check for casadi_energy_freetf: does the
% free-t_f / two-primary-clock NLP construct and step from a real backbone?
% Not a convergence run -- small iteration caps, we only certify "builds + moves".
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
E   = load(fullfile(cfg.dirs.energy, cfg.fname('energy', 1.20)));
fprintf('backbone f1.20: N=%d  tf_ws=%.4f  tauf0=%.4f\n', numel(E.sigma)-1, E.X(8,end), E.tauf0);

% --- Test A: moonZone=0 (single-primary recovery) + free t_f ----------------
% Same clock as the backbone, so the warm start (cScale=1) should be nearly
% feasible; this validates the slack-state free-t_f reproduces the original.
oA = casadi_energy_freetf(E.sigma, E.rv0, E.rvf, p.Tmax, p.c, p.muStar, ...
        E.X, E.U, E.tauf0, struct('moonZone',0,'muGain',1,'maxIter',40));
fprintf('\n[A single-primary, free tf] ok=%d status=%s\n', oA.success, oA.ipoptStatus);
fprintf('    maxDefect=%.2e  tf=%.4f  cScale=%.4f  edge=%.1f%%\n', ...
        oA.maxDefect, oA.tf, oA.cScale, 100*oA.edge);

% --- Test B: moonZone=0.15 (two-primary clock) + free t_f -------------------
% New clock -> warm start is NOT feasible (mesh<->time relation changed), so
% expect nonzero defect; we only certify the two-primary NLP constructs & steps.
oB = casadi_energy_freetf(E.sigma, E.rv0, E.rvf, p.Tmax, p.c, p.muStar, ...
        E.X, E.U, E.tauf0, struct('moonZone',0.15,'muGain',1,'maxIter',40));
fprintf('\n[B two-primary, free tf] ok=%d status=%s\n', oB.success, oB.ipoptStatus);
fprintf('    maxDefect=%.2e  tf=%.4f  cScale=%.4f  edge=%.1f%%\n', ...
        oB.maxDefect, oB.tf, oB.cScale, 100*oB.edge);

% --- Test C: muGain=0 (well-less root) — the easy end of the gravity homotopy
oC = casadi_energy_freetf(E.sigma, E.rv0, E.rvf, p.Tmax, p.c, p.muStar, ...
        E.X, E.U, E.tauf0, struct('moonZone',0.15,'muGain',0,'maxIter',40));
fprintf('\n[C two-primary, muGain=0 well-less] ok=%d status=%s\n', oC.success, oC.ipoptStatus);
fprintf('    maxDefect=%.2e  tf=%.4f  cScale=%.4f  edge=%.1f%%\n', ...
        oC.maxDefect, oC.tf, oC.cScale, 100*oC.edge);

fprintf('\nSMOKE DONE (machinery constructs + steps on all three).\n');
