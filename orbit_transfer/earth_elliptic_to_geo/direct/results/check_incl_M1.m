S = load('results/M1_3d_fixedLf.mat');
res = S.res;
Xf = res.fuel.X(1:6,end);
r = Xf(1:3); v = Xf(4:6);
h = cross(r,v);
incl = acosd(h(3)/norm(h));
fprintf('r_end = [%.6f %.6f %.6f]\n', r);
fprintf('v_end = [%.6f %.6f %.6f]\n', v);
fprintf('h = [%.6e %.6e %.6e]\n', h);
fprintf('incl_terminal_deg = %.6f\n', incl);
fprintf('X_rows = %d\n', size(res.fuel.X,1));
fprintf('cfg: thrustN=%g ctf=%g hx0=%g term=%s\n', res.cfg.thrustN, res.cfg.ctf, res.cfg.hx0, res.cfg.term);
fprintf('report: certified=%d defect=%.4e edge=%.4f revs=%.4f switches=%d apoBurnRatio=%.4f mf=%.4f dV=%.4f\n', ...
  res.report.certified, res.report.defect, res.report.edge, res.report.revs, res.report.switches, res.report.apoBurnRatio, res.report.m_f_kg, res.report.dV_kms);
