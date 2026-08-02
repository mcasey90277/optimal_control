% DIAG_DUALPOLISH  Re-solve M1 at eps=0 warmTight from its own solution and
% check whether freshly-polished duals restore the primer law (campaign
% protocol: the eps=0 re-solve that regenerates consistent duals).
cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo');
S = load('results/M1_3d_fixedLf.mat');  res = S.res;
p  = kepler_lt_params(res.cfg.thrustN, 1500, 2000);
term = geo_terminal('fixed', p, res.Lf);
o = casadi_lt_2body(res.sg, res.fuel.X, res.fuel.U, res.fuel.tauf0, term, ...
    struct('par',p,'mode','fixedtf','eps',0,'tfTarget',res.tf,'rv0',res.rv0, ...
           'maxIter',3000,'warmTight',true,'printLevel',3));
fprintf('POLISH: status=%s defect=%.2e mf=%.2f sw=%d edge=%.3f primerAlignDeg=%.4f\n', ...
    o.ipoptStatus, o.maxDefect, o.m_f_kg, o.switches, o.edge, o.primerAlignDeg);
ver = verify_pmp_2body(o, p);
save('results/M1_dualpolish.mat', 'o', 'ver');
fprintf('VERIFY(polished): primer=%.3f burn=%.1f coast=%.1f lamMrel=%.1e pass=%d\n', ...
    ver.primerMeanDeg, ver.burnSignPct, ver.coastSignPct, ver.lamMendRel, ver.pass);
