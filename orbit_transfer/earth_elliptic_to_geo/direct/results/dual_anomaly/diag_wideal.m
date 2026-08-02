% DIAG_WIDEAL  Test: do the +-1.01 alpha box bounds contaminate the duals?
% Re-solve M1 at eps=0 warmTight with alpha bounds widened to +-2 (the unit
% equality already constrains alpha) and check the primer alignment.
SC = '/private/tmp/claude-501/-Users-msc-Desktop-optimal-control/7953f18e-8fd4-48b1-a47f-b565d00b2ad7/scratchpad';
addpath(SC);
cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo');
S = load('results/M1_3d_fixedLf.mat');  res = S.res;
p  = kepler_lt_params(res.cfg.thrustN, 1500, 2000);
term = geo_terminal('fixed', p, res.Lf);
o = casadi_lt_2body_wideal(res.sg, res.fuel.X, res.fuel.U, res.fuel.tauf0, term, ...
    struct('par',p,'mode','fixedtf','eps',0,'tfTarget',res.tf,'rv0',res.rv0, ...
           'maxIter',3000,'warmTight',true,'printLevel',3));
fprintf('WIDEAL: status=%s defect=%.2e mf=%.2f sw=%d primerAlignDeg=%.4f\n', ...
    o.ipoptStatus, o.maxDefect, o.m_f_kg, o.switches, o.primerAlignDeg);
ver = verify_pmp_2body(o, p);
save('results/M1_wideal_probe.mat', 'o', 'ver');
fprintf('VERIFY(wideal): primer=%.3f burn=%.1f coast=%.1f lamMrel=%.1e betaSpread=%.1f pass=%d\n', ...
    ver.primerMeanDeg, ver.burnSignPct, ver.coastSignPct, ver.lamMendRel, ver.betaSpreadPct, ver.pass);
