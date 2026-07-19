function tol = sosc_defaults()
% SOSC_DEFAULTS  Single source of SOSC certificate tolerances (canonical units,
% magnitudes O(1)). See process/DESIGN_sosc.md sec 6.
% OUTPUTS: tol - struct of scalar thresholds.
tol = struct( ...
    'recon',       1e-6, ...  % rebuild reproduces saved primal
    'drift',       1e-6, ...  % warm-resolve drift (report/warn, not fail)
    'stat',        1e-6, ...  % stationarity ||grad L||_inf
    'feas',        1e-8, ...  % equality residual / inequality violation
    'dual',        1e-8, ...  % inequality dual-sign violation
    'comp',        1e-6, ...  % complementarity max|lam*slack|
    'active',      1e-7, ...  % inequality slack -> active
    'mu',          1e-6, ...  % relative multiplier -> strongly-active
    'inertiaZero', 1e-6);     % relative pivot magnitude -> zero eigenvalue
% inertiaZero CALIBRATED on the 10 N row (DESIGN sec 6: "the 10 N integration
% test is where they are confirmed/tightened"). The bang-bang min-fuel KKT is
% near-singular with a ~270-dim genuine null space (reduced Hessian eigenvalues
% clustered at ~1e-10..1e-6; verified independently by dense eig of Z'HZ and of
% the full KKT, both reproducing DESIGN sec 11.1: KKT inertia (1865,1749,271),
% reduced (116,0,270)). Its ldl D-pivots show a WIDE, STABLE plateau: negative
% pivots with |.|>1e-5 number exactly 1749 (= true nneg) and stay 1749 through
% |.|>1e-3, while all spurious near-zero pivots sit below ~1e-5 (269th smallest
% |pivot|=1.23e-5, 270th=2.03e-4). Any zt in [1e-5,1e-3] recovers the true
% inertia (red.nneg=0, nFlat~265-269). With scale=normest(K)~125, inertiaZero=
% 1e-6 -> zt~1.25e-4, mid-plateau. The old 1e-9 (zt~1.25e-7) sat BELOW the
% spurious band and mis-signed ~56 noise pivots as negative curvature -> a
% spurious FAIL. 1e-6 errs conservatively (masks only truly-negligible
% curvature; genuine negatives here are O(1e-3+), well above threshold).
end
