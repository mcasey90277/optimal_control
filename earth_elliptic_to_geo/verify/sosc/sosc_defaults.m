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
    'inertiaZero', 1e-9);     % relative pivot magnitude -> zero eigenvalue
end
