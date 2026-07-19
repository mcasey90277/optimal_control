function tol = sosc_defaults()
% SOSC_DEFAULTS  Single source of SOSC certificate tolerances (canonical units,
% magnitudes O(1)). See process/DESIGN_sosc.md sec 6, 11.4.
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
    'inertiaZero', 1e-9, ...  % relative eigenvalue magnitude -> zero eigenvalue
    'maxNullDim',  10000);    % n at/below which the dense null-space is formed
% inertiaZero REVERTED to 1e-9 (DESIGN sec 11.4, 2026-07-19): this value is
% correct for gold-standard dense `eig(full(K))`, which is now the primary
% inertia method (sosc_inertia). On the 10 N row, `eig` gives the true reduced
% inertia (116,0,270) robustly for relative zt in [1e-9,1e-8] -- a clean
% spectral gap (~270 eigenvalues at 1e-10, nothing until ~1e-4) makes the
% classification insensitive to zt across that whole window. The earlier
% Amendment-B raise to 1e-6 was a workaround for `ldl`-pivot-sign inertia
% (count_inertia), which is UNRELIABLE on this near-singular KKT: at zt=1e-9
% ldl mis-signs 56 spurious negative pivots, and ldl's own correct window
% (zt in [1e-7,1e-6]) is disjoint from eig's -- there is no single zt that
% works for both methods, so tuning zt to ldl's window is not robust across
% rungs. `eig` is therefore the gold standard and 1e-9 is restored.
% maxNullDim=10000 is the FINAL-method size guard (DESIGN sec 12.1): the direct
% reduced-Hessian method forms Z=null(full(A)) (dense), tractable only for
% n<=maxNullDim (covers 10/5/2.5 N). Above it (1 N, 0.5 N) the dense null-space
% is intractable, so sosc_inertia returns IN.robust=false, IN.method='scale-skip'
% and the verdict is INCONCLUSIVE-by-scale rather than an unvalidated inertia.
end
