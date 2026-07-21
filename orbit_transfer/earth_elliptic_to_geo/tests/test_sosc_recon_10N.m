% TEST_SOSC_RECON_10N  Task-1 gate: sosc_load_row normalizes the 10 N MEE_M2
% row, and rebuilding+re-solving the NLP from it reproduces the saved primal
% to tol.recon. If this fails, res.fp is insufficient and must be fixed
% before any SOSC work is trusted.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

matPath = fullfile(module_root(),'results','MEE_M2_10N.mat');
assert(isfile(matPath), 'need the certified 10 N cache MEE_M2_10N.mat');
saved = sosc_load_row(matPath);

assert(isequal(size(saved.X),[7, numel(saved.sigma)]), 'X shape');
assert(saved.thrustN==10 && saved.xf(1)==1, 'thrustN/xf');
assert(saved.tfTarget>0, 'tfTarget resolved');

par  = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
opts = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',saved.maxIter, ...
    'warmTight',true,'printLevel',0);
o = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts);
drift = max(abs(o.X(:) - saved.X(:)));
fprintf('recon drift ||x_rebuilt - x_saved||_inf = %.3e\n', drift);
assert(o.success, 'rebuild re-solve did not converge');
assert(drift < 1e-6, sprintf('recon drift %.3e >= 1e-6 tol.recon', drift));

% Normalizer correctness on the real xf-less MEE_M2 row (5 N) and a PSR-final
% row -- both lack fp.xf; expect the GEO default. No re-solve (cheap checks).
s5 = sosc_load_row(fullfile(module_root(),'results','MEE_M2_5N.mat'));
assert(isequal(s5.xf,[1;0;0;0;0]), '5 N: xf must default to GEO when fp.xf absent');
assert(abs(s5.tfTarget - 67.0194) < 1e-3 && s5.thrustN==5, '5 N: tfTarget/thrustN');
sP = sosc_load_row(fullfile(module_root(),'results','MEE_M2_1N_PSR_psr_final.mat'));
assert(strcmp(sP.kind,'PSR') && sP.thrustN==1, 'PSR: kind/thrustN');
assert(isequal(sP.xf,[1;0;0;0;0]), 'PSR: xf defaults to GEO');
assert(abs(sP.tfTarget - 335.7122) < 1e-3, 'PSR: tfTarget = fpFinal.tf (already resolved)');
assert(isequal(size(sP.X),[7,numel(sP.sigma)]), 'PSR: X/sigma shapes consistent');
fprintf('test_sosc_recon_10N PASSED\n');
