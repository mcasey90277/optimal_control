% REFINE_1N  Plan-required refinement spot-check for the 1 N fuel solve:
% re-solve at 1.5x the certified N, warm-started by mesh-refining the
% certified solution, and confirm |delta m_f| < 1 kg.
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here); addpath(here);
resDir = fullfile(here, 'results');

S = load(fullfile(resDir, 'MEE_M2_1N.mat'));
res = S.res;
p   = kepler_lt_params(1, 1500, 2000);

sigmaSrc = res.sigma;
Xsrc = res.fuel.X;  Usrc = res.fuel.U;  dLsrc = res.fuel.dL;
Nsrc = numel(sigmaSrc) - 1;
Ndst = round(1.5 * Nsrc);
fprintf('Refinement spot-check: Nsrc=%d -> Ndst=%d (1.5x)\n', Nsrc, Ndst);

sigmaDst = linspace(0, 1, Ndst+1).';
W = interp_warmstart(Xsrc, Usrc, dLsrc, sigmaSrc, sigmaDst);

x0 = W.X(:,1);
opts = struct('par', p, 'mode', 'fixedtf', 'eps', 0, 'tfTarget', res.tf, ...
    'x0', x0, 'maxIter', 1500, 'warmTight', true, 'printLevel', 0);

t0 = tic;
out = casadi_lt_mee(sigmaDst, W.X, W.U, W.dL, opts);
wallSec = toc(t0);

dmf = abs(out.m_f_kg - res.report.m_f_kg);
fprintf(['REFINE 1.5xN: success=%d defect=%.3e termErr=%.3e tfErr=%.3e ' ...
         'mf=%.4f kg (base %.4f kg, |dmf|=%.4f kg) sw=%d revs=%.4f wall=%.1fs\n'], ...
        out.success, out.maxDefect, out.termErr, out.tfErr, out.m_f_kg, ...
        res.report.m_f_kg, dmf, out.switches, out.dL/(2*pi), wallSec);
fprintf('GATE |dmf|<1 kg: %s\n', mat2str(dmf < 1));

save(fullfile(resDir, 'MEE_M2_1N_refine1p5x.mat'), 'out', 'dmf', 'Ndst');
