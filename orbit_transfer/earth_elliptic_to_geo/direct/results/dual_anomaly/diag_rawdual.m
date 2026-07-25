function D = diag_rawdual(matPath)
% DIAG_RAWDUAL  Characterize the opti.dual vs raw lam_g disagreement on the
% defect rows, and test whether the raw multipliers fix the primer anomaly.
%
% CONTEXT. diag_t1_beta.m established, on the certified 10 N row:
%   - the RAW multiplier vector opti.lam_g satisfies KKT stationarity to
%     1.5e-14, and tangential dL/dbeta on burn nodes is 8e-17 (machine zero),
%     so beta-stationarity holds EXACTLY and the analytic primer derivation in
%     mee_primer_switch.m is not the problem;
%   - the Ldot-guard multipliers (the only control-dependent constraint that
%     could rotate the primer) are all zero -- eliminated;
%   - but opti.dual(conDef{k}) and the corresponding rows of opti.lam_g differ
%     by 76% relative.
% That isolates the anomaly to the DUAL EXTRACTION -- DESIGN_dual_map.md's
% H2(new) -- and this script characterizes the difference and tests the fix.
%
% WHAT IT MEASURES.
%   (1) Structure of the disagreement: per-state-row and per-node ratio of
%       opti.dual to raw lam_g, so a uniform scale, a per-row scale, a
%       permutation and an outright mismatch are distinguishable.
%   (2) The fix test: rebuild the nodal costate from the RAW lam_g defect rows
%       and re-run the primer/switching reconstruction. If the primer
%       misalignment collapses, the root cause is confirmed and the repair is
%       to source lamDef from opti.lam_g rather than opti.dual.
%
% INPUTS:
%   matPath - certified fuel row [char]; default results/MEE_M2_10N.mat
% OUTPUTS:
%   D - struct: .ratioMedByRow[7x1] .ratioMedAll .ratioSpread
%       .primerDeg_optidual .primerDeg_rawlamg .signPct_optidual
%       .signPct_rawlamg .tanAN_raw_med
%
% REFERENCES:
%   [1] diag_t1_beta.m (the T1 result this follows up).
%   [2] process/DESIGN_dual_map.md sec 3 T1 (H2-new branch: chase the
%       convention against low-level nlpsol lam_g).
import casadi.*
if nargin < 1 || isempty(matPath)
    matPath = fullfile(module_root(), 'results', 'MEE_M2_10N.mat');
end
saved = sosc_load_row(matPath);
par   = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
opts  = struct('par', par, 'mode', 'fixedtf', 'eps', 0, ...
    'tfTarget', saved.tfTarget, 'x0', saved.X(:,1), 'xf', saved.xf, ...
    'maxIter', saved.maxIter, 'warmTight', true, 'printLevel', 0, ...
    'returnModel', true);
fprintf('=== diag_rawdual: %s ===\n', saved.tag);
o = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts);
fprintf('warm re-solve: %s | maxDefect %.3e\n', o.ipoptStatus, o.maxDefect);

opti = o.model.opti;  sol = opti.debug;
lamv = full(sol.value(opti.lam_g));
sg   = saved.sigma(:);  Nseg = numel(sg) - 1;

di        = find(strcmp({o.model.creg.label}, 'defect'), 1);
lamDefRaw = reshape(lamv(o.model.creg(di).rows), 7, Nseg);
lamDefOpt = o.lamDef;

% --- (1) structure of the disagreement --------------------------------------
rat = lamDefOpt ./ lamDefRaw;
rat(~isfinite(rat)) = nan;
D.ratioMedByRow = median(rat, 2, 'omitnan');
D.ratioMedAll   = median(rat(:), 'omitnan');
D.ratioSpread   = iqr(rat(:));
fprintf('opti.dual / raw lam_g ratio, median by state row:\n');
lbl = {'P','ex','ey','hx','hy','m','t'};
for r = 1:7
    fprintf('   %-3s %+.6f\n', lbl{r}, D.ratioMedByRow(r));
end
fprintf('   all rows: median %+.6f | IQR %.3e\n', D.ratioMedAll, D.ratioSpread);
fprintf('   max|opt| %.3e   max|raw| %.3e\n', ...
    max(abs(lamDefOpt(:))), max(abs(lamDefRaw(:))));

% is it a per-node scale? report the node-wise ratio spread
ratNode = median(rat, 1, 'omitnan');
fprintf('   per-node ratio: median %+.6f  min %+.6f  max %+.6f\n', ...
    median(ratNode,'omitnan'), min(ratNode), max(ratNode));

% --- (2) the fix test --------------------------------------------------------
verOpt = verify_pmp_mee(o, par, sg, struct('eps', 0));
D.primerDeg_optidual = verOpt.primerMedianDeg;
D.signPct_optidual   = verOpt.overallSignPct;

oRaw = o;  oRaw.lamDef = lamDefRaw;
verRaw = verify_pmp_mee(oRaw, par, sg, struct('eps', 0));
D.primerDeg_rawlamg = verRaw.primerMedianDeg;
D.signPct_rawlamg   = verRaw.overallSignPct;

fprintf('\n--- FIX TEST ---------------------------------------------------\n');
fprintf('primer median: opti.dual %8.3f deg   raw lam_g %8.3f deg\n', ...
    D.primerDeg_optidual, D.primerDeg_rawlamg);
fprintf('sign agree   : opti.dual %8.2f %%     raw lam_g %8.2f %%\n', ...
    D.signPct_optidual, D.signPct_rawlamg);
fprintf('switchAlign  : opti.dual %8.3e      raw lam_g %8.3e\n', ...
    verOpt.maxSwitchAlignErr, verRaw.maxSwitchAlignErr);
fprintf('transversality lamM rel: opti.dual %.2e   raw lam_g %.2e\n', ...
    verOpt.lamMendRel, verRaw.lamMendRel);
fprintf('pass: opti.dual %d | raw lam_g %d\n', verOpt.pass, verRaw.pass);

D.tanAN_raw_med = verRaw.tangentialResidNormRelMedian;

outFile = fullfile(module_root(), 'results', 'dual_anomaly', ...
                   sprintf('diag_rawdual_%s.mat', saved.tag));
save(outFile, 'D');
fprintf('saved %s\n', outFile);
end
