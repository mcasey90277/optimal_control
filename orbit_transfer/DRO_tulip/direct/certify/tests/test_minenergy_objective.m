function test_minenergy_objective()
% TEST_MINENERGY_OBJECTIVE  casadi_mintime_dro's fixed-tf min-ENERGY mode.
%
% The costate pipeline's Step 1 for a min-energy catalog is the same
% transcription with (a) t_f pinned and (b) the objective replaced by the
% Bertrand-Epenoy energy endpoint Int s^2 dt. Two things must hold:
%
%   1. ENERGY MODE SOLVES AND MEANS WHAT IT SAYS. On a flagship 12x12 cell,
%      warm-started from the min-time solution stretched to gamma*tf, with
%      opts.objective = 'energy' and opts.tfFix = gamma*tf: success, the
%      returned tf equals tfFix to round-off (all lifted copies), the
%      throttle is INTERIOR somewhere (thrMin well below 1 -- else it is
%      min-time in disguise), and out.J equals the Simpson quadrature of
%      s^2 recomputed from the returned nodes/midpoints.
%   2. DEFAULT PATH BYTE-IDENTICAL. Without the new options a min-time solve
%      (trapezoid AND Hermite-Simpson, N = 60, same seed) reproduces the
%      pre-change reference stored beside this test to the last bit -- the
%      option must not perturb the existing NLP.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL; errors on failure)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'lib'));
S  = load(fullfile(here, '..', '..', 'results', 'dsweep_12x12_cells.mat'));
cc = S.CELLS{6,6};
nN = 61;  N = nN - 1;
tS = linspace(0, cc.tf, nN);
Xs = interp1(cc.tNodes, cc.X', tS, 'pchip')';
Us = interp1(cc.tNodes, cc.U', tS, 'pchip')';
Us(1:3,:) = Us(1:3,:) ./ vecnorm(Us(1:3,:));
Us(4,:)   = min(max(Us(4,:), 0), 1);

% --- 1: energy mode ---------------------------------------------------------
gam = 1.2;  tfFix = gam*cc.tf;
oE = casadi_mintime_dro(cc.rv0, cc.rvf, cc.Tmax, cc.c, cc.muStar, N, Xs, Us, ...
        tfFix, struct('scheme','hermite-simpson', 'minAltKm',500, ...
                      'objective','energy', 'tfFix',tfFix, 'returnModel',true));
check('energy: success', oE.success, oE.ipoptStatus);
check('energy: tf == tfFix', abs(oE.tf - tfFix) < 1e-12, sprintf('%.3e', oE.tf - tfFix));
check('energy: lifted copies equal', oE.tfSpread < 1e-12);
check('energy: throttle interior', oE.thrMin < 0.9, sprintf('thrMin %.3f', oE.thrMin));
h  = tfFix/N;
J  = sum((h/6)*(oE.U(4,1:end-1).^2 + 4*oE.Um(4,:).^2 + oE.U(4,2:end).^2));
check('energy: out.J = Simpson(s^2)', abs(oE.J - J) < 1e-12*max(1,J), ...
      sprintf('J %.9f vs %.9f', oE.J, J));
check('energy: lamDef 7 x N', isequal(size(oE.lamDef), [7 N]));
check('energy: J < tf (s <= 1)', oE.J < tfFix);

% --- 2: default path byte-identical ---------------------------------------
ref = load(fullfile(here, 'mintime_ref_N60.mat'));
o1 = casadi_mintime_dro(cc.rv0, cc.rvf, cc.Tmax, cc.c, cc.muStar, N, Xs, Us, ...
        cc.tf, struct('scheme','trapezoid','minAltKm',500,'returnModel',true));
o2 = casadi_mintime_dro(cc.rv0, cc.rvf, cc.Tmax, cc.c, cc.muStar, N, Xs, Us, ...
        cc.tf, struct('scheme','hermite-simpson','minAltKm',500,'returnModel',true));
check('default trap: X, U, tf, lamDef bitwise', isequal(o1.X, ref.X1) && ...
      isequal(o1.U, ref.U1) && o1.tf == ref.tf1 && isequal(o1.lamDef, ref.lam1));
check('default HS: X, U, Um, tf, lamDef bitwise', isequal(o2.X, ref.X2) && ...
      isequal(o2.U, ref.U2) && isequal(o2.Um, ref.Um2) && o2.tf == ref.tf2 && ...
      isequal(o2.lamDef, ref.lam2));
fprintf('test_minenergy_objective: ALL PASS\n');
end

function check(label, cond, detail)
if nargin < 3, detail = ''; end
if cond
    fprintf('  [PASS] %s %s\n', label, detail);
else
    fprintf('  [FAIL] %s %s\n', label, detail);
    error('test_minenergy_objective:fail', '%s', label);
end
end
