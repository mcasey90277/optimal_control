function test_cr3bp_minenergy_pmp()
% TEST_CR3BP_MINENERGY_PMP  The min-energy PMP vector field, Jacobian, and
% STM propagator, checked against pumpkyn's min-time field and against
% finite differences.
%
% FOUR TESTS
%   1. SATURATION = MIN-TIME. With a large |lambda_v| the energy throttle
%      s* = clip((T/2)(|lambda_v|/m + lambda_m/c), 0, 1) sits at 1, and the
%      14-state field AND its Jacobian must equal pumpkyn.cr3bp.tfMinEoM's
%      (which has u = 1 there) to round-off. This pins every sign, frame
%      and costate convention to the reference the catalogs are built on.
%   2. INTERIOR THROTTLE. Scale lambda_v down so s* is interior; s must equal
%      the formula, and the Jacobian must match a central-difference
%      Jacobian of F (this is where the energy field differs from
%      min-time: s depends on lambda_v, m, lambda_m).
%   3. HAMILTONIAN CONSERVED. Propagate 0.4 ND from an interior-throttle
%      state; H = s^2 + lambda'f is autonomous and must be constant.
%   4. STM = d(flow)/dy0. Central differences of the propagated endpoint
%      must match the returned PHI.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per test; errors on failure)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

muStar = 0.012150585609624;
Tmax   = 0.1756418;
c      = 8.673746;
I14    = reshape(eye(14), [], 1);

% a benign state away from the Moon; costates of catalog-like magnitude
rv  = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02];
m   = 0.97;
lr  = [15.6; 32.9; -0.09];
lv  = [-0.10; 0.045; -0.00015];
lm  = 0.13;

% --- 1: saturation --------------------------------------------------------
ySat = [rv; m; lr; 400*lv; lm];            % (T/2)|lambda_v|/m ~ 4 > 1 -> s = 1
[F, A, aux] = cr3bp_minenergy_pmp(ySat, Tmax, c, muStar);
Fref = pumpkyn.cr3bp.tfMinEoM(0, [ySat; I14], Tmax, c, muStar);
Aref = reshape(Fref(15:210), 14, 14);
check('sat: s = 1', abs(aux.s - 1) < 1e-15, sprintf('s = %.6g', aux.s));
check('sat: F == tfMinEoM F', max(abs(F - Fref(1:14))) < 1e-12, ...
      sprintf('max diff %.2e', max(abs(F - Fref(1:14)))));
relSat = max(abs(A(:) - Aref(:))) / max(abs(Aref(:)));   % entries reach 5e3
check('sat: A == tfMinEoM A', relSat < 1e-11, sprintf('rel %.2e', relSat));

% --- 2: interior throttle -------------------------------------------------
yInt = [rv; m; lr; 40*lv; lm];             % s* ~ 0.4
[F, A, aux] = cr3bp_minenergy_pmp(yInt, Tmax, c, muStar);
sExp = (Tmax/2)*(40*sqrt(sum(lv.^2))/m + lm/c);
check('int: s interior', aux.s > 0.02 && aux.s < 0.98, sprintf('s = %.4f', aux.s));
check('int: s formula', abs(aux.s - sExp) < 1e-14);
Afd = zeros(14);  h = 1e-6;
for k = 1:14
    e = zeros(14,1);  e(k) = h;
    Afd(:,k) = (cr3bp_minenergy_pmp(yInt+e, Tmax, c, muStar) ...
              - cr3bp_minenergy_pmp(yInt-e, Tmax, c, muStar)) / (2*h);
end
relA = max(abs(A(:) - Afd(:))) / max(abs(Afd(:)));
check('int: A == FD Jacobian', relA < 1e-6, sprintf('rel %.2e', relA));
% mass row: mdot = -s T/c
check('int: mdot = -sT/c', abs(F(7) + aux.s*Tmax/c) < 1e-15);

% --- 3: H conserved along the flow ---------------------------------------
[yh, ~, T, Y] = cr3bp_minenergy_prop(0.4, yInt, false, Tmax, c, muStar);
Hs = zeros(numel(T),1);
for k = 1:numel(T)
    [~, ~, ax] = cr3bp_minenergy_pmp(Y(k,:)', Tmax, c, muStar);
    Hs(k) = ax.H;
end
check('prop: H conserved', max(abs(Hs - Hs(1))) < 1e-9, ...
      sprintf('drift %.2e (H = %.6f)', max(abs(Hs - Hs(1))), Hs(1)));
check('prop: endpoint = last row', max(abs(yh - Y(end,:)')) == 0);

% --- 4: STM by finite differences ----------------------------------------
[~, PHI] = cr3bp_minenergy_prop(0.4, yInt, true, Tmax, c, muStar);
PHIfd = zeros(14);  h = 1e-6;
for k = 1:14
    e = zeros(14,1);  e(k) = h;
    yp = cr3bp_minenergy_prop(0.4, yInt+e, false, Tmax, c, muStar);
    ym = cr3bp_minenergy_prop(0.4, yInt-e, false, Tmax, c, muStar);
    PHIfd(:,k) = (yp - ym)/(2*h);
end
relP = max(abs(PHI(:) - PHIfd(:))) / max(abs(PHIfd(:)));
check('prop: STM == FD flow Jacobian', relP < 1e-6, sprintf('rel %.2e', relP));
fprintf('test_cr3bp_minenergy_pmp: ALL PASS\n');
end

function check(label, cond, detail)
if nargin < 3, detail = ''; end
if cond
    fprintf('  [PASS] %s %s\n', label, detail);
else
    fprintf('  [FAIL] %s %s\n', label, detail);
    error('test_cr3bp_minenergy_pmp:fail', '%s', label);
end
end
