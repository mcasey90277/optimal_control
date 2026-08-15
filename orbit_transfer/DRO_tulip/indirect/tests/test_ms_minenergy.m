function test_ms_minenergy()
% TESTS_MS_MINENERGY  The fixed-tf min-energy multiple-shooting binding.
%
% A SYNTHETIC PROBLEM WITH A KNOWN ANSWER: pick lambda_0 at a benign CR3BP
% state, tune lambda_m(0) by a scalar root-find so that lambda_m(tf) = 0,
% propagate the min-energy PMP for tf, and take the endpoint (r, v) as rvf.
% Then (rv0, rvf, tf) is a min-energy BVP whose exact solution is that
% lambda_0. Perturb the exact trajectory (as golden_cells does) and demand:
%
%   1. RECOVERY: ms_minenergy converges (||R|| < 1e-10) to lambda_0 within
%      1e-8, in a few iterations, from a 1e-3-relative perturbed seed.
%   2. HAMILTONIAN: info.H is constant across junctions (fixed tf, autonomous
%      -> H is a first integral, NOT zero).
%   3. ACCEPTANCE: ss_bvp_accept on the SAME closures returns it unchanged.
%   4. INTERIOR: the throttle along the exact arc is interior somewhere (so
%      the test exercises the energy law, not the saturated min-time one).
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL; errors on failure)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                                  % indirect/
addpath(fullfile(fileparts(fileparts(fileparts(here))), 'costate_common'));

muStar = 0.012150585609624;  Tmax = 0.1756418;  c = 8.673746;
rv0 = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02];
lr  = 0.1*[15.6; 32.9; -0.09];  lv = 40*[-0.10; 0.045; -0.00015];   % s* ~ 0.4-0.7 on the arc
tf  = 0.5;
% tune lambda_m(0) so lambda_m(tf) = 0
lmEnd = @(lm0) endval(cr3bp_minenergy_prop(tf, [rv0; 1; lr; lv; lm0], false, Tmax, c, muStar), 14);
lm0 = fzero(lmEnd, [0, 5]);
y0  = [rv0; 1; lr; lv; lm0];
[yh, ~, T, Y] = cr3bp_minenergy_prop(tf, y0, false, Tmax, c, muStar);
rvf = yh(1:6);
sArc = zeros(numel(T),1);
for k = 1:numel(T)
    [~, ~, ax] = cr3bp_minenergy_pmp(Y(k,:)', Tmax, c, muStar);  sArc(k) = ax.s;
end
check('interior throttle on arc', min(sArc) > 0.02 && max(sArc) < 0.98, ...
      sprintf('s in [%.3f, %.3f]', min(sArc), max(sArc)));

% perturbed seed on K = 6 junctions
K = 6;  tG = linspace(0, tf, K+1);
Yg = interp1(T, Y, tG, 'pchip')';
rng(7);
Yp = Yg .* (1 + 1e-3*randn(size(Yg)));
seed = struct('tf', tf, 'tGrid', tG, 'Y', Yp);
[z, info] = ms_minenergy(rv0, rvf, tf, seed, Tmax, c, muStar, struct('tolR', 1e-10));
check('recovery: converged', info.converged, sprintf('normR %.1e', info.normR));
check('recovery: lambda0 within 1e-8', max(abs(z - y0(8:14))) < 1e-8, ...
      sprintf('max |dz| %.1e', max(abs(z - y0(8:14)))));
check('recovery: few iterations', info.iters <= 8, sprintf('%d', info.iters));
check('recovery: z is 7x1', isequal(size(z), [7 1]));
check('H constant across junctions', max(abs(info.H - info.H(1))) < 1e-8, ...
      sprintf('drift %.1e (H = %.4f)', max(abs(info.H - info.H(1))), info.H(1)));

% acceptance on the same closures
[zA, infoA] = ms_minenergy(rv0, rvf, tf, seed, Tmax, c, muStar, ...
                           struct('tolR', 1e-10, 'accept', true));
check('accept: ss gate accepted', infoA.accept.accepted, ...
      sprintf('dz %.1e normR0 %.1e', infoA.accept.dz, infoA.accept.normR0));
check('accept: same z', max(abs(zA - z)) == 0);
fprintf('test_ms_minenergy: ALL PASS\n');
end

function v = endval(y, k)
v = y(k);
end

function check(label, cond, detail)
if nargin < 3, detail = ''; end
if cond
    fprintf('  [PASS] %s %s\n', label, detail);
else
    fprintf('  [FAIL] %s %s\n', label, detail);
    error('test_ms_minenergy:fail', '%s', label);
end
end
