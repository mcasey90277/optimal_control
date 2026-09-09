function ok = test_ms_tfmin_hom()
%% Purpose:
%
%   Tests ms_tfmin_hom -- the HOMOGENEOUS min-time binding: objective
%   multiplier rho free, H = rho + lam'f = 0 at t_f, normalisation
%   rho^2 + |lam_0|^2 = 1. Built because the normal chart (rho = 1) runs to
%   infinity on the DRO -> tulip fast family as thrust falls toward 72 mN
%   (|lam_0| 46 -> 1449 along the arc, every arclength step pure costate):
%   in that chart there is no finite corner to round. On the sphere the
%   multipliers stay bounded and rho -> 0 is visible as the branch losing
%   normality.
%
%   Fixtures: the golden DRO cell (a certified 1 N min-time root, rho = 1).
%   (1) Re-normalised, it must be a root of the homogeneous system with
%       rho = 1/sqrt(1 + |lam_0|^2) and the SAME trajectory (t_f equal).
%   (2) The homogeneous Jacobian must agree with finite differences.
%   (3) The solve from a mildly perturbed seed must return to the same t_f.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
G = load(fullfile(here, 'golden_cells_data.mat'));
c = G.cells(1);                                   % dro
z8 = c.z8;  tf = z8(8);
rv0 = c.rv0(:);  rvf = c.rvf(:);  Tnd = c.Tnd;  cnd = c.cnd;  mu = c.muStar;
K = 12;
seed = seed_from_z8(z8, rv0, K, Tnd, cnd, mu);
% seed_from_z8 FLIES z8 to cut the junctions, so they carry integrator
% error (~1e-7 in the ms residual). The fixture for "a normal root maps to
% a homogeneous root" must be a CONVERGED normal root: polish first.
[z0, i0] = ms_tfmin(rv0, rvf, seed, Tnd, cnd, mu, struct('tolR', 1e-11));
assert(i0.converged, 'fixture: normal polish did not converge');
seed.Y = [i0.Y, i0.Y(:,end)];  seed.tf = z0(8);  seed.tGrid = linspace(0, z0(8), K+1);

% (1) the normal root, re-normalised onto the sphere, is a homogeneous root
a = 1/sqrt(1 + norm(seed.Y(8:14,1))^2);
seedH = seed;  seedH.Y(8:14, :) = a*seedH.Y(8:14, :);  seedH.extra = a;   % rho = a
[~, ih] = ms_tfmin_hom(rv0, rvf, seedH, Tnd, cnd, mu, struct('assembleOnly', true));
ok = chk(ok, norm(ih.R, inf) < 1e-9, ...
         sprintf('re-normalised normal root is a homogeneous root: |R| = %.1e', norm(ih.R, inf)));
ok = chk(ok, abs(ih.p(end) - a) < 1e-14 && abs(ih.p(end-1) - z0(8)) < 1e-14, ...
         sprintf('packing: rho = p(end) = %.6f, t_f = p(end-1) = %.6f', ih.p(end), ih.p(end-1)));

% (2) Jacobian vs central finite differences (a handful of columns)
Rf = ih.residual;  p = ih.p;  n = numel(p);
cols = unique([1 4 7 8 n-1 n]);
errs = zeros(size(cols));
for k = 1:numel(cols)
    h = 1e-6*max(abs(p(cols(k))), 1);
    e = zeros(n,1);  e(cols(k)) = h;
    fd = (Rf(p + e) - Rf(p - e))/(2*h);
    errs(k) = norm(fd - ih.J(:, cols(k))) / max(norm(ih.J(:, cols(k))), 1);
end
ok = chk(ok, max(errs) < 1e-5, sprintf('Jacobian vs FD on columns %s: max rel err %.1e', mat2str(cols), max(errs)));

% (3) solve from a perturbed seed -> same trajectory (t_f), bounded rho
seedP = seedH;  seedP.Y(8:14, :) = seedP.Y(8:14, :)*(1 + 1e-3);  seedP.extra = a*(1 - 2e-3);
[zh, ihs] = ms_tfmin_hom(rv0, rvf, seedP, Tnd, cnd, mu, struct('tolR', 1e-10));
ok = chk(ok, ihs.converged && abs(zh(8) - z0(8)) < 1e-8, ...
         sprintf('solve from perturbed seed: converged = %d, |tf - tf_ref| = %.1e', ihs.converged, abs(zh(8) - z0(8))));
ok = chk(ok, abs(zh(1:7)'*zh(1:7) + ihs.rho^2 - 1) < 1e-10 && ihs.rho > 0, ...
         sprintf('on the sphere: rho^2 + |lam0|^2 - 1 = %.1e, rho = %.6f', zh(1:7)'*zh(1:7) + ihs.rho^2 - 1, ihs.rho));
% the de-normalised costates are the normal ones
ok = chk(ok, norm(zh(1:7)/ihs.rho - z8(1:7)) < 1e-6*norm(z8(1:7)), ...
         sprintf('lam/rho recovers the rho = 1 costates: rel err %.1e', norm(zh(1:7)/ihs.rho - z8(1:7))/norm(z8(1:7))));

if ok, fprintf('TEST_MS_TFMIN_HOM: ALL PASS\n');
else,  fprintf('TEST_MS_TFMIN_HOM: FAILURE (see lines above)\n');
end
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
