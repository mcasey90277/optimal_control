% MINTIME_MS_GATE  Validate the min-time MS machinery before the real solve:
%   (1) single-shooting tulip min-time converges (seed source),
%   (2) MS seed has ~0 continuity,
%   (3) the analytic block Jacobian matches central finite differences.
% Small M=4 keeps the dense FD check cheap.

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
[rv0, rvf, P] = mintime_params();  Tmax = P.Tmax25;  c = P.c;  mu = P.muStar;
zSeed = [ 190.476497248065; -79.7064866984696; -0.430399154713168; ...
            0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
            4.32931089137559; 6.29081541876621];

[zt, rn, o] = mintime_solve(rv0, rvf, zSeed, Tmax, c, mu, 1500);
fprintf('single-shoot tulip min-time: ||R||=%.3e  nSwitch=%d  tf=%.4f\n', rn, o.nSwitch, zt(8));

M = 4;
[z, prob, si] = mintime_ms_seed(zt(1:7), zt(8), rv0, rvf, Tmax, c, mu, M);
fprintf('MS seed M=%d: maxCont=%.2e  termErr=%.2e\n', M, si.maxCont, si.termErr);

[R, J, info] = mintime_ms_residual(z, prob, true);
nZ = numel(z);  Jfd = zeros(numel(R), nZ);
for iZ = 1:nZ
    h = 1e-5 * max(1, abs(z(iZ)));
    zp = z;  zp(iZ) = zp(iZ) + h;
    zm = z;  zm(iZ) = zm(iZ) - h;
    Rp = mintime_ms_residual(zp, prob, false);
    Rm = mintime_ms_residual(zm, prob, false);
    Jfd(:,iZ) = (Rp - Rm) / (2*h);
end
Jf = full(J);
D = abs(Jf - Jfd);
relerr = norm(Jf - Jfd, 'fro') / max(norm(Jfd, 'fro'), 1);
[mx, lin] = max(D(:));  [ri, ci] = ind2sub(size(D), lin);
% RELATIVE per-entry error, ignoring tiny entries (FD noise floor)
scale = max(abs(Jf), abs(Jfd));  mask = scale > 1e-3;
relEntry = D ./ max(scale, 1e-12);
maxRelEntry = max(relEntry(mask));
fprintf('J vs FD: fro rel err = %.3e   max abs diff = %.2e at (r=%d,c=%d) |J|=%.2e |Jfd|=%.2e\n', ...
        relerr, mx, ri, ci, abs(Jf(ri,ci)), Jfd(ri,ci));
fprintf('         max RELATIVE per-entry error (|J|>1e-3) = %.3e\n', maxRelEntry);
fprintf('         tf-column (last col) rel err = %.3e\n', ...
        norm(Jf(:,nZ)-Jfd(:,nZ))/max(norm(Jfd(:,nZ)),1e-12));
fprintf('system size: nR=%d nZ=%d (square=%d)\n', numel(R), nZ, numel(R)==nZ);
if maxRelEntry < 1e-4, fprintf('GATE: PASS (analytic J == FD to FD-noise floor)\n'); else, fprintf('GATE: FAIL\n'); end
