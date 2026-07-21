% Z0_ACCEPT_75MN  THE Z0 acceptance test (pre-registered, Z0_BUILD.md par.8):
% Newton with the exact variational-STM Jacobian, from the banked P0i floor
% iterate (75 mN energy problem, ||R|| = 1.56e-2 where CS and FD Jacobians
% both failed). GATE: ||R|| <= 1e-8 within 30 iterations.
%
%   PASS -> the 75 mN LADDER ANCHOR is converged (first converged indirect
%           solve of the campaign at any multi-rev thrust); Z3 marches from it.
%   FAIL -> the wall is not derivative quality; pivot to direct-side seed.
%
% Requires: results/p0i_fd_finish.mat, Z0 gates green (test_ztl_z0).
% Output:   results/z0_anchor_75mN.mat (on PASS) / z0_accept_trace.mat (always)

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');

A = load(fullfile(resDir, 'p0i_fd_finish.mat'));  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;                              % 75 mN
tfL  = a.tf;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-12, 'odeAbsTol', 1e-14);

kMap = 2*Tmax/P0.c;                              % legacy -> BE costate map
lam  = kMap*a.lam0(:);

fprintf('=== Z0 ACCEPTANCE: variational-STM Newton @ 75 mN, eps=1 ===\n');

% Sanity gate A: mapped seed must reproduce the P0i floor (~1.5e-2)
[R, ~] = resjac(lam, false, rv0, rvf, tfL, P);
rn = norm(R);
fprintf('sanity: ||R(mapped seed)|| = %.4e  (P0i floor was 1.562e-2)\n', rn);
assert(rn > 1e-3 && rn < 1e-1, ...
    'mapped-seed residual far from the P0i floor -- map or EOM wrong, STOP');

% --- Newton with Armijo backtracking -----------------------------------------
maxIter = 30;
trace = nan(maxIter, 3);
tStart = tic;
for it = 1:maxIter
    [R, J] = resjac(lam, true, rv0, rvf, tfL, P);
    rn = norm(R);

    cs = max(abs(J), [], 1);  cs(cs == 0) = 1;
    dz = -(J ./ cs) \ R;
    dz = dz ./ cs(:);

    alpha = 1;  accepted = false;
    for ls = 1:14
        Rt = resjac(lam + alpha*dz, false, rv0, rvf, tfL, P);
        if all(isfinite(Rt)) && norm(Rt) < rn*(1 - 1e-4*alpha)
            accepted = true;  break
        end
        alpha = alpha/2;
    end
    if accepted
        lam = lam + alpha*dz;  rn = norm(Rt);
    end
    trace(it, :) = [rn, alpha*accepted, cond(J ./ cs)];
    fprintf('  it %2d: ||R|| = %.6e  alpha = %.3g  cond(Jeq) = %.2e  (%.0f s)\n', ...
            it, rn, alpha*accepted, cond(J ./ cs), toc(tStart));
    if rn < 1e-8, break; end
    if ~accepted
        fprintf('  Newton step rejected at all alpha -- stopping.\n');
        break
    end
end
save(fullfile(resDir, 'z0_accept_trace.mat'), 'trace', 'lam', 'rn');

% --- verdict + anchor product -------------------------------------------------
if rn < 1e-8
    o = ztl_flow([rv0(:); 1; lam], [0 tfL], P, false);
    uT = ramp_u(o.y, P);
    mF = o.yf(7);
    anchor = struct('Tmax_mN', 75, 'tf', tfL, 'eps', 1, 'lam0_BE', lam, ...
        'resNorm', rn, 'iters', it, 'kMap', kMap, ...
        'mProp_kg', P0.m0kg*(1-mF), 'dV_kms', P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'uMin', min(uT), 'uMax', max(uT), 'nSegs', o.nSegs, ...
        'rv0', rv0, 'rvf', rvf, 'P', P, 'P0', P0);
    save(fullfile(resDir, 'z0_anchor_75mN.mat'), 'anchor');
    fprintf(['\nGATE Z0: PASS -- LADDER ANCHOR CONVERGED @ 75 mN in %d iters.\n' ...
             '  ||R|| = %.3e  prop = %.4f kg  dV = %.4f km/s  u in [%.3f, %.3f]\n' ...
             '  saved %s\n'], it, rn, anchor.mProp_kg, anchor.dV_kms, ...
            anchor.uMin, anchor.uMax, fullfile(resDir, 'z0_anchor_75mN.mat'));
else
    fprintf('\nGATE Z0: FAIL -- floors at %.3e with the EXACT variational J.\n', rn);
    fprintf('The wall is not derivative quality; pivot per Z0_BUILD.md par.8.\n');
end

% ---------------------------------------------------------------------------
function [R, J, o] = resjac(lamv, withJ, rv0, rvf, tfL, P)
% Shooting residual [rv(tf)-rvf; lam_m(tf)] and (optionally) its exact
% Jacobian from the variational STM's costate columns.
o = ztl_flow([rv0(:); 1; lamv], [0 tfL], P, withJ);
R = [o.yf(1:6) - rvf(:); o.yf(14)];
if withJ
    J = [o.PHI(1:6, 8:14); o.PHI(14, 8:14)];
else
    J = [];
end
end

function u = ramp_u(yGrid, P)
% Clamped ramp throttle along a stored trajectory grid [Nx14].
lamv = sqrt(sum(yGrid(:, 11:13).^2, 2));
S = 1 - lamv*P.c./yGrid(:, 7) - yGrid(:, 14);
u = min(max(0.5 - S/(2*P.eps), 0), 1);
end
