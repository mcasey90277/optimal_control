function traj = ms_traj(Z, prob)
% MS_TRAJ  Stitched trajectory + PMP diagnostics from an MS solution.
%
% Integrates every arc from its joint value and concatenates. Joint
% mismatches of a CONVERGED solution are ~1e-12; maxJointDefect reports
% them so callers can verify. S is the min-fuel switching function
% S = 1 - ||lambda_v||c/m - lambda_m; u the smoothed throttle at prob.epsSmooth.
%
% INPUTS:
%   Z    - unknowns [(14M-7)x1]
%   prob - problem struct with tJ [1x(M+1)]
%
% OUTPUTS:
%   traj - struct: t [1xL], Y [14xL], S [1xL], u [1xL], mf, dV_kms,
%          prop_kg, switches (S zero crossings), bangFrac (fraction of
%          time points with u>0.95 or u<0.05), maxJointDefect

M = numel(prob.tJ) - 1;
[~, yJ] = ms_unpack(Z, prob);
t = []; Y = []; maxJD = 0;
for k = 1:M
    [tk, Yk] = ode113(@(tt, y) lt_pmp_eom_minfuel(tt, y, prob.Tmax, prob.c, ...
                      prob.muStar, prob.epsSmooth), ...
                      [prob.tJ(k) prob.tJ(k+1)], yJ(:, k), prob.odeOpts);
    if k < M
        maxJD = max(maxJD, max(abs(Yk(end, :).' - yJ(:, k+1))));
        t = [t, tk(1:end-1).'];  Y = [Y, Yk(1:end-1, :).'];  %#ok<AGROW>
    else
        t = [t, tk.'];           Y = [Y, Yk.'];              %#ok<AGROW>
    end
end

nLamV = sqrt(sum(Y(11:13, :).^2, 1));
S     = 1 - nLamV.*prob.c./Y(7, :) - Y(14, :);
u     = (1 - tanh(S/(2*prob.epsSmooth)))/2;
mf    = Y(7, end);

traj = struct('t', t, 'Y', Y, 'S', S, 'u', u, 'mf', mf, ...
    'dV_kms',  prob.c*log(1/mf)*prob.p.lStar/prob.p.tStar, ...
    'prop_kg', prob.p.m0kg*(1 - mf), ...
    'switches', nnz(diff(sign(S)) ~= 0), ...
    'bangFrac', mean(u > 0.95 | u < 0.05), ...
    'maxJointDefect', maxJD);
end
