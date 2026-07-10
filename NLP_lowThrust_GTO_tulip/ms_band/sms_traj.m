function traj = sms_traj(Z, prob)
% SMS_TRAJ  Stitched trajectory + PMP diagnostics from a Sundman-MS solution.
%
% Integrates every arc from its joint value in sigma and concatenates.
% Physical time is the carried state t = Y(8,:), so all diagnostics/plots
% can stay in physical time. S is the min-fuel switching function
% S = 1 - ||lamV||c/m - lamM; u the smoothed throttle at prob.epsSmooth
% (pointwise logic unchanged from MS_TRAJ). maxJointDefect reports joint
% mismatches (~1e-12 for a converged solution).
%
% INPUTS:
%   Z    - unknowns [(16M-8)x1]
%   prob - problem struct with sJ [1x(M+1)]
%
% OUTPUTS:
%   traj - struct: sig [1xL], t [1xL] (physical time = Y(8,:)), Y [16xL],
%          S [1xL], u [1xL], mf, dV_kms, prop_kg, switches (S zero
%          crossings), bangFrac (fraction of samples with u>0.95 or
%          u<0.05), maxJointDefect

M = numel(prob.sJ) - 1;
[~, yJ] = sms_unpack(Z, prob);
sig = []; Y = []; maxJD = 0;
for k = 1:M
    [sk, Yk] = ode113(@(ss, y) sms_eom(ss, y, prob.Tmax, prob.c, ...
                      prob.muStar, prob.epsSmooth, prob.pSund), ...
                      [prob.sJ(k) prob.sJ(k+1)], yJ(:, k), prob.odeOpts);
    if k < M
        maxJD = max(maxJD, max(abs(Yk(end, :).' - yJ(:, k+1))));
        sig = [sig, sk(1:end-1).'];  Y = [Y, Yk(1:end-1, :).'];  %#ok<AGROW>
    else
        sig = [sig, sk.'];           Y = [Y, Yk.'];              %#ok<AGROW>
    end
end

t     = Y(8, :);
nLamV = sqrt(sum(Y(12:14, :).^2, 1));
S     = 1 - nLamV.*prob.c./Y(7, :) - Y(15, :);
u     = (1 - tanh(S/(2*prob.epsSmooth)))/2;
mf    = Y(7, end);

traj = struct('sig', sig, 't', t, 'Y', Y, 'S', S, 'u', u, 'mf', mf, ...
    'dV_kms',  prob.c*log(1/mf)*prob.p.lStar/prob.p.tStar, ...
    'prop_kg', prob.p.m0kg*(1 - mf), ...
    'switches', nnz(diff(sign(S)) ~= 0), ...
    'bangFrac', mean(u > 0.95 | u < 0.05), ...
    'maxJointDefect', maxJD);
end
