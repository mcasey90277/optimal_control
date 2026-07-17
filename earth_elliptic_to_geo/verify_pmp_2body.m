function ver = verify_pmp_2body(out, par)
% VERIFY_PMP_2BODY  First-order PMP consistency from the NLP's KKT duals.
%
% Costates = defect-constraint multipliers out.lamDef (discrete adjoint up to a
% positive scale and a global sign). Checks:
%   (1) primer direction: thrust || -lam_v on burn arcs (scale-invariant);
%   (2) switching sign law via the empirical-beta route: the min-fuel switching
%       function S = 1 - beta*W with W = (c/m)*||lam_v|| + lam_m; beta pinned
%       as the ROBUST MEDIAN of 1/W over switch-adjacent intervals (absorbs the
%       covector scale); require S<0 on >=98% of burn nodes, S>0 on >=98% of
%       coast nodes;
%   (3) relative transversality |lam_m(end)| / max|lam_m| <= 1e-3 (free mass).
%
% INPUTS:
%   out - casadi_lt_2body result struct [fields .U .X .lamDef]
%   par - params struct (kepler_lt_params output) [scalar struct, uses .c]
%
% OUTPUTS:
%   ver - struct with fields:
%           .primerMeanDeg  - mean primer-vector misalignment on burns [deg]
%           .beta           - robust-median switching-function scale [scalar]
%           .betaSpreadPct  - MAD spread of beta candidates [%]
%           .burnSignPct    - pct of burn nodes with S<0 [%]
%           .coastSignPct   - pct of coast nodes with S>0 [%]
%           .lamMendRel     - |lam_m(end)| / max|lam_m|, relative transversality
%           .pass           - true if all gates satisfied [logical]
%
% REFERENCES:
%   [1] NLP_lowThrust_GTO_tulip/sundman_minfuel/verify_tf_front.m (empirical-beta).
%   [2] HONEST_EVALUATION_DV_TF_FRONT.md (robust-beta lesson; relative gate).
N    = size(out.lamDef, 2);
ss   = out.U(4, 1:N);
burn = ss > 0.5;
lamV = out.lamDef(4:6, :);  lamM = out.lamDef(7, :);
mAvg = 0.5*(out.X(7,1:N) + out.X(7,2:N+1));
% global sign: primer alignment must be < 90 deg on burns
lVn = sqrt(sum(lamV.^2, 1));
cosb = zeros(1, N);
for k = 1:N
    cosb(k) = dot(out.U(1:3,k), -lamV(:,k)) / max(lVn(k), 1e-30);
end
if mean(cosb(burn)) < 0, lamV = -lamV; lamM = -lamM; cosb = -cosb; end
ver.primerMeanDeg = mean(real(acosd(max(-1, min(1, cosb(burn))))));
W = (par.c ./ mAvg) .* sqrt(sum(lamV.^2,1)) + lamM;
sw = find(diff(burn) ~= 0);                     % switch-adjacent intervals
bcand = 1 ./ W([sw, min(sw+1, N)]);
ver.beta = median(bcand);
ver.betaSpreadPct = 100 * median(abs(bcand - ver.beta)) / abs(ver.beta);  % MAD, no toolbox
S = 1 - ver.beta * W;
ver.burnSignPct  = 100 * mean(S(burn)  < 0);
ver.coastSignPct = 100 * mean(S(~burn) > 0);
ver.lamMendRel   = abs(lamM(end)) / max(abs(lamM));
ver.pass = ver.primerMeanDeg < 1.0 && ver.burnSignPct >= 98 && ...
           ver.coastSignPct >= 98 && ver.lamMendRel <= 1e-3;
fprintf(['verify_pmp_2body: primer %.3f deg | burn %.1f%% coast %.1f%% ' ...
         '(beta spread %.1f%%) | lamM rel %.1e | pass=%d\n'], ver.primerMeanDeg, ...
         ver.burnSignPct, ver.coastSignPct, ver.betaSpreadPct, ver.lamMendRel, ver.pass);
end
