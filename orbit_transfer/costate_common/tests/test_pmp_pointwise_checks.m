function ok = test_pmp_pointwise_checks()
%% Purpose:
%
%   Tests pmp_pointwise_checks -- the first-order checks that live on the
%   FLIGHT rather than in the shooting residual: the Hamiltonian, the
%   adjoint equations against a finite-difference dH/dx, the EXACT
%   minimum-principle gap of the control the propagator actually applied
%   (not a sampled restatement of the analytic minimiser), the weak
%   throttle condition Q_mt >= 0, and transversality lambda_m(t_f) = 0.
%   A flight with the costate sign corrupted must show a positive gap.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  z = E.z(:);
[t, Y] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6); 1; z(1:7)], B.Tnd, B.cnd, B.mu);

P = pmp_pointwise_checks(t, Y, B.Tnd, B.cnd, B.mu);
ok = chk(ok, P.Hmax < 1e-6, sprintf('Hamiltonian max|H| = %.2e', P.Hmax));
ok = chk(ok, P.adjErr < 1e-7, sprintf('adjoint equations, relative error %.2e', P.adjErr));
ok = chk(ok, P.fdAgree < 1e-5, sprintf('two FD step sizes agree to %.2e', P.fdAgree));
ok = chk(ok, P.dirGap <= 1e-12, sprintf('exact minimum-principle gap of the APPLIED control %.2e', P.dirGap));
ok = chk(ok, isfield(P, 'fullGap') && P.fullGap <= 1e-12, ...
         sprintf('the FULL control gap (direction + throttle) %.2e', P.fullGap));
ok = chk(ok, isfield(P, 'gapMin') && isfield(P, 'gapMax') && P.gapMin <= P.gapMax && abs(P.gapMin) <= 1e-12, ...
         sprintf('signed extremes kept: %+.2e .. %+.2e', P.gapMin, P.gapMax));
ok = chk(ok, P.throttleErr < 1e-10, sprintf('applied throttle is 1 to %.2e', P.throttleErr));
ok = chk(ok, isfield(P, 'throttleMassErr') && P.throttleMassErr < 1e-10 && P.throttleAccErr < 1e-10, ...
         sprintf('on BOTH rows: acceleration %.2e, mass %.2e', P.throttleAccErr, P.throttleMassErr));
ok = chk(ok, P.minQmt >= 0, sprintf('weak throttle condition min Q_mt = %.3f >= 0', P.minQmt));
ok = chk(ok, P.lamMf < 1e-6, sprintf('transversality |lambda_m(t_f)| = %.2e', P.lamMf));
ok = chk(ok, numel(P.tSample) == P.nSample && P.tSample(end) <= t(end), ...
         sprintf('%d samples, coverage reported as t/t_f in [%.3f, %.3f]', ...
                 P.nSample, P.tSample(1)/t(end), P.tSample(end)/t(end)));

% MUTATION: a vector field that thrusts AGAINST the primer. The check must
% see it -- this is the bug class N6 exists for (a sign error in the
% propagator's control law), and it is invisible to the shooting residual,
% which would happily converge to the wrong problem.
Pc = pmp_pointwise_checks(t, Y, B.Tnd, B.cnd, B.mu, struct('rhs', @wrongSignField));
ok = chk(ok, Pc.dirGap > 1e-2, sprintf('a wrong-sign thrust field opens the gap: %.2e', Pc.dirGap));
ok = chk(ok, Pc.adjErr > 1e-3, sprintf('and breaks the adjoint check too: %.2e', Pc.adjErr));

% MUTATION 2 (Astra review 2026-09-11): a field whose MASS row burns at half
% throttle while the acceleration is correct. Invisible to the direction gap
% and to the acceleration-side throttle; the mass-row throttle must see it.
Pm = pmp_pointwise_checks(t, Y, B.Tnd, B.cnd, B.mu, struct('rhs', @halfMassFlowField));
ok = chk(ok, Pm.dirGap <= 1e-12 && Pm.throttleAccErr < 1e-10, ...
         sprintf('half mass flow: direction gap %.1e and acceleration throttle %.1e stay clean', ...
                 Pm.dirGap, Pm.throttleAccErr));
ok = chk(ok, Pm.throttleMassErr > 0.49 && Pm.throttleErr > 0.49, ...
         sprintf('but the MASS-row throttle reads it: error %.3f', Pm.throttleMassErr));
ok = chk(ok, Pm.fieldGap > 1e-3 && Pm.fullGap <= 1e-12, ...
         sprintf('and the FIELD''s own gap opens (%.2e) while the reconstructed one stays closed (%.1e)', Pm.fieldGap, Pm.fullGap));

% MUTATION 3: an OVER-UNIT thrust along the minimiser. The old gap started
% at zero and took a max, so the negative gap this produces was clipped to
% a perfect zero; the signed minimum and |full gap| must both show it.
Po = pmp_pointwise_checks(t, Y, B.Tnd, B.cnd, B.mu, struct('rhs', @overThrustField));
ok = chk(ok, Po.gapMin < -1e-6, sprintf('over-unit thrust: signed minimum gap %.2e < 0 (was clipped)', Po.gapMin));
ok = chk(ok, Po.fullGap > 1e-6 && Po.throttleAccErr > 1e-4, ...
         sprintf('|full gap| %.2e and acceleration throttle error %.1e flag it', Po.fullGap, Po.throttleAccErr));

if ok, fprintf('TEST_PMP_POINTWISE_CHECKS: ALL PASS\n'); else, fprintf('TEST_PMP_POINTWISE_CHECKS: FAIL\n'); end
end

function F = wrongSignField(y, Tmax, c, mu)
% WRONGSIGNFIELD  The min-time field with the thrust direction flipped.
% INPUTS: y; Tmax; c; mu.  OUTPUTS: F [14x1].
F  = mintime_rhs_point(y, Tmax, c, mu);
F0 = mintime_rhs_point(y, 0, c, mu);
F(4:6) = 2*F0(4:6) - F(4:6);
end

function F = halfMassFlowField(y, Tmax, c, mu)
% HALFMASSFLOWFIELD  Correct acceleration, mass row at half throttle.
% INPUTS: y; Tmax; c; mu.  OUTPUTS: F [14x1].
F = mintime_rhs_point(y, Tmax, c, mu);
F(7) = 0.5*F(7);
end

function F = overThrustField(y, Tmax, c, mu)
% OVERTHRUSTFIELD  Thrust acceleration scaled by 1.001 along the minimiser.
% INPUTS: y; Tmax; c; mu.  OUTPUTS: F [14x1].
F  = mintime_rhs_point(y, Tmax, c, mu);
F0 = mintime_rhs_point(y, 0, c, mu);
F(4:6) = F0(4:6) + 1.001*(F(4:6) - F0(4:6));
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
