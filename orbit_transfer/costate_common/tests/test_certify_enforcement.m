function ok = test_certify_enforcement()
%% Purpose:
%
%   END-TO-END ENFORCEMENT test of certify_root: every quantity the
%   certifier was widened to require on 2026-09-11 must, when out of
%   tolerance, turn C.ok false with a reason that NAMES it. Instrument-level
%   mutation tests (test_pmp_pointwise_checks and friends) prove the
%   instruments can see a defect; they do not prove the caller reads them.
%   Astra review #2: "the production certifier still enforces less than the
%   study" -- the computed-not-enforced defect class of FINDINGS 36/41/47.
%
%   Two kinds of mutation:
%     - TOLERANCE SQUEEZE: the anchor is left alone and one tolerance is set
%       below its measured value -- proves the gate is WIRED to the number;
%     - FIELD INJECTION: a wrong field is handed to one instrument through
%       the pwOpts / gatesOpts hooks -- proves the number DETECTS a defect
%       and the certifier refuses on it.
%   The anchor must certify untouched.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);  z = E.z(:);
seed = seed_from_z8(z, rv0(1:6), 24, B.Tnd, B.cnd, B.mu);
base = struct('sA', E.sA, 'sD', E.sD);

C0 = certify_root(seed, rv0(1:6), rvf(1:6), B, base);
ok = chk(ok, C0.ok, sprintf('the anchor certifies untouched (%s)', C0.reason));
ok = chk(ok, all(isfinite([C0.fullGap C0.throttleErr C0.fieldErr C0.adjErrRef C0.nullResid C0.Hresid C0.liftMargin])), ...
         sprintf('and carries the new numbers: |full gap| %.1e, throttle %.1e, X2 %.1e/%.1e, lift %.1e, |H+1| %.1e, margin %.0fx', ...
                 C0.fullGap, C0.throttleErr, C0.fieldErr, C0.adjErrRef, C0.nullResid, C0.Hresid, C0.liftMargin));
if isstruct(C0.conjDense), D0 = C0.conjDense; else, D0 = struct('clear', false, 'nUnresolved', NaN, 'nNearMiss', NaN, 'nZero', NaN); end
ok = chk(ok, isstruct(C0.conjDense) && D0.clear && D0.nUnresolved == 0, ...
         sprintf('and the dense scan: %d near-miss cleared, %d unresolved, %d zero', ...
                 D0.nNearMiss, D0.nUnresolved, D0.nZero));

% ---- tolerance squeezes: each gate is WIRED ------------------------------
sq = {'tolGap',        1e-20, 'full gap';
      'tolThrottle',   1e-20, 'throttle';
      'tolField',      1e-20, 'X2';
      'tolLift',       1e-12, 'backward error';
      'liftMarginMin', 1e6,   'lift_margin'};
for k = 1:size(sq, 1)
    o = base;  o.(sq{k,1}) = sq{k,2};
    Ck = certify_root(seed, rv0(1:6), rvf(1:6), B, o);
    ok = chk(ok, ~Ck.ok && contains(lower(Ck.reason), lower(sq{k,3})), ...
             sprintf('%s = %g refuses, naming it: %s', sq{k,1}, sq{k,2}, Ck.reason));
end

% ---- field injections: each number DETECTS ---------------------------------
% (a) a wrong-sign thrust field handed to the pointwise checks: the
%     Hamiltonian is the first gate to see it, and the certifier must refuse
o = base;  o.pwOpts = struct('rhs', @wrongSignField);
Ca = certify_root(seed, rv0(1:6), rvf(1:6), B, o);
ok = chk(ok, ~Ca.ok && Ca.fullGap > 1e-2, ...
         sprintf('wrong-sign field in the pointwise checks: refused (%s); |full gap| stored %.2e', Ca.reason, Ca.fullGap));
% (b) half mass flow: acceleration rows clean, mass row at u = 0.5
o = base;  o.pwOpts = struct('rhs', @halfMassFlowField);
Cb = certify_root(seed, rv0(1:6), rvf(1:6), B, o);
ok = chk(ok, ~Cb.ok && Cb.throttleErr > 0.49, ...
         sprintf('half mass flow: refused (%s); throttle error stored %.3f', Cb.reason, Cb.throttleErr));
% (c) a wrong-sign field handed to the GATES' X2 comparison only: nothing
%     else changes, so the refusal must be X2's
o = base;  o.gatesOpts = struct('rhs', @wrongSignField);
Cc = certify_root(seed, rv0(1:6), rvf(1:6), B, o);
ok = chk(ok, ~Cc.ok && contains(Cc.reason, 'X2') && Cc.fieldErr > 1e-3, ...
         sprintf('wrong field in X2 only: refused by X2 (%s)', Cc.reason));
% (d) the dense scan with an impossible clear factor: the anchor's endpoint
%     near-miss becomes UNRESOLVED, and unresolved must block
o = base;  o.conjOpts = struct('clearFactor', 1e12);
Cd = certify_root(seed, rv0(1:6), rvf(1:6), B, o);
ok = chk(ok, ~Cd.ok && contains(Cd.reason, 'UNRESOLVED') && isstruct(Cd.conjDense) && Cd.conjDense.nUnresolved >= 1, ...
         sprintf('an unresolved scan candidate blocks: %s', Cd.reason));
% (e) the dense scan switched off is honoured but visible
o = base;  o.conjSpectrum = false;
Ce = certify_root(seed, rv0(1:6), rvf(1:6), B, o);
ok = chk(ok, Ce.ok && isempty(Ce.conjDense), 'conjSpectrum = false: certifies, and conjDense is empty (not a silent pass)');

if ok, fprintf('TEST_CERTIFY_ENFORCEMENT: ALL PASS\n'); else, fprintf('TEST_CERTIFY_ENFORCEMENT: FAIL\n'); end
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

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
