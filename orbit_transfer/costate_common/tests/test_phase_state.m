function ok = test_phase_state()
%% Purpose:
%
%   Tests phase_state -- one periodic orbit in, two closures out: the state
%   at a PHASE FRACTION and its derivative with respect to that fraction.
%   It is what every campaign writes by hand today (a dozen interp1
%   'spline' sites, plus the private copies in transfer_study and
%   arclength_arrival).
%
%   Checks:
%     1. wrapping: s, s+1 and s-1 give the SAME state, bitwise;
%     2. shape: a column, one per state component;
%     3. the derivative is with respect to the FRACTION (it carries the
%        period), verified against a central difference;
%     4. EQUIVALENCE GATE: on the real tau = 1 DRO and 7-petal tulip it
%        reproduces, BITWISE, what arclength_arrival's private makePP /
%        ppDer produced before the move (values captured 2026-09-11, at
%        four phases including one next to the seam).
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
addpath(here);

th = linspace(0, 2*pi, 61).';
[x, seam, dx] = phase_state(th, [cos(th), sin(th)]);

ok = chk(ok, isequal(x(0.25), x(1.25)) && isequal(x(0.25), x(-0.75)), ...
         'the phase wraps: s, s+1 and s-1 give the same state bitwise');
ok = chk(ok, isequal(size(x(0.3)), [2 1]), sprintf('a column comes back: %s', mat2str(size(x(0.3)))));
ok = chk(ok, seam.value < 1e-12 && seam.deriv < 1e-12, ...
         sprintf('the seam comes with it: value %.1e, derivative %.1e', seam.value, seam.deriv));
h = 1e-6;  s0 = 0.317;
fd = (x(s0+h) - x(s0-h))/(2*h);
ok = chk(ok, norm(dx(s0) - fd) < 1e-6*norm(fd), ...
         sprintf('the derivative is per unit PHASE: rel err %.1e vs a central difference', norm(dx(s0) - fd)/norm(fd)));

% ---- 4. equivalence with the pre-move implementations -------------------
muStar = 0.012150585609624;
[tD, rvD] = get_family_orbit('dro',   struct('tau', 1.0, 'muStar', muStar));
[tT, rvT] = get_family_orbit('tulip', struct('Np', 7, 'pm', -1, 'muStar', muStar));
[xD, ~, ~]   = phase_state(tD, rvD);
[xA, ~, dxA] = phase_state(tT, rvT);

gold = struct( ...
 'stateA_00754', [0.97717298240350747 -0.0018391491735214541 0.012893090747814712 -0.41929028240655702 -0.97510584090005825 -0.28545908666055719].', ...
 'stateA_09991', [1.0437913453349701 -0.00062868892698376782 -0.079852201373800211 0.0014185191390348618 0.13340529070816065 -0.005245501582510477].', ...
 'dstateA_03137', [0.88959023401290604 0.19607095692062596 1.0682308179699813 -1.1140335301783042 4.056725786950329 8.6827746724276036].', ...
 'stateD_09991', [0.91452677659067028 -0.00044294528436846859 1.1432921558637711e-33 -0.0027076148153366841 0.49215585906683712 2.6749347279890007e-33].');

pairs = {xA(0.0754), gold.stateA_00754, 'arrival state at 0.0754'; ...
         xA(0.9991), gold.stateA_09991, 'arrival state at 0.9991 (next to the seam)'; ...
         dxA(0.3137), gold.dstateA_03137, 'arrival phase DERIVATIVE at 0.3137'; ...
         xD(0.9991), gold.stateD_09991, 'departure state at 0.9991'};
for k = 1:size(pairs, 1)
    d = max(abs(pairs{k,1} - pairs{k,2}));
    ok = chk(ok, isequal(pairs{k,1}, pairs{k,2}), ...
             sprintf('%s reproduced bitwise (max diff %.1e)', pairs{k,3}, d));
end

if ok, fprintf('TEST_PHASE_STATE: ALL PASS\n'); else, fprintf('TEST_PHASE_STATE: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
