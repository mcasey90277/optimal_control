function ok = test_certify_crossing()
%% Purpose:
%
%   Tests certify_crossing -- the gate stack applied to ONE arrival-phase
%   candidate (a homogeneous-chart ms root at a grid phase) before it may
%   enter the sheet: normal-chart polish, flown position AND velocity,
%   pumpkyn tfMin witness, conjugate test, min-time hypothesis gates. The
%   fixture is the certified 70 mN anchor itself (17.798 d at sA = 0.0754),
%   which must come back certified with its known numbers, and a corrupted
%   copy (costates scaled by 3) which must be REFUSED with a named reason.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

[B, anc] = arclength_arrival('setup');
C = certify_crossing(anc.p, anc.sA, B, anc);
ok = chk(ok, C.ok, sprintf('anchor certifies (reason: %s)', C.reason));
ok = chk(ok, abs(C.tfDays - 17.798) < 2e-3, sprintf('t_f = %.4f d (17.798)', C.tfDays));
ok = chk(ok, abs(C.dvKms - 0.7485) < 1e-3 && abs(C.mfKg - 12.20) < 0.02, ...
         sprintf('dV = %.4f km/s, fuel = %.2f kg (0.7485, 12.20)', C.dvKms, C.mfKg));
ok = chk(ok, C.flyKm < 1 && C.flyVms < 1, sprintf('flown miss %.3f km, %.3f m/s', C.flyKm, C.flyVms));
ok = chk(ok, C.dz < 1e-6, sprintf('tfMin witness |dz| = %.1e', C.dz));
ok = chk(ok, C.conj == 1 && C.g.dimS == 1 && C.g.minLamV > 0 && C.g.minQmt > 0, ...
         sprintf('conj %d, dim S %d, min|lam_v| %.2e, min Q %.2e', C.conj, C.g.dimS, C.g.minLamV, C.g.minQmt));

% the homogeneous -> normal chart conversion at a rho far from the anchor's.
% Scaling every multiplier (lam_0, the junction costates, rho) by a common
% factor gives the SAME extremal in a different place on the sphere, so the
% conversion must return the same z8. This exercises the division at
% rho = 0.5 rather than only at the anchor's 0.053, and it is where a wrong
% row range would show: the state rows must NOT be scaled.
s = 0.5/C.rho;
pS = anc.p;  K = anc.K;
pS(1:7) = s*pS(1:7);
Yj = reshape(pS(8:8+14*(K-1)-1), 14, K-1);
Yj(8:14,:) = s*Yj(8:14,:);
pS(8:8+14*(K-1)-1) = Yj(:);
pS(end) = s*pS(end);
Cs = certify_crossing(pS, anc.sA, B, anc);
ok = chk(ok, Cs.ok && abs(Cs.rho - 0.5) < 1e-12 && norm(Cs.z - C.z) < 1e-6*norm(C.z), ...
         sprintf('chart conversion at rho = %.3f returns the same extremal (|dz| = %.1e)', ...
                 Cs.rho, norm(Cs.z - C.z)));

% a corrupted candidate must be refused, not silently accepted
pBad = anc.p;  pBad(1:7) = 3*pBad(1:7);
Cb = certify_crossing(pBad, anc.sA, B, anc, struct('wallSec', 60));
ok = chk(ok, ~Cb.ok && ~isempty(Cb.reason), sprintf('corrupted candidate refused: %s', Cb.reason));

if ok, fprintf('TEST_CERTIFY_CROSSING: ALL PASS\n'); else, fprintf('TEST_CERTIFY_CROSSING: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
