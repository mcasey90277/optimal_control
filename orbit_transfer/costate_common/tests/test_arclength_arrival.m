function ok = test_arclength_arrival()
%% Purpose:
%
%   Tests arclength_arrival -- the ARRIVAL-PHASE pseudo-arclength binding
%   of the min-time DRO -> tulip problem (the sheet's spine engine). Checks:
%
%   (1) the analytic parameter derivative. Arrival phase sA enters ONLY the
%       terminal state-matching rows, R_sA = [0; -x_A'(sA); 0; 0] with
%       x_A'(sA) = T_period * f_cr3bp(x_A(sA)); it must agree with a central
%       finite difference of the full residual in sA.
%   (2) the ballistic field is consistent with the propagated orbit.
%   (3) a SHORT arc from the certified 70 mN anchor toward larger sA
%       reproduces the independently measured walk (FINDINGS 36: t_f falls
%       from 17.798 d toward 17.42 d by sA = 0.0899) and stays certified.
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
tStar = 382981.289129055;

[B, anc] = arclength_arrival('setup');           % problem closures + anchor

% (4) a SECOND seed at another grid phase: the certified 26.44 d root at
%     cell (1,11), sA = 0.0754 + 10/12. Its .mat stores the root as
%     best.z / best.it.Y and carries no Tnd/cnd, so setup must take the
%     phase from opts and guard the operating point by the re-solve.
mat2 = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', 'mintime_70mN_certified.mat');
try
    [~, anc2] = arclength_arrival('setup', struct('sA0', 0.0754 + 10/12, 'anchorMat', mat2));
    ok2 = abs(anc2.p(anc2.ctf) - 5.963936) < 5e-3 && abs(anc2.sA - (0.0754 + 10/12)) < 1e-12;
    msg2 = sprintf('second seed set up at sA = %.4f: t_f = %.4f ND (stored 5.9639)', anc2.sA, anc2.p(anc2.ctf));
catch ME
    ok2 = false;  msg2 = ['second seed setup threw: ' ME.message];
end
ok = chk(ok, ok2, msg2);

% (2) ballistic field vs the propagated periodic orbit
s = 0.3137;
xA = B.stateA(s);
% central difference along the ORBIT (second order); a forward difference
% at h = 5e-4 ND carries an O(h) error of ~1e-3 and fails a 1e-6 test
h = 1e-6*B.tauA;
[~, Yp] = pumpkyn.cr3bp.prop(h, xA, B.mu);  [~, Ym] = pumpkyn.cr3bp.prop(-h, xA, B.mu);
fdx = (Yp(end,1:6)' - Ym(end,1:6)') / (2*h);
ok = chk(ok, norm(fdx - cr3bp_field(xA, B.mu)) < 1e-6*norm(fdx), ...
         sprintf('cr3bp_field matches the propagated orbit: rel err %.1e', norm(fdx - cr3bp_field(xA, B.mu))/norm(fdx)));

% (1) analytic R_sA vs central FD of the full residual
p0 = anc.p;  sA0 = anc.sA;
Ran = B.dRdq(p0, sA0);
h = 1e-6;
Rp = B.res(sA0 + h);  Rp = Rp(p0);
Rm = B.res(sA0 - h);  Rm = Rm(p0);
Rfd = (Rp - Rm)/(2*h);
ok = chk(ok, norm(Ran - Rfd) < 1e-5*norm(Rfd), ...
         sprintf('analytic R_sA vs FD: rel err %.1e (|R_sA| = %.3g)', norm(Ran - Rfd)/norm(Rfd), norm(Rfd)));
nz = find(abs(Ran) > 0);
ok = chk(ok, numel(nz) == 6 && all(nz == anc.termRows(1:6)'), ...
         sprintf('R_sA is nonzero ONLY on the six terminal state rows (%d nonzeros)', numel(nz)));

% (3) a short arc reproduces the measured walk and stays certified
A = arclength_arrival(anc, struct('direction', +1, 'sAStop', sA0 + 0.015, ...
                                  'nStep', 60, 'logFile', ''));
ok = chk(ok, numel(A.q) > 3 && A.q(end) >= sA0 + 0.014, ...
         sprintf('arc advanced from sA = %.4f to %.4f in %d roots', sA0, A.q(end), numel(A.q)));
tfEnd = A.p{end}(anc.ctf)*tStar/86400;
ok = chk(ok, tfEnd < 17.798 && tfEnd > 17.3, ...
         sprintf('t_f fell along the arc: %.4f d at sA = %.4f (walk gave 17.42 d at 0.0899)', tfEnd, A.q(end)));
ok = chk(ok, max(A.normR) < 1e-8, sprintf('every root converged: max |R| = %.1e', max(A.normR)));
rho = cellfun(@(pp) pp(end), A.p);
ok = chk(ok, all(rho > 0) && max(rho)/min(rho) < 3, ...
         sprintf('rho stays positive and bounded (%.4f .. %.4f): no normality loss', min(rho), max(rho)));

if ok, fprintf('TEST_ARCLENGTH_ARRIVAL: ALL PASS\n');
else,  fprintf('TEST_ARCLENGTH_ARRIVAL: FAILURE (see lines above)\n');
end
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
