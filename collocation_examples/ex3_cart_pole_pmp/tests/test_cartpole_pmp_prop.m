function ok = test_cartpole_pmp_prop()
%% Purpose:
%
%   The propagator and the STM it hands the shooting engine:
%     1. the STM matches central finite differences of the propagated map,
%        column by column (the STM is what the Newton step is built from, so
%        a wrong one shows up as slow convergence, not as a wrong answer);
%     2. PHI(0) = I;
%     3. the group property: propagating dt in one call equals two calls of
%        dt/2, and the STMs multiply;
%     4. it THROWS rather than returning junk when the flow blows up (the
%        engine's contract);
%     5. with no STM requested the state still matches the STM run.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
y0 = [0; 0; 0; 0; 0.8; -1.5; 0.4; 0.2];
dt = 0.6;

[yE, PHI] = cartpole_pmp_prop(dt, y0, true, p);
ok = chk(ok, isequal(size(PHI), [8 8]) && numel(yE) == 8, 'shapes: yEnd [8 x 1], PHI [8 x 8]');

PHIfd = zeros(8);
h = 1e-6;
for col = 1:8
    yp = y0;  yp(col) = yp(col) + h;
    ym = y0;  ym(col) = ym(col) - h;
    PHIfd(:,col) = (cartpole_pmp_prop(dt, yp, false, p) - cartpole_pmp_prop(dt, ym, false, p))/(2*h);
end
rel = max(abs(PHI(:) - PHIfd(:)))/max(abs(PHIfd(:)));
ok = chk(ok, rel < 1e-6, sprintf('STM matches finite differences: relative %.1e', rel));

[~, PHI0] = cartpole_pmp_prop(0, y0, true, p);
ok = chk(ok, max(abs(PHI0(:) - reshape(eye(8), 64, 1))) < 1e-12, 'PHI(0) = I');

[yH, PHIa] = cartpole_pmp_prop(dt/2, y0, true, p);
[yF, PHIb] = cartpole_pmp_prop(dt/2, yH, true, p);
ok = chk(ok, max(abs(yF - yE)) < 1e-9, sprintf('one step equals two half steps (%.1e)', max(abs(yF - yE))));
ok = chk(ok, max(abs(reshape(PHIb*PHIa - PHI, 64, 1))) < 1e-7, 'and their STMs multiply');

yNo = cartpole_pmp_prop(dt, y0, false, p);
ok = chk(ok, max(abs(yNo - yE)) < 1e-9, 'the state is the same with and without the STM');

threw = false;
try
    cartpole_pmp_prop(50, [0; 0; 0; 0; 1e6; 1e6; 1e6; 1e6], false, p);
catch
    threw = true;
end
ok = chk(ok, threw, 'a blown-up flow THROWS rather than returning junk');

if ok, fprintf('TEST_CARTPOLE_PMP_PROP: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP_PROP: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
