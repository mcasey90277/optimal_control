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
%     4. THE COLLAPSE CONTRACT, by identifier: a blown-up flow throws
%        cartpole_pmp_prop:collapse on BOTH paths (with and without the STM);
%        a coding error keeps its OWN identifier and is never relabelled (so
%        the engine can tell a rejected iterate from a bug); malformed
%        requests are refused as :input, even at dt = 0;
%     5. with no STM requested the state still matches the STM run, and PHI
%        comes back empty;
%     6. SHORT-TIME ORACLE: over a tiny step the flow reproduces the field.
%        The semigroup checks alone would pass for a consistently wrong flow
%        (even the identity map); this one ties the flow to cartpole_pmp_rhs.
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
ok = chk(ok, isequal(size(PHI), [8 8]) && isequal(size(yE), [8 1]), 'shapes: yEnd [8 x 1] (a COLUMN), PHI [8 x 8]');

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

%% The collapse contract, by IDENTIFIER and on BOTH paths. "Any exception"
%  is not a pass: a coding bug throws too, and the point of the contract is
%  that the engine can tell the two apart.
%  The fixture is a state whose DERIVATIVE overflows on the first evaluation
%  (q2dot = 1e200 is finite; its square is not), so the integrator cannot
%  take a step and the detector fires in ~0.05 s on both paths. A merely
%  "large costate" state is NOT usable here: under the corrected dynamics
%  (2026-09-17) such a flow oscillates rather than diverging, and ode113
%  grinds for minutes at tiny steps -- a unit test must not depend on how
%  long an integrator takes to give up.
yBlow = [0; 0; 0; 1e200; 1; 1; 1; 1];
warnWas = warning('off', 'all');  restoreW = onCleanup(@() warning(warnWas));
idNo  = thrownId(@() cartpole_pmp_prop(1, yBlow, false, p));
idStm = thrownId(@() cartpole_pmp_prop(1, yBlow, true,  p));
ok = chk(ok, strcmp(idNo, 'cartpole_pmp_prop:collapse'), ...
         sprintf('a blown-up flow throws the COLLAPSE identifier, state-only path (%s)', idNo));
ok = chk(ok, strcmp(idStm, 'cartpole_pmp_prop:collapse'), ...
         sprintf('and on the STM path (%s)', idStm));

% a coding error must NOT be relabelled: a parameter struct with no .g
pBad = rmfield(p, 'g');
idBug = thrownId(@() cartpole_pmp_prop(dt, y0, false, pBad));
ok = chk(ok, ~isempty(idBug) && ~strcmp(idBug, 'cartpole_pmp_prop:collapse'), ...
         sprintf('a coding error keeps its own identifier (%s), not :collapse', idBug));

% malformed requests are refused by name, before the dt == 0 shortcut
ok = chk(ok, strcmp(thrownId(@() cartpole_pmp_prop(0, [y0(1:7); NaN], false, p)), 'cartpole_pmp_prop:input') ...
          && strcmp(thrownId(@() cartpole_pmp_prop(0, y0(1:7), false, p)), 'cartpole_pmp_prop:input') ...
          && strcmp(thrownId(@() cartpole_pmp_prop(dt, y0 + 1i*1e-20, true, p)), 'cartpole_pmp_prop:input'), ...
         'non-finite, mis-sized and complex y0 are refused as :input (even at dt = 0)');

% the no-STM path returns an EMPTY second output, as the contract says
[~, PHIno] = cartpole_pmp_prop(dt, y0, false, p);
ok = chk(ok, isempty(PHIno), 'needSTM = false returns PHI = []');

% short-time oracle: over a tiny step the flow IS the field (independent of
% the semigroup checks above, which a consistently wrong flow would satisfy)
hS = 1e-5;
yS = cartpole_pmp_prop(hS, y0, false, p);
fS = cartpole_pmp_rhs(y0, p);
ok = chk(ok, max(abs((yS - y0)/hS - fS)) < 1e-3*max(1, max(abs(fS))), ...
         sprintf('short-time flow matches the field: %.1e', max(abs((yS - y0)/hS - fS))));

if ok, fprintf('TEST_CARTPOLE_PMP_PROP: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP_PROP: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function id = thrownId(fn)
%% Purpose:
%
%   The identifier of the error fn() throws, or '' if it returns normally.
%
id = '';
try
    fn();
catch err
    id = err.identifier;
    if isempty(id), id = '(no identifier)'; end
end
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
