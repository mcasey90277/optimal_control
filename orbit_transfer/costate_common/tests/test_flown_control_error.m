function ok = test_flown_control_error()
%% Purpose:
%
%   Tests flown_control_error (the G1b gate) and the two helpers it flies
%   through, ctrl_quad and cr3bp_thrust_rhs, on a synthetic "direct
%   solution" whose nodes are the TRUE flight of a known control, so the
%   right answer is zero error:
%     - constant thrust direction, throttle th(t) = 0.3 + 0.8 t - 0.6 t^2
%       (a global quadratic, so the Hermite-Simpson node/midpoint/node
%       reconstruction is EXACT on every interval, and a linear one is not).
%   Checks:
%     1. ctrl_quad reproduces a quadratic exactly at the node, midpoint,
%        node and an interior point.
%     2. EXACT NODES + midpoint controls: erNd, evNd at integrator level.
%     3. INJECTED MISS: displacing the terminal node by d in position (and
%        separately in velocity) reads back as d.
%     4. RECONSTRUCTION MATTERS: dropping the midpoint controls (linear
%        interpolation of the same quadratic throttle) gives a clearly
%        larger error -- the gate sees a control that was not flown.
%     5. TIME GRID: tNodes = [] with s * tf gives the same numbers.
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

%% 1. ctrl_quad on a quadratic:
qf = @(w) [1 - w; 2*w.^2; w; 0.5 + w - w.^2];
w  = 0.37;
ok = chk(ok, max(abs(ctrl_quad(qf(0), qf(0.5), qf(1), 0)   - qf(0)))   < 1e-15 && ...
             max(abs(ctrl_quad(qf(0), qf(0.5), qf(1), 0.5) - qf(0.5))) < 1e-15 && ...
             max(abs(ctrl_quad(qf(0), qf(0.5), qf(1), 1)   - qf(1)))   < 1e-15 && ...
             max(abs(ctrl_quad(qf(0), qf(0.5), qf(1), w)   - qf(w)))   < 1e-15, ...
         'ctrl_quad reproduces a quadratic at node, midpoint, node and w = 0.37');

%% Synthetic flight: CR3BP + thrust, known control:
mu   = 0.012150585609624;
Tmax = 0.05;  c = 3;
dirU = [0.3; 0.9; -0.2] / sqrt(0.94);
th   = @(t) 0.3 + 0.8*t - 0.6*t.^2;                % stays in [0.3, 0.57] on [0, 1]
uOf  = @(t) [dirU; th(t)];
z0   = [0.85; 0.02; 0.01; 0.01; 0.45; 0.0; 1];
tf   = 1;  N = 10;
tN   = linspace(0, tf, N+1);
oo   = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
[~, Z] = ode113(@(t, z) cr3bp_thrust_rhs(z, uOf(t), mu, Tmax, c), tN, z0, oo);
X  = Z.';
U  = cell2mat(arrayfun(uOf, tN, 'UniformOutput', false));
tM = tN(1:end-1) + diff(tN)/2;
Um = cell2mat(arrayfun(uOf, tM, 'UniformOutput', false));
o  = struct('X', X, 'U', U, 'Um', Um, 'tNodes', tN, 'tf', tf);

%% 2. Exact nodes:
[er, ev] = flown_control_error(o, mu, Tmax, c);
ok = chk(ok, er < 1e-10 && ev < 1e-10, sprintf('exact nodes: er %.1e, ev %.1e ND', er, ev));

%% 3. Injected miss:
d  = 1e-6;
oP = o;  oP.X(1:3,end) = oP.X(1:3,end) + d*[0.6; 0; 0.8];
[erP, evP] = flown_control_error(oP, mu, Tmax, c);
oV = o;  oV.X(4:6,end) = oV.X(4:6,end) + d*[0; 1; 0];
[erV, evV] = flown_control_error(oV, mu, Tmax, c);
ok = chk(ok, abs(erP - d) < 1e-9 && evP < 1e-10, sprintf('position miss d reads back as %.4e (ev %.1e)', erP, evP));
ok = chk(ok, abs(evV - d) < 1e-9 && erV < 1e-10, sprintf('velocity miss d reads back as %.4e (er %.1e)', evV, erV));

%% 4. Reconstruction matters:
oL = o;  oL.Um = [];
[erL, evL] = flown_control_error(oL, mu, Tmax, c);
ok = chk(ok, erL > 1e3*max(er, 1e-13) && erL > 1e-7, ...
         sprintf('linear reconstruction of a quadratic throttle misses: er %.1e vs %.1e', erL, er));

%% 5. Time grid from s * tf:
oS = rmfield(o, 'tNodes');  oS.tNodes = [];  oS.s = tN / tf;
[erS, evS] = flown_control_error(oS, mu, Tmax, c);
ok = chk(ok, abs(erS - er) < 1e-12 && abs(evS - ev) < 1e-12, ...
         sprintf('tNodes = [] with s * tf gives the same errors (%.1e, %.1e)', erS, evS));

if ok, fprintf('TEST_FLOWN_CONTROL_ERROR: ALL PASS\n');
else,  fprintf('TEST_FLOWN_CONTROL_ERROR: FAILURE (see lines above)\n');
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
