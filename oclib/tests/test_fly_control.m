function ok = test_fly_control()
%% Purpose:
%
%   Unit test for oc.fly_control (the G1b / G2 flown-control engine) on
%   problems with exact answers:
%     1. EXACT FLOW, BOTH MODES: a damped oscillator dz = A z flown
%        'perInterval' and 'span' lands on expm(A t_f) z0 to < 1e-10.
%     2. TRAJECTORY: out.t runs t_0 -> t_f, out.Z matches it row for row
%        and ends at zEnd; in 'perInterval' every node time is present.
%     3. ONE OUTPUT: zEnd is the same whether or not out is requested.
%     4. CONTROL FROM GLOBAL TIME: a piecewise-linear control looked up
%        from t (the closure contract in the header) reproduces its exact
%        integral, a piecewise quadratic, in both modes.
%     5. SOLVER OPTION: @ode45 at its own tolerances is honoured and still
%        lands within them.
%     6. REFUSAL: an unknown mode throws oc:fly_control:mode.
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

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                         % oclib root -> +oc visible
ok = true;

%% 1. Exact flow, both modes:
A    = [0 1; -2 -0.3];                            % damped oscillator
rhs  = @(t, z) A*z;
z0   = [1; 0];
tG   = linspace(0, 3, 9);
zRef = expm(A*tG(end)) * z0;
[zP, outP] = oc.fly_control(z0, tG, rhs);
[zS, outS] = oc.fly_control(z0, tG, rhs, struct('mode', 'span'));
ok = check(sprintf('perInterval lands on the exact flow (%.1e)', max(abs(zP - zRef))), ...
           max(abs(zP - zRef)) < 1e-10) && ok;
ok = check(sprintf('span lands on the exact flow (%.1e)', max(abs(zS - zRef))), ...
           max(abs(zS - zRef)) < 1e-10) && ok;

%% 2. Trajectory output:
ok = check('perInterval trajectory spans t0..tf and ends at zEnd', ...
           outP.t(1) == tG(1) && outP.t(end) == tG(end) ...
           && size(outP.Z, 1) == numel(outP.t) && isequal(outP.Z(end,:).', zP)) && ok;
ok = check('perInterval trajectory holds every node time', ...
           all(ismember(tG, outP.t))) && ok;
ok = check('span trajectory spans t0..tf and ends at zEnd', ...
           outS.t(1) == tG(1) && outS.t(end) == tG(end) && isequal(outS.Z(end,:).', zS)) && ok;

%% 3. One output gives the same terminal state:
ok = check('zEnd identical with and without the trajectory output', ...
           isequal(oc.fly_control(z0, tG, rhs), zP)) && ok;

%% 4. Control reconstructed from global time:
tC   = [0 0.7 1.5 2.6 3];                         % non-uniform nodes
uC   = [1 -2 0.5 3 -1];                           % nodal control values
rhsU = @(t, z) interp1(tC, uC, t, 'linear');      % dz = u(t)
zEx  = trapz(tC, uC);                             % exact: piecewise-linear integral
zPu  = oc.fly_control(0, tC, rhsU);
zSu  = oc.fly_control(0, tC, rhsU, struct('mode', 'span'));
ok = check(sprintf('interval-lookup control, perInterval (%.1e)', abs(zPu - zEx)), ...
           abs(zPu - zEx) < 1e-10) && ok;
ok = check(sprintf('interval-lookup control, span (%.1e)', abs(zSu - zEx)), ...
           abs(zSu - zEx) < 1e-6) && ok;

%% 5. Solver option:
o45 = struct('solver', @ode45, 'RelTol', 1e-9, 'AbsTol', 1e-11);
z45 = oc.fly_control(z0, tG, rhs, o45);
ok = check(sprintf('@ode45 honoured, within its tolerance (%.1e)', max(abs(z45 - zRef))), ...
           max(abs(z45 - zRef)) < 1e-7 && ~isequal(z45, zP)) && ok;

%% 6. Refusal:
ok = check('unknown mode refused', ...
           throwsId(@() oc.fly_control(z0, tG, rhs, struct('mode', 'rk4')), ...
                    'oc:fly_control:mode')) && ok;

fprintf('TEST_FLY_CONTROL: %s\n', passText(ok));
end

% ------------------------------------------------------------------------
function ok = check(name, cond)
%% Purpose:
%
%   Print one gate line and return its verdict.
%
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end

function tf = throwsId(fn, id)
%% Purpose:
%
%   True when fn() throws an error whose identifier is id.
%
tf = false;
try
    fn();
catch ME
    tf = strcmp(ME.identifier, id);
end
end

function s = passText(ok)
%% Purpose:
%
%   'ALL PASS' or 'FAILED'.
%
if ok, s = 'ALL PASS'; else, s = 'FAILED'; end
end
