function ok = test_drone_params()
%% Purpose:
%
%   drone_params returns the specified constants for each model, with the
%   dimensions the field will read, and refuses an unknown model by
%   identifier. The constants are the spec's verbatim (Luis & Le Ny 2016,
%   Table 1; Forster 2015); aMax is derived from them.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

k = drone_params('kinematic');
ok = chk(ok, strcmp(k.model, 'kinematic') && k.nx == 3 && k.nu == 3, 'kinematic: nx = 3, nu = 3');
ok = chk(ok, isequal(sort(fieldnames(k)), sort({'model'; 'nx'; 'nu'})), 'kinematic: exactly model, nx, nu');

q = drone_params('pointmass');
ok = chk(ok, q.nx == 6 && q.nu == 3 && q.g == 9.81, 'pointmass: nx = 6, nu = 3, g = 9.81');
ok = chk(ok, abs(q.aMax - (4*0.15/0.027 - 9.81)) < 1e-12, 'pointmass: aMax = 4 fMax / m - g');

r = drone_params('rigid6dof');
want = struct('m', 0.027, 'l', 0.046, 'kappa', 0.005964, 'g', 9.81, 'fMax', 0.15, 'thetaMax', pi/3);
f = fieldnames(want);
for n = 1:numel(f)
    ok = chk(ok, isfield(r, f{n}) && r.(f{n}) == want.(f{n}), sprintf('rigid6dof: %s = %g', f{n}, want.(f{n})));
end
ok = chk(ok, isequal(size(r.I), [3 1]) && isequal(r.I, [1.395e-5; 1.436e-5; 2.173e-5]), 'rigid6dof: I is the 3x1 Crazyflie diagonal');
ok = chk(ok, r.nx == 12 && r.nu == 4, 'rigid6dof: nx = 12, nu = 4');
ok = chk(ok, strcmp(r.model, 'rigid6dof') && strcmp(drone_params('RIGID6DOF').model, 'rigid6dof'), 'model name is case-insensitive on input, lower-case on output');

try
    drone_params('hexacopter');
    ok = chk(ok, false, 'unknown model refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_params:model'), 'unknown model refused by identifier drone_params:model');
end

if ok, fprintf('TEST_DRONE_PARAMS: ALL PASS\n'); else, fprintf('TEST_DRONE_PARAMS: FAILURE\n'); end
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
