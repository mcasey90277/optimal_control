% TEST_MEE_XF  Structural (no-IPOPT-solve) test for casadi_lt_mee's opts.xf
% terminal-target option: confirms the default resolves to GEO [1;0;0;0;0]
% and a custom xf is honored, via the opts.selftest early-return hook (Task 1).
%
% REFERENCES: [1] .superpowers/sdd/task-1-brief.md Steps 1-2.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
par = kepler_lt_params(10,1500,2000);
% default -> GEO
o = casadi_lt_mee((0:1).', zeros(7,2), zeros(4,2), 1, ...
    struct('par',par,'x0',zeros(7,1),'selftest',true));
assert(isequal(o.xf,[1;0;0;0;0]), 'default xf must be GEO [1;0;0;0;0]');
% custom -> honored
xf = [0.9; 0.01; -0.02; 0.05; 0];
o2 = casadi_lt_mee((0:1).', zeros(7,2), zeros(4,2), 1, ...
    struct('par',par,'x0',zeros(7,1),'xf',xf,'selftest',true));
assert(isequal(o2.xf,xf), 'custom xf must be honored');
fprintf('test_mee_xf PASSED\n');
