root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
t = sosc_defaults();
for f = {'recon','drift','stat','feas','dual','comp','active','mu','inertiaZero'}
    assert(isfield(t,f{1}) && t.(f{1})>0, sprintf('missing/nonpos tol.%s',f{1}));
end
assert(t.feas==1e-8, 'feas must match the existing maxDefect<1e-8 gate');
fprintf('test_sosc_defaults PASSED\n');
