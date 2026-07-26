root = fileparts(fileparts(mfilename('fullpath'))); cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
T = recertify_table3([10 5 2.5 1 0.5]);
save('/tmp/task10_T.mat','T');
disp('=== struct dump ===');
for k=1:numel(T)
    fprintf('T(%d): thrustN=%g tag=%s verdict=%s status=%s nFlat=%g method=%s robust=%d red=[%g %g %g] drift=%.3e stat=%.3e\n', ...
        k, T(k).thrustN, T(k).tag, T(k).verdict, T(k).status, T(k).nFlat, T(k).method, T(k).robust, ...
        T(k).red.npos, T(k).red.nneg, T(k).red.nzero, T(k).drift, T(k).stat);
end
