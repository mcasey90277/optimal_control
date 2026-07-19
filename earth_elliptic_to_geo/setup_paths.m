function setup_paths()
% SETUP_PATHS  Put all earth_elliptic_to_geo code subfolders on the MATLAB path.
% Call ONCE per session before using the module (tests, the shell watchdogs, and
% interactive users all call it first). setup_paths.m and module_root.m stay at
% the module root; every other .m lives in a functional subfolder. isfolder-
% guarded so it is a safe no-op if a subfolder is absent.
r = fileparts(mfilename('fullpath'));
subs = {'core','coords','drivers','psr','verify',fullfile('verify','sosc'), ...
        'frontdoor','reproduce','viz','cartesian_legacy','lib','tests','attic'};
addpath(r);
for k = 1:numel(subs)
    d = fullfile(r, subs{k});
    if isfolder(d), addpath(d); end
end
end
