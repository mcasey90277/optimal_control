function r = module_root()
% MODULE_ROOT  Absolute path to the earth_elliptic_to_geo module root, so code
% can locate the shared results/ cache dir regardless of which subfolder the
% calling file lives in (this file stays at the module top level through the
% Stage-B subfolder reorg; fileparts of its own path is the root).
% OUTPUTS: r - module root directory [char]
r = fileparts(mfilename('fullpath'));
end
