function setup_paths()
% SETUP_PATHS  min_time (PMP min-time root) paths: self + shared CR3BP lib
% (gto_tulip_endpoints + pumpkyn via cr3bp_common).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
