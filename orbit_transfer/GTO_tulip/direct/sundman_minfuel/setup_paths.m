function setup_paths()
% SETUP_PATHS  Direct Sundman min-fuel engine paths: self + shared CR3BP lib
% (cr3bp_common: params/config/endpoints + pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
