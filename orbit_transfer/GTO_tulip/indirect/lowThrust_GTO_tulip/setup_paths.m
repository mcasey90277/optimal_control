function setup_paths()
% SETUP_PATHS  Base indirect GTO->tulip campaign paths: self + shared CR3BP lib
% (cr3bp_common brings pumpkyn, which this campaign's solvers use).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
