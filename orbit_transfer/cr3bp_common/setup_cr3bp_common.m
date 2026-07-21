function setup_cr3bp_common()
% SETUP_CR3BP_COMMON  Add the shared CR3BP GTO-transfer library + pumpkyn.
%
% Single source of the cross-module problem definition (cr3bp_lt_params,
% minfuel_config, gto_tulip_endpoints, gto_elfo_endpoints) and the pumpkyn
% toolbox path. Called by every GTO module's setup_paths.m.
%
% INPUTS:  (none)
% OUTPUTS: (none) - modifies the MATLAB path in-place
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-21-gto-direct-indirect-restructure-design.md
here = fileparts(mfilename('fullpath'));
addpath(here);
pumpkynSrc = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', ...
                      'pumpkyn', 'src');
assert(exist(fullfile(pumpkynSrc, '+pumpkyn'), 'dir') == 7, ...
    'setup_cr3bp_common:missing', 'pumpkyn not found at %s', pumpkynSrc);
addpath(pumpkynSrc);
end
