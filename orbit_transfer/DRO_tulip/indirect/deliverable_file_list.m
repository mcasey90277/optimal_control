function [L, roots] = deliverable_file_list()
%% Purpose:
%
%   The single list of files the 70 mN costate-library deliverable ships
%   (the catalog itself excluded: the builder exports it with
%   export_deliverable_catalog). Consumed by the rehearsal test and by the
%   deliverable builder, so the list is defined exactly once. Repo-only:
%   this file is NOT shipped.
%
%  ASSUMPTIONS / NOTES:
%
% • Source folders are computed from this file's own location, never typed.
% • Files under pumpkyn, pumpkynPie, movie2gifLite or matlabroot are not
%   listed: the recipient has them.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  L                        struct array [N x 1]    .file (bare file name)
%                                                   .from ('indirect' |
%                                                   'costate_common' |
%                                                   'examples')
%
%  roots                    struct                  .indirect
%                                                   .costate_common
%                                                   .examples: the
%                                                   ABSOLUTE source folders
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
roots = struct('indirect', here, ...
               'costate_common', fullfile(fileparts(fileparts(here)), 'costate_common'), ...
               'examples', fullfile(here, 'examples'));

cc  = {'get_family_orbit', 'phase_state', 'periodic_pp', 'validate_flight', ...
       'scalar_verdict', 'catalog_schema', 'fly_transfer', ...
       'run_capped', 'current_pool'};
ind = {'costate_catalog_pick', 'transfers_at', 'status_key', 'catalog_endpoints', ...
       'fly_catalog_entry', 'tfmin_accepts', 'entry_thrust_program', 'costate_lib_describe', ...
       'costate_catalog_extremes', 'plot_status_torus', 'costate_catalog_extremes_movies', ...
       'print_transfer_summary', 'plot_transfer_3d', 'sample_status_pair'};
ex  = {'costate_catalog_example', 'status_layer_example'};

files = [strcat(cc, '.m'), strcat(ind, '.m'), strcat(ex, '.m')];
from  = [repmat({'costate_common'}, 1, numel(cc)), repmat({'indirect'}, 1, numel(ind)), ...
         repmat({'examples'}, 1, numel(ex))];
L = struct('file', files(:), 'from', from(:));
end
