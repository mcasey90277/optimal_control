function [L, roots] = deliverable_file_list(profile)
%% Purpose:
%
%   The single list of files a 70 mN costate-library deliverable ships
%   (the catalog itself excluded: the builder exports it with
%   export_deliverable_catalog). Consumed by the rehearsal test and by the
%   deliverable builder, so each list is defined exactly once. Repo-only:
%   this file is NOT shipped.
%
%   Two profiles:
%
%   * 'developer' (deliverable 8): the status-layer tools WITH the flight
%     layer and the tfMin acceptance check (25 files).
%   * 'user' (deliverable 9): what a user of the library needs and nothing
%     that builds, flies or re-checks entries -- describe, pick, list the
%     transfers at a cell, the status legend, extremes, the status torus,
%     the movie, the one orbit helper the extremes and the movie draw with,
%     and one worked example (9 files).
%
%  ASSUMPTIONS / NOTES:
%
% • Source folders are computed from this file's own location, never typed.
% • Files under pumpkyn, pumpkynPie, movie2gifLite or matlabroot are not
%   listed: the recipient has them.
%
%% Inputs:
%
%  profile                  char (optional)         'developer' [default] |
%                                                   'user'
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
%  M. Casey  'user' profile (deliverable 9)                    (c) 10/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1 || isempty(profile), profile = 'developer'; end
here = fileparts(mfilename('fullpath'));
roots = struct('indirect', here, ...
               'costate_common', fullfile(fileparts(fileparts(here)), 'costate_common'), ...
               'examples', fullfile(here, 'examples'));

switch profile
    case 'developer'
        cc  = {'get_family_orbit', 'phase_state', 'periodic_pp', 'validate_flight', ...
               'scalar_verdict', 'catalog_schema', 'fly_transfer', ...
               'run_capped', 'current_pool'};
        ind = {'costate_catalog_pick', 'transfers_at', 'status_key', 'catalog_endpoints', ...
               'fly_catalog_entry', 'tfmin_accepts', 'entry_thrust_program', 'costate_lib_describe', ...
               'costate_catalog_extremes', 'plot_status_torus', 'costate_catalog_extremes_movies', ...
               'print_transfer_summary', 'plot_transfer_3d', 'sample_status_pair'};
        ex  = {'costate_catalog_example', 'status_layer_example'};
    case 'user'
        cc  = {'get_family_orbit'};
        ind = {'costate_lib_describe', 'costate_catalog_pick', 'transfers_at', 'status_key', ...
               'costate_catalog_extremes', 'plot_status_torus', 'costate_catalog_extremes_movies'};
        ex  = {'costate_library_example'};
    otherwise
        error('deliverable_file_list:profile', 'unknown profile %s (developer | user)', profile);
end

files = [strcat(cc, '.m'), strcat(ind, '.m'), strcat(ex, '.m')];
from  = [repmat({'costate_common'}, 1, numel(cc)), repmat({'indirect'}, 1, numel(ind)), ...
         repmat({'examples'}, 1, numel(ex))];
L = struct('file', files(:), 'from', from(:));
end
