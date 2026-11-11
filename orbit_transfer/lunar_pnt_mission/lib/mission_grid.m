function g = mission_grid(name)
%% Purpose:
%
%   The lunar south-pole surface grids of the PNT mission study, as
%   latitude chunks for proj7's gdop_grid_from_states (which evaluates
%   lat_range(1):lat_res:lat_range(2) x lon_range(1):lon_res:lon_range(2)
%   per call; the chunks bound its [Nt x 3 x M x Npts] working array).
%
%   'paper'  : the run_requirements_table grid behind the paper's Table 1:
%              five chunks with lat_range pairs [-89.9 -87] [-86.5 -84]
%              [-83.5 -81] [-80.5 -78] [-77.5 -75], lat_res 0.5, lon 0:1:358
%              (30 latitudes x 359 longitudes = 10,770 points; chunk 1
%              evaluates -89.9:0.5:-87.4, NOT -89.9, -89, -88.5, ...),
%              cadence 348.1 s.
%   'coarse' : a pole chunk [-89.9 -89.9] plus [-89.5 -75.5] at lat_res 2,
%              lon 0:2:358 (9 x 180 = 1,620 points), cadence 300 s.
%   'pole'   : ONE point, the lunar south pole: latitude -90 exactly,
%              longitude 0, cadence 300 s (pumpkyn fromLLA puts it at
%              1737.1 km straight below the Moon's centre; checked in
%              test_validate_pole_series).
%
%% Inputs:
%
%  name                     char                    'paper' | 'coarse' | 'pole'
%
%% Outputs:
%
%  g.name                   char                    The grid name
%  g.chunks                 {1 x C} cell of [1 x 2] lat_range pairs [deg]
%  g.lat_res                scalar                  Latitude step [deg]
%  g.lon_range              [1 x 2]                 Longitude range [deg]
%  g.lon_res                scalar                  Longitude step [deg]
%  g.dtSec                  scalar                  Epoch cadence [s]
%  g.mask                   scalar                  Elevation mask [deg]
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  'pole' grid (one point at -90, 0)                (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
    for nm = {'paper', 'coarse', 'pole'}
        gg = mission_grid(nm{1});
        nP = 0;
        for cc = 1:numel(gg.chunks)
            nP = nP + numel(gg.chunks{cc}(1):gg.lat_res:gg.chunks{cc}(2)) * ...
                      numel(gg.lon_range(1):gg.lon_res:gg.lon_range(2));
        end
        fprintf('%-6s: %d chunks, %d points, %.1f s cadence\n', gg.name, ...
                numel(gg.chunks), nP, gg.dtSec);
    end
    return;
end

switch lower(char(name))
    case 'paper'
        g.name      = 'paper';
        g.chunks    = {[-89.9 -87], [-86.5 -84], [-83.5 -81], [-80.5 -78], [-77.5 -75]};
        g.lat_res   = 0.5;
        g.lon_range = [0 358];
        g.lon_res   = 1;
        g.dtSec     = 348.1;
        g.mask      = 0;
    case 'coarse'
        g.name      = 'coarse';
        g.chunks    = {[-89.9 -89.9], [-89.5 -75.5]};
        g.lat_res   = 2;
        g.lon_range = [0 358];
        g.lon_res   = 2;
        g.dtSec     = 300;
        g.mask      = 0;
    case 'pole'
        g.name      = 'pole';
        g.chunks    = {[-90 -90]};
        g.lat_res   = 1;
        g.lon_range = [0 0];
        g.lon_res   = 1;
        g.dtSec     = 300;
        g.mask      = 0;
    otherwise
        error('mission_grid:name', 'unknown grid "%s" (paper | coarse | pole)', char(name));
end
end
