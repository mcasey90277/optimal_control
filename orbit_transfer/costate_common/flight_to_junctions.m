function [Y, tGrid] = flight_to_junctions(tj, yj, K, opts)
%% Purpose:
%
%   A flown trajectory cut into the K+1 JUNCTION STATES a multiple-shooting
%   seed needs: interpolate the flight onto a uniform grid in normalized
%   time, and -- when the seed is for a DIFFERENT thrust or time of flight
%   than the flight it came from -- replace the mass row by the all-burn
%   identity rather than rescaling the old one.
%
%   Six engines wrote these two lines inline (extend_thrust_ladder,
%   densify_ladder, lowthrust_ladder, probe_deep_rungs, probe_abstract_case
%   and the direct warm start), and `seed_from_z8` held a seventh copy.
%
%  ASSUMPTIONS / NOTES:
%
% • The query is in NORMALIZED time (`tu/tu(end)`), the form six of the
%   seven sites used. `seed_from_z8` queried in ABSOLUTE time, which differs
%   in the last bits: measured 3e-13 absolute / 5.5e-14 relative on a real
%   70 mN flight, all of it in the costate rows. That is a seed
%   perturbation, not a solution change -- the root is set by the
%   boundary-value problem, not by the guess -- and `golden_cells` gates it
%   (iteration counts and residuals included).
%
% • Duplicate propagator time samples are dropped before interpolation:
%   `interp1` refuses a repeated site, and a variable-step integrator can
%   emit one.
%
% • The mass law is DERIVED, never rescaled from the source flight: a
%   rescale once carried a scaling mistake through a rung change (review
%   finding recorded in extend_thrust_ladder).
%
%% Inputs:
%
%  tj                       [m x 1]                 Flight times, ascending
%  yj                       [m x n]                 Flight states, one ROW
%                                                   per time (n >= 7 for a
%                                                   mass law: row 7 is mass)
%  K                        double                  Segment count; K+1
%                                                   junctions
%  opts                     struct (optional)
%   .tf [tj(end)] the SEED's time of flight, which may differ from the
%   flight's own (a rung-scaled guess)
%   .massLaw [] struct('Tnd', , 'cnd', ): replace row 7 with the all-burn
%   mass m = 1 - Tnd*t/cnd on the seed's own grid
%   .rows [1:min(14, n)] which columns of yj to carry
%
%% Outputs:
%
%  Y                        [numel(rows) x K+1]     Junction states
%  tGrid                    [1 x K+1]               linspace(0, tf, K+1)
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
             tD = linspace(0, 2, 201).';
             yD = [sin(tD), cos(tD), tD, tD.^2, sin(2*tD), cos(2*tD), 1 - 0.05*tD];
    [YD, tgD] = flight_to_junctions(tD, yD, 8, struct('rows', 1:7));
    fprintf('8 segments -> Y is %d x %d, grid 0 .. %.3f\n', size(YD, 1), size(YD, 2), tgD(end));
    [Y, tGrid] = deal(YD, tgD);
    return;
end

if nargin < 4, opts = struct(); end
             tj = tj(:);

if size(yj, 1) ~= numel(tj)
    error('flight_to_junctions:size', ...
          'yj has %d rows but tj has %d times: yj must carry one ROW per time', ...
          size(yj, 1), numel(tj));
end
if ~(isscalar(K) && isfinite(K) && K >= 1 && K == round(K))
    error('flight_to_junctions:segments', 'K must be a positive whole number of segments, got %g', K);
end
           rows = fieldd(opts, 'rows', 1:min(14, size(yj, 2)));
             tf = fieldd(opts, 'tf', tj(end));
        massLaw = fieldd(opts, 'massLaw', []);

%% Uniform grid in normalized time; duplicate sites are dropped because
%% interp1 refuses them and an adaptive integrator can emit one:
          sGrid = linspace(0, 1, K+1);
       [tu, iu] = unique(tj);
              Y = interp1(tu/tu(end), yj(iu, rows), sGrid, 'pchip')';
          tGrid = sGrid*tf;

%% The all-burn mass row, DERIVED on the seed's own grid:
if ~isempty(massLaw)
    if numel(rows) < 7
        error('flight_to_junctions:massRow', 'a mass law needs at least 7 rows, got %d', numel(rows));
    end
         Y(7,:) = 1 - massLaw.Tnd*tGrid/massLaw.cnd;
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
