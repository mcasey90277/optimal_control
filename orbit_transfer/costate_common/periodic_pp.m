function [pp, seam, dpp] = periodic_pp(t, y, opts)
%% Purpose:
%
%   The C1-PERIODIC cubic interpolant through ONE PERIOD of a closed orbit,
%   with its derivative and its seam mismatch.
%
%   Why it is not `spline`: an ordinary not-a-knot cubic is not C1 across
%   the seam at s = 0 -- exactly where a phase near zero is evaluated, and
%   exactly where a continuation needs d x / d(phase). On the 7-petal tulip
%   the not-a-knot derivative jumps by ~1e-2 there while the periodic form
%   closes to ~1e-15. In VALUE the two agree to well under a millimetre
%   except within a fraction of a percent of the seam, where the difference
%   reaches metres (measured 2026-09-11: 3.2 km/1000 on the DRO at s =
%   0.999) -- so this matters most for derivatives, and for anything that
%   evaluates near the seam.
%
%   Extracted 2026-09-11 from two private copies that disagreed on policy:
%   transfer_study's periodicPP (errored without the Curve Fitting Toolbox)
%   and arclength_arrival's makePP (fell back to an ordinary spline in
%   silence). The policy is now an OPTION, and the fallback says so.
%
%  ASSUMPTIONS / NOTES:
%
% • Samples run DOWN THE ROWS: y is [m x n] for m sample times.
% • `seam.value` is the DATA's own closure |y(end,:) - y(1,:)| -- the orbit's
%   periodicity error, which no interpolant can improve. `seam.deriv` is the
%   INTERPOLANT's derivative mismatch across the seam, which is what the
%   scheme controls.
% • `dpp` differentiates the interpolant's OWN coefficients, so it is exactly
%   the derivative of what ppval evaluates (not a re-fit).
%
%% Inputs:
%
%  t                        [m x 1] or [1 x m]      sample times of one
%                                                   period, strictly
%                                                   increasing, m >= 4
%
%  y                        [m x n]                 samples, one row per time
%
%  opts                     struct (optional)
%   .scheme ['periodic'] or 'notaknot' (the ordinary spline -- for
%   comparison and for callers that do not touch the seam)
%   .onMissing ['error'] or 'notaknot': what to do when scheme is
%   'periodic' and the Curve Fitting Toolbox (csape) is absent
%
%% Outputs:
%
%  pp                       struct                  piecewise polynomial,
%                                                   evaluate with ppval
%
%  seam                     struct                  .value (the data's
%                                                   closure) .deriv (the
%                                                   interpolant's derivative
%                                                   mismatch at the seam)
%                                                   .scheme (what was built)
%
%  dpp                      struct                  derivative of pp
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
            thD = linspace(0, 2*pi, 61).';
     [ppD, sD] = periodic_pp(thD, [cos(thD), sin(thD)]);
    [~, sNK]   = periodic_pp(thD, [cos(thD), sin(thD)], struct('scheme', 'notaknot'));
    fprintf('circle: periodic seam derivative %.2e, not-a-knot %.2e\n', sD.deriv, sNK.deriv);
    [pp, seam, dpp] = deal(ppD, sD, []);
    return;
end

if nargin < 3, opts = struct(); end
           scheme = lower(fieldd(opts, 'scheme', 'periodic'));
        onMissing = lower(fieldd(opts, 'onMissing', 'error'));
                t = t(:).';

%% Refuse malformed input by name -- a silently transposed table would
%% interpolate the wrong thing:
if size(y, 1) ~= numel(t)
    error('periodic_pp:size', ...
          'y has %d rows but t has %d samples: y must be [m x n] with one ROW per time', ...
          size(y, 1), numel(t));
end
if numel(t) < 4
    error('periodic_pp:tooFew', 'a cubic through one period needs at least 4 samples, got %d', numel(t));
end
if any(diff(t) <= 0)
    error('periodic_pp:tIncreasing', 'the sample times must be strictly increasing');
end

%% Build the interpolant:
                Y = y.';                        % n x m, as csape/spline want
switch scheme
    case 'periodic'
        if exist('csape', 'file') == 2
               pp = csape(t, Y, 'periodic');
             used = 'periodic';
        else
            switch onMissing
                case 'notaknot'
                       pp = spline(t, Y);
                     used = 'notaknot (Curve Fitting Toolbox absent)';
                    warning('periodic_pp:noToolbox', ...
                            ['csape is unavailable: falling back to an ordinary spline, whose ' ...
                             'derivative jumps at the seam']);
                otherwise
                    error('periodic_pp:noToolbox', ...
                          ['csape (Curve Fitting Toolbox) is required for a periodic interpolant. ' ...
                           'Pass opts.onMissing = ''notaknot'' to accept an ordinary spline, whose ' ...
                           'derivative jumps at the seam']);
            end
        end
    case 'notaknot'
               pp = spline(t, Y);
             used = 'notaknot';
    otherwise
        error('periodic_pp:scheme', 'unknown scheme "%s": use ''periodic'' or ''notaknot''', scheme);
end

%% Its derivative, from its own coefficients:
              dpp = ppder(pp);

%% The two seams, which measure different things (see the notes above):
       seam.value = norm(y(end,:) - y(1,:));
       seam.deriv = norm(ppval(dpp, t(1)) - ppval(dpp, t(end)));
      seam.scheme = used;
end

% ------------------------------------------------------------------------
function dpp = ppder(pp)
% PPDER  Derivative of a piecewise polynomial, by differentiating its own
% coefficients -- so it is exactly the derivative of what ppval evaluates.
% INPUTS: pp.  OUTPUTS: dpp.
[br, co, ~, or, dm] = unmkpp(pp);
if or == 1
    dpp = mkpp(br, zeros(size(co, 1), 1), dm);  return
end
              w = (or-1):-1:1;                  % powers of the derivative
            dpp = mkpp(br, co(:, 1:or-1) .* w, dm);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
