function s = fmt_num(v, w, p, unit)
%% Purpose:
%
%   Format a number at a FIXED WIDTH, rendering the missing case as dashes
%   of that same width instead of throwing.
%
%   The two shapes that broke a campaign print: a NaN where a certified
%   transfer time was expected (char(string(NaN)) is <missing>, and char()
%   refuses it), and an empty value where a count was expected. Both are
%   ordinary states of a sparse grid, not errors.
%
%   Fixed width also keeps a movie title or a table column from moving
%   frame to frame, which is the other reason to have one of these.
%
%% Inputs:
%
%  v                        any                     the value
%  w                        double                  total field width
%  p                        double                  decimals (default 4)
%  unit                     char (optional)         appended verbatim
%
%% Outputs:
%
%  s                        char                    exactly w chars, plus
%                                                   any unit. A value too
%                                                   wide for w at precision
%                                                   p is shown in %g at
%                                                   width w, and as '#'s if
%                                                   even that does not fit.
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(p), p = 4; end
if nargin < 4, unit = ''; end
if isempty(v) || ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
    s = repmat('-', 1, w);
else
    s = sprintf('%*.*f', w, p, v);
    if numel(s) > w                       % too wide for the column: say so
        s = sprintf('%*.*g', w, max(1, w - 6), v);   % at width, fewer digits
        if numel(s) > w, s = repmat('#', 1, w); end  % cannot be shown at all
    end
end
if ~isempty(unit), s = [s ' ' unit]; end
end
