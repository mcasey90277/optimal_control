function m = status_meet(a, b)
%% Purpose:
%
%   The MEET (greatest common lower bound) of two optimality-status codes
%   on the lattice ordered by what was ESTABLISHED, not by the number:
%
%           4 (sufficient)     2 (conjugate point found)
%                     \           /
%                      3 (necessary only)
%                           |
%                      1 (neither)
%
%   so 1 < 3 < 2 and 1 < 3 < 4, with 2 and 4 incomparable. Status 2 is a
%   STRONGER claim than 3: it needs everything 3 needs plus a reproducible
%   refutation. The meet of two runs' statuses is the claim both runs
%   support -- the common ground relabel_borderline keeps, and the test
%   audit_status_layer uses for a lower bound (s is reached by s' iff
%   status_meet(s, s') == s).
%
%  ASSUMPTIONS / NOTES:
%
% • meet(x, x) = x; meet(x, 1) = 1; meet(2, 3) = meet(4, 3) = 3;
%   meet(2, 4) = 3. Symmetric.
% • Both inputs must be real finite integer codes in {1, 2, 3, 4} of a
%   numeric class (char, logical, complex, NaN, empty, arrays refused).
%
%% Inputs:
%
%  a                        numeric scalar          status code 1..4
%
%  b                        numeric scalar          status code 1..4
%
%% Outputs:
%
%  m                        double scalar           the meet, a code 1..4
%
%% Revision History:
%  M. Casey                                                   (c) 10/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

assert(isCode(a) && isCode(b), 'status_meet:code', 'status codes must be integers in {1,2,3,4}');
T = [1 1 1 1; ...
     1 2 3 3; ...
     1 3 3 3; ...
     1 3 3 4];
m = T(double(a), double(b));
end

% ---------------------------------------------------------------------------
function tf = isCode(v)
% ISCODE  v is a real finite integer status code in 1..4 of a numeric class.
% INPUTS: v (any).  OUTPUTS: tf [logical].
tf = isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v == round(v) && v >= 1 && v <= 4;
end
