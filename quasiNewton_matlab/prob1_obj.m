function [f, g] = prob1_obj(x, A, e)
% PROB1_OBJ  Objective and gradient for Kanamori-Ohara Problem 1 (quadratic).
%
%   f(x) = 1/2 x'A x - e'x,   grad f(x) = A x - e.
% Gradient is computed only when requested (nargout > 1), so the same handle
% serves both line-search (f only) and update (f and g) calls.
%
% INPUTS:
%   x - point                          [n x 1]
%   A - SPD tridiagonal matrix         [n x n]
%   e - linear-term vector (ones)      [n x 1]
%
% OUTPUTS:
%   f - objective value                [scalar]
%   g - gradient (if requested)        [n x 1]
%
% REFERENCES:
%   Kanamori & Ohara, arXiv:1010.2846 (2010), Problem 1.

    f = 0.5 * (x' * A * x) - e' * x;
    if nargout > 1
        g = A * x - e;
    end
end
