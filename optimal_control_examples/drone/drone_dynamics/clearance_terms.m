function [h, rho, n, hBall] = clearance_terms(X, pairs, geom)
%% Purpose:
%
%   Pairwise clearance of a set of agent configurations: for every pair
%   (a, b) and every column of X, the gap h = |x_a - x_b| - 2 r, the
%   distance rho, and the unit vector n = (x_a - x_b)/rho. With geom.ball
%   set, also the obstacle gap |x_i - c| - R per agent. Complex-step safe:
%   distances are sqrt of sums of squares, never a magnitude.
%
%  ASSUMPTIONS / NOTES:
%
% • rho = 0 (coincident agents) gives n = NaN; the barrier and the solver
%   never evaluate there (the NLP uses the squared form, the barrier's
%   quadratic branch is only reached at h <= l > -2r).
%
%% Inputs:
%
%  X                        [3M x N]                Agent-major stacked
%                                                   positions, one
%                                                   configuration per
%                                                   column
%
%  pairs                    [2 x P]                 From pair_list
%
%  geom                     struct                  .r, .ball ([] or .c .R)
%
%% Outputs:
%
%  h                        [P x N]                 rho - 2r
%
%  rho                      [P x N]                 Pair distances
%
%  n                        [3 x P x N]             Unit vectors a -> b
%
%  hBall                    [M x N] or []           |x_i - c| - R
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 circle's initial clearances:
     [a, ~, g] = drone_layouts('paper', 4);
     hh = clearance_terms(a(:), pair_list(4), g);
     fprintf('h at t = 0: [%s]\n', sprintf('%.3f ', hh));
     return
end

M = size(X, 1)/3;
N = size(X, 2);
P = size(pairs, 2);
h   = zeros(P, N, 'like', X);
rho = zeros(P, N, 'like', X);
n   = zeros(3, P, N, 'like', X);
for c = 1:N
    Pc = reshape(X(:, c), 3, M);
    for q = 1:P
        dv = Pc(:, pairs(1, q)) - Pc(:, pairs(2, q));
        dd = sqrt(sum(dv.^2));
        rho(q, c) = dd;
        h(q, c) = dd - 2*geom.r;
        n(:, q, c) = dv/dd;
    end
end
hBall = [];
if isfield(geom, 'ball') && ~isempty(geom.ball)
    hBall = zeros(M, N, 'like', X);
    for c = 1:N
        Pc = reshape(X(:, c), 3, M);
        for a = 1:M
            dv = Pc(:, a) - geom.ball.c;
            hBall(a, c) = sqrt(sum(dv.^2)) - geom.ball.R;
        end
    end
end
end
