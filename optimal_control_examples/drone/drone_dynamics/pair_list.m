function pairs = pair_list(M)
%% Purpose:
%
%   Every unordered pair of M agents, as a 2 x P index array with the
%   smaller index first, lexicographic order. P = M (M - 1) / 2.
%
%% Inputs:
%
%  M                        double                  Number of agents
%
%% Outputs:
%
%  pairs                    [2 x P]                 [a; b] with a < b
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo:
     disp(pair_list(4));
     return
end

P = M*(M - 1)/2;
pairs = zeros(2, P);
n = 0;
for a = 1:M-1
    for b = a+1:M
        n = n + 1;
        pairs(:, n) = [a; b];
    end
end
end
