function [X, W, tf] = unpack_z(Z, N)
% UNPACK_Z  Split the NLP decision vector into states, controls, and tf.
%
% INPUTS:
%   Z  - decision vector [10*(N+1)+1 x 1]: [X(:); W(:); tf]
%   N  - number of trapezoidal segments [scalar]
%
% OUTPUTS:
%   X  - states at the N+1 nodes [7x(N+1)]
%   W  - unit thrust-direction controls at the N+1 nodes [3x(N+1)]
%   tf - transfer time (ND) [scalar]

nNodes = N + 1;
X  = reshape(Z(1:7*nNodes), 7, nNodes);
W  = reshape(Z(7*nNodes + (1:3*nNodes)), 3, nNodes);
tf = Z(end);
end
