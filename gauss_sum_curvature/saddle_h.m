function [z, J, Hcell] = saddle_h(x)
% SADDLE_H  Saddle embedding h: R^2 -> R^3 with indefinite measurement curvature.
%
% Models a nonlinear measurement z = h(x) whose embedded surface is a hyperbolic
% paraboloid (a saddle), so the curvature of the third component is INDEFINITE
% (curves up along x1, down along x2). Used to test curvature estimation and
% curvature-driven Gaussian-sum splitting.
%
%   h(x) = [ x1 ; x2 ; (x1^2 - x2^2)/2 ]
%
% INPUTS:
%   x - state                                              [2 x 1]
%
% OUTPUTS:
%   z     - measurement h(x)                               [3 x 1]
%   J     - Jacobian dh/dx                                 [3 x 2]
%   Hcell - 1x3 cell of component Hessians d^2 h_k/dx^2,
%           each [2 x 2]; H{1}=H{2}=0, H{3}=diag(1,-1)     {1 x 3}
%
% NOTE:
%   The "Hessian of h" is a 3rd-order tensor -- one 2x2 matrix per output
%   component. Only the third component is curved, and it is indefinite, which
%   is the point: measurement-model curvature is generally not positive definite.

    x = x(:);
    z = [x(1); x(2); 0.5*(x(1)^2 - x(2)^2)];
    J = [1, 0;
         0, 1;
         x(1), -x(2)];
    if nargout > 2
        Z = zeros(2,2);
        Hcell = {Z, Z, [1 0; 0 -1]};
    end
end
