function A = ztl_A_sun(Y, P, regime)
% ZTL_A_SUN  Exact 15x15 Jacobian A = d(dY/dtau)/dY of the Sundman-regularized
% field, by complex step of the field (never through an integrator).
%
% Within a fixed regime ztl_eom_sun is analytic in Y, so the complex-step
% derivative is exact to machine precision (Z0's derivative rule).
%
% INPUTS:
%   Y      - augmented Sundman state [15x1], REAL
%   P      - struct: .muStar .c .Tmax .eps .pSund
%   regime - 'on' | 'medium' | 'off'
% OUTPUTS:
%   A - d(dY/dtau)/dY [15x15]

h = 1e-50;
A = zeros(15);
for k = 1:15
    Yp = complex(Y(:));
    Yp(k) = Yp(k) + 1i*h;
    A(:, k) = imag(ztl_eom_sun(Yp, P, regime))/h;
end
end
