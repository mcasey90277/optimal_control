function A = ztl_A(y, P, regime)
% ZTL_A  Exact 14x14 Jacobian A = df/dy of ztl_eom's field at fixed regime,
% by complex step OF THE FIELD (never through an integrator).
%
% Within a fixed regime the ztl_eom field is analytic in y (no abs/min/max/
% norm/real -- see its coding rules), so the complex-step derivative is
% exact to machine precision. This is the load-bearing derivative rule of
% the ZTL build (Z0_BUILD.md par.3): differencing the FLOW is forbidden;
% probing the FIELD is exact.
%
% INPUTS:
%   y      - augmented state [14x1], REAL
%   P      - struct: .muStar .c .Tmax .eps
%   regime - 'on' | 'medium' | 'off'
%
% OUTPUTS:
%   A - df/dy [14x14]
%
% REFERENCES:
%   [1] Squire & Trapp, SIAM Review 40(1), 1998 (complex-step derivative).

h = 1e-50;
A = zeros(14);
for k = 1:14
    yp = complex(y(:));
    yp(k) = yp(k) + 1i*h;
    A(:, k) = imag(ztl_eom(yp, P, regime))/h;
end
end
