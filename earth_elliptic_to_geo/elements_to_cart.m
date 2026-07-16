function [r, v] = elements_to_cart(P, ex, ey, hx, hy, L, mu)
% ELEMENTS_TO_CART  Paper/MEE-style elements -> inertial Cartesian state.
%
% Elements per Haberkorn-Martinon-Gergaud 2004: P (semi-latus rectum),
% (ex,ey) = e*(cos,sin)(Om+om), (hx,hy) = tan(i/2)*(cos,sin)(Om), L = Om+om+theta.
%
% INPUTS:  P,ex,ey,hx,hy,L - elements [scalars];  mu - grav parameter [scalar]
% OUTPUTS: r, v - inertial position/velocity [3x1 each]
%
% REFERENCES:
%   [1] Walker/Betts modified-equinoctial <-> Cartesian formulas.
w  = 1 + ex*cos(L) + ey*sin(L);
s2 = 1 + hx^2 + hy^2;
a2 = hx^2 - hy^2;
rm = P / w;
r  = (rm/s2) * [cos(L) + a2*cos(L) + 2*hx*hy*sin(L);
                sin(L) - a2*sin(L) + 2*hx*hy*cos(L);
                2*(hx*sin(L) - hy*cos(L))];
sq = sqrt(mu/P);
v  = (1/s2) * [-sq*( sin(L) + a2*sin(L) - 2*hx*hy*cos(L) + ey - 2*ex*hx*hy + a2*ey);
               -sq*(-cos(L) + a2*cos(L) + 2*hx*hy*sin(L) - ex + 2*ey*hx*hy + a2*ex);
                2*sq*(hx*cos(L) + hy*sin(L) + ex*hx + ey*hy)];
end
