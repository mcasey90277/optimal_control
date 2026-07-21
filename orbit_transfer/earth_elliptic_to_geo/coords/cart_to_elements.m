function el = cart_to_elements(r, v, mu)
% CART_TO_ELEMENTS  Inertial Cartesian state -> paper/MEE-style elements.
%
% INPUTS:  r, v - inertial position/velocity [3x1];  mu - grav parameter
% OUTPUTS: el - struct .P .ex .ey .hx .hy .L   (L in (-pi, pi])
%
% REFERENCES: inverse of elements_to_cart (roundtrip-tested).
r = r(:);  v = v(:);
hv = cross(r, v);
el.P  = dot(hv,hv)/mu;
hn = hv/norm(hv);
el.hx = -hn(2)/(1+hn(3));                    % tan(i/2)cos(Om)
el.hy =  hn(1)/(1+hn(3));                    % tan(i/2)sin(Om)
ev = cross(v, hv)/mu - r/norm(r);            % Laplace eccentricity vector
s2 = 1 + el.hx^2 + el.hy^2;
fh = [1+el.hx^2-el.hy^2;  2*el.hx*el.hy;      -2*el.hy] / s2;   % equinoctial basis
gh = [2*el.hx*el.hy;      1-el.hx^2+el.hy^2;   2*el.hx] / s2;
el.ex = dot(ev, fh);
el.ey = dot(ev, gh);
el.L  = atan2(dot(r,gh), dot(r,fh));
end
