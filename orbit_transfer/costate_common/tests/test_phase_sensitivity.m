function ok = test_phase_sensitivity()
% TEST_PHASE_SENSITIVITY  phase_sensitivity and phase_edge_residuals against
% ORACLES THAT ARE NOT THE CODE UNDER TEST:
%
%   phase_sensitivity   the minimum-time single integrator xdot = u, |u| <= 1
%       in R^3, whose answer is known in closed form: T = |x_f - x_0|, the
%       costate is the constant unit vector lambda = -(x_f - x_0)/T (so that
%       lambda . f = -1 along the extremal, the library's convention). With
%       both endpoints sliding along curves, dT/ds is differentiated
%       NUMERICALLY from T itself and compared with the costate formula.
%       Nothing orbital is involved, so a sign or a transpose cannot hide
%       behind the CR3BP.
%   phase_edge_residuals   a smooth periodic T(s) with its exact derivative
%       (the residual must sit at the trapezoid rule's third-order bound), a
%       planted jump (it must be flagged at exactly its two edges), both
%       array dimensions, the wrap across the seam, a non-uniform list.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));           % costate_common

% ---- the single-integrator oracle ------------------------------------------
x0 = @(s) [cos(2*pi*s); sin(2*pi*s); 0.3*s];                   % departure curve
xf = @(s) [4 + 0.5*cos(2*pi*s); 1 + 2*sin(2*pi*s); -1 + s^2];  % arrival curve
T  = @(s0, sf) sqrt(sum((xf(sf) - x0(s0)).^2));
s0 = 0.31;  sf = 0.74;  e = 1e-6;
lam = -(xf(sf) - x0(s0))/T(s0, sf);                            % constant along the extremal
dx0 = (x0(s0 + e) - x0(s0 - e))/(2*e);   dxf = (xf(sf + e) - xf(sf - e))/(2*e);
[g0, gf] = phase_sensitivity(lam, lam, dx0, dxf);
n0 = (T(s0 + e, sf) - T(s0 - e, sf))/(2*e);   nf = (T(s0, sf + e) - T(s0, sf - e))/(2*e);
ok = chk(ok, abs(g0 - n0) < 1e-8, sprintf('dT/ds_departure: costate formula %+.8f, numerical %+.8f', g0, n0));
ok = chk(ok, abs(gf - nf) < 1e-8, sprintf('dT/ds_arrival  : costate formula %+.8f, numerical %+.8f', gf, nf));
% the sign, stated in words: carry the arrival point directly AWAY and the time grows
away = (xf(sf) - x0(s0))/T(s0, sf);
[~, gAway] = phase_sensitivity(lam, lam, dx0, away);
ok = chk(ok, abs(gAway - 1) < 1e-12, 'moving the arrival point straight away costs one unit of time per unit distance');
ok = chk(ok, throws(@() phase_sensitivity(lam(1:2), lam, dx0, dxf)), 'mismatched sizes are refused');
ok = chk(ok, throws(@() phase_sensitivity([lam(1:2); NaN], lam, dx0, dxf)), 'a non-finite costate is refused');

% ---- edge residuals: smooth, exact gradient -> the trapezoid bound ----------
n = 24;  s = (0:n-1)/n;  h = 1/n;
Tn = 2 + 0.3*sin(2*pi*s);  Gn = 0.3*2*pi*cos(2*pi*s);
r = phase_edge_residuals(Tn, Gn, s, 2);
bound = h^3/12 * 0.3*(2*pi)^3;                                 % (h^3/12) max|T'''|
ok = chk(ok, isequal(size(r), size(Tn)) && max(abs(r)) <= 1.02*bound && max(abs(r)) > 0.5*bound, ...
         sprintf('smooth T: max residual %.2e sits at the third-order bound %.2e', max(abs(r)), bound));
% ---- a planted jump is flagged at its two edges, and only there -------------
Tj = Tn;  Tj(10) = Tj(10) + 0.05;
rj = phase_edge_residuals(Tj, Gn, s, 2);
big = find(abs(rj) > 0.01);
ok = chk(ok, isequal(big, [9 10]), sprintf('a 0.05 jump at cell 10 flags edges %s (9 and 10)', mat2str(big)));
% ---- the wrap: the last edge joins cell n to cell 1 -------------------------
Tw = Tn;  Tw(1) = Tw(1) + 0.05;
rw = phase_edge_residuals(Tw, Gn, s, 2);
ok = chk(ok, isequal(find(abs(rw) > 0.01), [1 n]), 'a jump at cell 1 flags edge 1 and the SEAM edge n -> 1');
% ---- along dimension 1 of a matrix, and NaN cells -----------------------------
M = repmat(Tn(:), 1, 3);  GM = repmat(Gn(:), 1, 3);  M(5, 2) = NaN;
rm = phase_edge_residuals(M, GM, s, 1);
ok = chk(ok, isequal(size(rm), [n 3]) && all(isnan(rm([4 5], 2))) && max(abs(rm(:, 1))) <= 1.02*bound, ...
         'dimension 1: same bound per column; an empty cell blanks its two edges');
% ---- a non-uniform list: each edge uses its own circular spacing -------------
su = sort(mod(0.013 + [0 0.03 0.1 0.18 0.3 0.41 0.5 0.62 0.7 0.81 0.9 0.97], 1));
Tu = 2 + 0.3*sin(2*pi*su);  Gu = 0.3*2*pi*cos(2*pi*su);
ru = phase_edge_residuals(Tu, Gu, su, 2);
hmax = max(mod(diff([su, su(1) + 1]), 1));
ok = chk(ok, max(abs(ru)) <= 1.02*hmax^3/12*0.3*(2*pi)^3, 'non-uniform list: every edge within the bound of the widest spacing');

if ok, fprintf('test_phase_sensitivity: ALL PASS\n'); else, fprintf('test_phase_sensitivity: FAIL\n'); end
end

function tf = throws(f)
% THROWS  Does calling f throw?  INPUTS: f (handle).  OUTPUTS: tf.
tf = false;
try, f(); catch, tf = true; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
