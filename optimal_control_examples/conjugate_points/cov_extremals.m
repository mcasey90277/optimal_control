function [E, scan] = cov_extremals(prob, pRange, nScan)
%% Purpose:
%
%   Find EVERY extremal of the two-point problem whose initial slope lies in
%   pRange, by shooting. The shooting function
%
%       r(p) = y(b; p) - yb
%
%   is sampled at nScan slopes; each sign change is refined with fzero. An
%   extremal is a root of r. A two-point problem can have several (the
%   catenary has two) or none, which is why the whole range is scanned
%   rather than one guess polished.
%
%   Note r'(p) = h(b): the slope of the shooting function at a root IS the
%   Jacobi field at the right end. Where r'(p) = 0, b is conjugate to a --
%   that is where extremals are born or die in pairs as b moves (a fold).
%
%% Inputs:
%
%  prob                     struct                  cov_problem output
%  pRange                   [1x2]                   slope interval to scan
%  nScan                    scalar                  scan samples [default 400]
%
%% Outputs:
%
%  E                        struct array            one per extremal, sorted
%                                                   by J: .p .S (cov_shoot to
%                                                   b) .J (the functional)
%                                                   .hb (= r'(p)) .resid
%                                                   (|y(b) - yb|)
%  scan                     struct                  .p [1xn] .r [1xn] (NaN
%                                                   where the flight failed)
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    %Demo: the minimal surface of revolution has two catenaries
    prob = cov_problem('y*sqrt(1+yp^2)', -0.5, 0.5, 1, 1);
    E = cov_extremals(prob, [-8 2]);
    for k = 1:numel(E)
        fprintf('extremal %d: y''(a) = %+.6f, J = %.6f\n', k, E(k).p, E(k).J);
    end
    return
end
if nargin < 3 || isempty(nScan), nScan = 400; end
assert(numel(pRange) == 2 && pRange(2) > pRange(1), 'cov_extremals:range', ...
       'pRange must be [pMin pMax] with pMax > pMin');

b = prob.b;
pv = linspace(pRange(1), pRange(2), nScan);
rv = nan(1, nScan);
for k = 1:nScan
    rv(k) = resid(pv(k));
end
scan = struct('p', pv, 'r', rv);

roots = [];
for k = 1:nScan-1
    r1 = rv(k);  r2 = rv(k+1);
    if ~isfinite(r1) || ~isfinite(r2), continue, end
    if r1 == 0
        roots(end+1) = pv(k);                                         %#ok<AGROW>
    elseif sign(r1) ~= sign(r2)
        roots(end+1) = fzero(@resid, pv([k k+1]), optimset('TolX', 1e-14)); %#ok<AGROW>
    end
end
if isfinite(rv(end)) && rv(end) == 0, roots(end+1) = pv(end); end
roots = uniquetol(roots, 1e-9);

E = struct('p', {}, 'S', {}, 'J', {}, 'hb', {}, 'resid', {});
for k = 1:numel(roots)
    S = cov_shoot(prob, roots(k), b);
    if ~S.ok, continue, end
    zb = deval(S.sol, b);
    E(end+1) = struct('p', roots(k), 'S', S, 'J', cov_functional(prob, S.sol), ...
                      'hb', zb(3), 'resid', abs(zb(1) - prob.yb));     %#ok<AGROW>
end
if ~isempty(E)
    [~, ord] = sort([E.J]);
    E = E(ord);
end

    function r = resid(p)
        % shooting residual y(b; p) - yb, NaN if the flight fails
        Sp = cov_shoot(prob, p, b);
        if Sp.ok, r = Sp.yEnd - prob.yb; else, r = NaN; end
    end
end
