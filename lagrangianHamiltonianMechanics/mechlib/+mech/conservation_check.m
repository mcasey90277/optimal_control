function chk = conservation_check(t, Y, sys, bars)
%% Purpose:
%
%   Energy drift and constraint residual along a solution of a system
%   built by mech.lagrange_dalembert, each against a bar, with one verdict
%   and a reason. Prints nothing; the caller prints.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (used for the
%                                                   length check only)
%
%  Y                        [2n x N]                Columns [q; qd]
%
%  sys                      struct                  From
%                                                   mech.lagrange_dalembert
%                                                   (.n, .energy,
%                                                   .constraint)
%
%  bars                     struct (optional)       .energy [1e-9] relative
%                                                   drift bar; .constraint
%                                                   [1e-9] max |A qd| bar
%
%% Outputs:
%
%  chk                      struct                  .E [1 x N]; .driftRel
%                                                   = max |E - E(1)| /
%                                                   max(|E(1)|, eps);
%                                                   .resid = max over the
%                                                   grid of max |A qd|;
%                                                   .residSeries [1 x N];
%                                                   .okEnergy;
%                                                   .okConstraint; .ok;
%                                                   .why (empty on a pass)
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a particle on a circle, exact, then with an injected drift
     syms x y xd yd real
     sys = mech.lagrange_dalembert([x; y], [xd; yd], (xd^2 + yd^2)/2, [x, y]);
     t = linspace(0, 2*pi, 201);
     Y = [cos(t); sin(t); -sin(t); cos(t)];
     chk = mech.conservation_check(t, Y, sys);
     fprintf('exact circle: drift %.1e, residual %.1e, ok = %d\n', chk.driftRel, chk.resid, chk.ok);
     Y(3:4, :) = Y(3:4, :).*(1 + 1e-6*t/t(end));          % a ramp: a uniform scale would rescale E(1) too and show nothing
     chk = mech.conservation_check(t, Y, sys);
     fprintf('ramped velocity: drift %.1e, ok = %d, why = %s\n', chk.driftRel, chk.ok, chk.why);
     return
end

if nargin < 4 || isempty(bars), bars = struct(); end
bE = getf(bars, 'energy', 1e-9);
bC = getf(bars, 'constraint', 1e-9);
n  = sys.n;
N  = size(Y, 2);
if size(Y, 1) ~= 2*n || numel(t) ~= N
    error('mech:conservation_check:shape', 'Y must be [2n x N] with N = numel(t) (n = %d; got Y %d x %d, numel(t) = %d)', ...
          n, size(Y, 1), size(Y, 2), numel(t));
end

E     = zeros(1, N);
resid = zeros(1, N);
for kk = 1:N
    qk  = Y(1:n, kk);
    qdk = Y(n+1:2*n, kk);
    E(kk) = sys.energy(qk, qdk);
    c = sys.constraint(qk, qdk);
    if ~isempty(c), resid(kk) = max(abs(c)); end
end

chk.E            = E;
chk.driftRel     = max(abs(E - E(1)))/max(abs(E(1)), eps);
chk.residSeries  = resid;
chk.resid        = max(resid);
chk.okEnergy     = chk.driftRel <= bE;
chk.okConstraint = chk.resid <= bC;
chk.ok           = chk.okEnergy && chk.okConstraint;
chk.why          = '';
if ~chk.okEnergy
    chk.why = sprintf('energy drift %.2e > %.1e', chk.driftRel, bE);
end
if ~chk.okConstraint
    if isempty(chk.why), sep = ''; else, sep = '; '; end
    chk.why = sprintf('%s%sconstraint residual %.2e > %.1e', chk.why, sep, chk.resid, bC);
end
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
%% Purpose:
%
%   Field with a default.
%
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
