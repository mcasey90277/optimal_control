function ok = test_local_residual()
% TEST_LOCAL_RESIDUAL  Unit test for oc.local_residual (G1-gate engine).
%
%   Three checks on a linear system dz = A z whose flow is exactly
%   expm(A dt):
%     1. EXACT NODES: nodes sampled from the true flow give per-interval
%        residuals at integrator tolerance (< 1e-10), NOT zero (that
%        would mean nothing was integrated).
%     2. INJECTED ERROR: corrupting one interior node by delta shows up
%        in the two adjacent intervals at amplitude ~delta, and nowhere
%        else -- the residual is LOCAL.
%     3. SHAPES: dX is [nx x N], one column per interval.
%
% INPUTS:  none
% OUTPUTS: ok - all checks pass [logical]

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                         % oclib root -> +oc visible
ok = true;

A = [0 1; -2 -0.3];                               % damped oscillator
N = 8;  tG = linspace(0, 3, N+1);
X = zeros(2, N+1);  X(:,1) = [1; 0];
for k = 1:N
    X(:,k+1) = expm(A*(tG(k+1)-tG(k))) * X(:,k);  % exact flow
end
rhs = @(t, z) A*z;

% 1. exact nodes -> tolerance-level residual:
dX = oc.local_residual(X, tG, rhs);
ok = check('shape [nx x N]', isequal(size(dX), [2 N])) && ok;
ok = check(sprintf('exact nodes: max|dX| %.1e < 1e-10', max(abs(dX(:)))), ...
           max(abs(dX(:))) < 1e-10) && ok;
ok = check('not identically zero (really integrated)', any(dX(:) ~= 0)) && ok;

% 2. injected error is local and has the right size:
kc = 4;  delta = 1e-6;
Xc = X;  Xc(1,kc+1) = Xc(1,kc+1) + delta;         % corrupt node kc+1
dXc = oc.local_residual(Xc, tG, rhs);
r = max(abs(dXc), [], 1);
ok = check(sprintf('interval kc residual ~ delta (%.2e)', r(kc)), ...
           abs(r(kc) - delta) < 0.01*delta) && ok;
ok = check('interval kc+1 sees the corrupted start', r(kc+1) > 0.1*delta) && ok;
far = r([1:kc-1, kc+2:end]);
ok = check(sprintf('all other intervals clean (max %.1e)', max(far)), ...
           max(far) < 1e-10) && ok;
fprintf('test_local_residual: %s\n', string(ok));
end

function ok = check(name, cond)
% CHECK  Print one gate line.  INPUTS: name;cond.  OUTPUTS: ok.
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
