function ok = test_cov_oscillator()
%% Purpose:
%
%   The harmonic-oscillator Lagrangian F = y'^2 - y^2 on [0, b], y(0) = 0,
%   y(b) = 1, checked against its CLOSED-FORM answers -- an oracle that
%   shares no code with what it tests:
%
%     extremal         y0(t) = sin(t)/sin(b),  y0'(0) = 1/sin(b)
%     functional       J = y0(b) y0'(b) = cot(b)   (integrate by parts)
%     Jacobi field     h(t) = sin(t): conjugate points at k pi
%     neighbours       linear problem: y0 + delta sin(t) exactly, so every
%                      neighbour crosses at k pi for ANY delta
%     second variation 2 int (eta'^2 - eta^2): eigenvalues 2((k pi/b)^2 - 1),
%                      negative count = #{k : k pi < b}
%     Delta J          exactly quadratic in eps (F is quadratic)
%
%   Run for b = 3 (no conjugate point), 4 (one) and 7 (two).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 all checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

addpath(fileparts(fileparts(mfilename('fullpath'))));
ok = true;
for b = [3 4 7]
    prob = cov_problem('yp^2 - y^2', 0, b, 0, 1);
    chk = oscillator_checks(prob);
    names = fieldnames(chk);
    for k = 1:numel(names)
        c = chk.(names{k});
        fprintf('  b = %d  %-34s %s   (%s)\n', b, names{k}, pf(c.ok), c.msg);
        ok = ok && c.ok;
    end
end
fprintf('test_cov_oscillator: %s\n', pf(ok));
end

% ---------------------------------------------------------------------------
function s = pf(ok)
% PF  PASS/FAIL text. INPUTS: ok logical. OUTPUTS: s char.
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
