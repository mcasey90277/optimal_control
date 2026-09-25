function ok = test_cov_mutation()
%% Purpose:
%
%   Deliberate-failure experiments: corrupt one piece of the oscillator
%   problem and require the check aimed at it to FAIL. (A check that
%   passes on a broken input is not a check -- repo rule from the cart-pole
%   mutation sweep, 2026-09-17.)
%
%     M1  Jacobi equation with the wrong sign of g_y (h'' = +h): h becomes
%         sinh, never zero -> the 'jacobi' check must fail.
%     M2  Legendre coefficient P quartered: the second variation gains
%         negative modes -> the 'second_variation' check must fail.
%     M3  Lagrangian integrand scaled by 1.01 in the functional only: the
%         'functional' check must fail (J no longer cot(b)).
%
%   Each mutant must also leave the checks it does NOT touch passing, so a
%   failure is attributed to the right instrument.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 every mutant was caught
%                                                   by its check and only by
%                                                   the checks it touches
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

addpath(fileparts(fileparts(mfilename('fullpath'))));
ok = true;
base = cov_problem('yp^2 - y^2', 0, 4, 0, 1);

m1 = base;  m1.gy = @(t, y, yp) -base.gy(t, y, yp);
m2 = base;  m2.P  = @(t, y, yp) base.P(t, y, yp)/4;
m3 = base;  m3.F  = @(t, y, yp) 1.01*base.F(t, y, yp);
mutants = {m1, 'jacobi',           {'extremal', 'functional', 'second_variation'}
           m2, 'second_variation', {'extremal', 'functional', 'jacobi', 'neighbours'}
           m3, 'functional',       {'extremal', 'jacobi', 'neighbours'}};
for k = 1:size(mutants, 1)
    chk = oscillator_checks(mutants{k,1});
    target = mutants{k,2};
    caught = isfield(chk, target) && ~chk.(target).ok;
    others = mutants{k,3};
    clean = all(cellfun(@(f) isfield(chk, f) && chk.(f).ok, others));
    fprintf('  M%d  %-18s caught: %s   others clean: %s\n', k, target, yn(caught), yn(clean));
    ok = ok && caught && clean;
end
fprintf('test_cov_mutation: %s\n', pf(ok));
end

% ---------------------------------------------------------------------------
function s = yn(b)
% YN  yes/no text. INPUTS: b logical. OUTPUTS: s char.
if b, s = 'yes'; else, s = 'NO'; end
end

function s = pf(ok)
% PF  PASS/FAIL text. INPUTS: ok logical. OUTPUTS: s char.
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
