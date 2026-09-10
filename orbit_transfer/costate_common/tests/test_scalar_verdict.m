function ok = test_scalar_verdict()
%% Purpose:
%
%   Tests scalar_verdict -- the guard that stops a MALFORMED external result
%   from passing a gate by accident.
%
%   Why it exists (Astra chain review, 2026-09-10): every one of these
%   slips through a plain `if` in MATLAB.
%     it.conj.pass = []      ->  `if C.conj ~= 1` is FALSE on empty, so an
%                                empty verdict behaves like a pass
%     g.minLamV = Inf        ->  `Inf > 0` is true, so a degenerate gate
%                                reads as satisfied
%     it.normR = NaN         ->  every comparison is false, so a NaN
%                                residual never trips a bound
%     a vector verdict       ->  `if` on a vector means ALL, which is not
%                                the test anyone intended
%   A gate must read a real finite SCALAR or refuse to interpret it at all.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

good = {1, 0, -2.5, true, int8(3), single(1.5)};
for k = 1:numel(good)
    [g, v] = scalar_verdict(good{k});
    ok = chk(ok, g && isa(v, 'double') && v == double(good{k}), ...
             sprintf('accepts %s %s -> %g', class(good{k}), mat2str(good{k}), v));
end

bad = {[], NaN, Inf, -Inf, [1 1], zeros(0,1), 'a', {1}, 1+2i, struct('a',1)};
names = {'empty', 'NaN', 'Inf', '-Inf', 'vector', '0x1', 'char', 'cell', 'complex', 'struct'};
for k = 1:numel(bad)
    [g, v] = scalar_verdict(bad{k});
    ok = chk(ok, ~g && isnan(v), sprintf('refuses %s', names{k}));
end

% the interface that matters: a refused verdict must not compare as a pass
[g, v] = scalar_verdict([]);
ok = chk(ok, ~(g && v == 1), 'an empty verdict cannot read as "pass"');
[g, v] = scalar_verdict(Inf);
ok = chk(ok, ~(g && v > 0), 'an infinite gate cannot read as "positive"');

if ok, fprintf('TEST_SCALAR_VERDICT: ALL PASS\n'); else, fprintf('TEST_SCALAR_VERDICT: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
