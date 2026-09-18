function ok = test_cartpole_params()
%% Purpose:
%
%   The plant constants have ONE home. Three examples (min-energy, min-time,
%   min-fuel) and their tests share this plant; literals copied into each
%   would be the same failure mode as the duplicated dynamics that hid a sign
%   error until 2026-09-17.
%
%   Checks: the specified values; a struct with exactly the four fields the
%   field function reads, and no more (a fifth field would mean some caller
%   is smuggling a problem-specific quantity through the plant); and that
%   cartpole_field actually accepts what this returns.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
p = cartpole_params();

ok = chk(ok, isstruct(p) && isscalar(p), 'returns a scalar struct');
ok = chk(ok, p.m1 == 5 && p.m2 == 1 && p.L == 2 && p.g == 9.8, ...
         sprintf('constants: m1 %g, m2 %g, L %g, g %g', p.m1, p.m2, p.L, p.g));
ok = chk(ok, isequal(sort(fieldnames(p)), sort({'m1'; 'm2'; 'L'; 'g'})), ...
         'exactly the four fields the plant reads, no more');

[F, G] = cartpole_field([0.1; 0.2; 0.3; 0.4], p);
ok = chk(ok, isequal(size(F), [4 1]) && isequal(size(G), [4 1]) && all(isfinite([F; G])), ...
         'cartpole_field accepts it and returns finite [4 x 1] columns');

if ok, fprintf('TEST_CARTPOLE_PARAMS: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PARAMS: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
