function ok = test_nd_propulsion()
%% Purpose:
%
%   Tests nd_propulsion -- the nondimensional propulsion conversion that
%   about twenty files each wrote out by hand:
%
%     g0  = 9.80665*tStar^2/(1000*lStar);
%     cnd = (ispS/tStar)*g0;
%     Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);
%
%   One file's comment already said it used "the same ND thrust/exhaust
%   conversion thrust_ladder_library uses", which is how a copy announces
%   itself.
%
%   Checks:
%     1. the three quantities equal the inline expressions BITWISE (the
%        migration must not move a single ulp anywhere);
%     2. the ndT closure equals Tnd at the same thrust, and scales linearly
%        across a ladder of rungs;
%     3. the physics is right: c = Isp*g0 gives the exhaust speed back in
%        km/s, and the ND thrust acceleration re-dimensionalises to T/m;
%     4. malformed input is refused by name.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
lStar = 389703.264829278;  tStar = 382981.289129055;

cases = [0.070 900 150; 15 1710 150; 0.025 1600 15];
worst = 0;
for k = 1:size(cases, 1)
    thrustN = cases(k,1);  ispS = cases(k,2);  m0kg = cases(k,3);
    g0  = 9.80665*tStar^2/(1000*lStar);
    cnd = (ispS/tStar)*g0;
    Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);
    p = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
    same = isequal(p.g0, g0) && isequal(p.cnd, cnd) && isequal(p.Tnd, Tnd);
    worst = max(worst, max(abs([p.g0 - g0, p.cnd - cnd, p.Tnd - Tnd])));
    ok = chk(ok, same, sprintf('%.3g N / Isp %g / %g kg reproduces the inline formulas bitwise', ...
             thrustN, ispS, m0kg));
end
ok = chk(ok, worst == 0, sprintf('worst difference over the cases: %g', worst));

p = nd_propulsion(0.070, 900, 150, lStar, tStar);
ok = chk(ok, isequal(p.ndT(0.070), p.Tnd), 'ndT at the built thrust equals Tnd');
% linearity, judged RELATIVELY: rescaling through a different expression
% costs rounding, ~1 ulp, so an absolute bar here would test floating-point
% exactness rather than linearity
rungs = [15 12 10 7 5 1];
a = arrayfun(@(T) p.ndT(T), rungs);   b = rungs*(p.Tnd/0.070);
relLin = max(abs(a - b)./abs(a));
ok = chk(ok, relLin < 4*eps, sprintf('ndT is linear in thrust across a rung ladder (rel %.1e < 4 eps)', relLin));

% the LADDER form: a rung sweep has an Isp and a mass but no single thrust,
% so thrustN may be omitted -- c and ndT are still defined, Tnd is not
q = nd_propulsion([], 900, 150, lStar, tStar);
ok = chk(ok, isequal(q.cnd, p.cnd) && isequal(q.ndT(0.070), p.Tnd) && isnan(q.Tnd), ...
         'with no thrust given: same c and ndT, Tnd is NaN');
ok = chk(ok, isnan(q.thrustN), 'and the echoed thrust is NaN, not a silent default');

% the Isp CLOSURE: probe_abstract_case continues over Isp, so c at an
% arbitrary Isp must come from the same place as everything else
ok = chk(ok, isequal(p.ndC(900), p.cnd), 'ndC at the built Isp equals cnd');
p1710 = nd_propulsion(0.07, 1710, 150, lStar, tStar);
ok = chk(ok, isequal(p.ndC(1710), p1710.cnd), 'ndC at another Isp equals a build at that Isp');

% NEITHER given: arclength_thrust receives c as an argument and only needs
% the thrust closure, so both may be omitted and both scalars are then NaN
r = nd_propulsion([], [], 150, lStar, tStar);
ok = chk(ok, isnan(r.Tnd) && isnan(r.cnd) && isequal(r.ndT(0.070), p.Tnd) && isequal(r.ndC(900), p.cnd), ...
         'with neither thrust nor Isp: both scalars NaN, both closures still exact');

cKmS = p.cnd*lStar/tStar;
ok = chk(ok, abs(cKmS - 900*9.80665/1000) < 1e-9, ...
         sprintf('c redimensionalises to Isp*g0 = %.4f km/s', cKmS));
aMs2 = p.Tnd*lStar*1000/tStar^2;
ok = chk(ok, abs(aMs2 - 0.070/150) < 1e-15, ...
         sprintf('T_nd redimensionalises to T/m = %.3e m/s^2', aMs2));

ok = chk(ok, refuses(@() nd_propulsion(-1, 900, 150, lStar, tStar), 'nd_propulsion:positive'), ...
         'a negative thrust is refused');
ok = chk(ok, refuses(@() nd_propulsion(0.07, 900, 0, lStar, tStar), 'nd_propulsion:positive'), ...
         'a zero mass is refused');
ok = chk(ok, refuses(@() nd_propulsion([0.07 0.08], 900, 150, lStar, tStar), 'nd_propulsion:scalar'), ...
         'a non-scalar thrust is refused (use ndT for a ladder)');

if ok, fprintf('TEST_ND_PROPULSION: ALL PASS\n'); else, fprintf('TEST_ND_PROPULSION: FAIL\n'); end
end

function r = refuses(fh, id)
% REFUSES  True when fh throws the named identifier.  INPUTS: fh; id.
% OUTPUTS: r [logical].
try
    fh();  r = false;
catch ME
    r = strcmp(ME.identifier, id);
    if ~r, fprintf('        (threw %s: %s)\n', ME.identifier, ME.message); end
end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
