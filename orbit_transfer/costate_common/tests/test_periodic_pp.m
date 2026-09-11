function ok = test_periodic_pp()
%% Purpose:
%
%   Tests periodic_pp -- the C1-PERIODIC cubic through one period of a
%   closed orbit, with its seam mismatch reported as a number. It replaces
%   two private copies that disagreed on policy: transfer_study's
%   periodicPP (errors without the Curve Fitting Toolbox) and
%   arclength_arrival's makePP (silently fell back to an ordinary spline,
%   whose derivative JUMPS at the seam -- exactly where dR/dsA is needed).
%
%   Checks, on an analytic circle where value and derivative are known:
%     1. it interpolates;
%     2. the derivative closure is the derivative;
%     3. the seam is continuous in value AND derivative;
%     4. scheme 'notaknot' visibly breaks the derivative seam -- the reason
%        the periodic form is required, and proof the option acts;
%     5. seam.value reports the DATA's own closure, which no interpolant
%        can improve;
%     6. malformed input is refused by name.
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

th = linspace(0, 2*pi, 61).';              % one period, first sample = last
y  = [cos(th), sin(th)];
[pp, seam, dpp] = periodic_pp(th, y);

% Interpolation error, checked by its ORDER rather than against a magic
% number: a cubic on m intervals is O(h^4), so halving h must cut the error
% by about 16. (The first version of this check asserted < 1e-7 at 60
% intervals, where the cubic's own error is ~3e-7 -- the threshold was
% wrong, not the interpolant.)
q  = [0.37, 1.9, 4.4, 6.1];                % off-node test points
e1 = max(max(abs(ppval(pp, q) - [cos(q); sin(q)])));
th2 = linspace(0, 2*pi, 121).';
pp2 = periodic_pp(th2, [cos(th2), sin(th2)]);
e2 = max(max(abs(ppval(pp2, q) - [cos(q); sin(q)])));
ok = chk(ok, e1 < 1e-6 && e2 < e1, sprintf('interpolates the circle: %.1e at 60 intervals, %.1e at 120', e1, e2));
ok = chk(ok, e1/e2 > 8 && e1/e2 < 32, sprintf('and converges at fourth order: error ratio %.1f (16 expected)', e1/e2));
d = ppval(dpp, q);
ok = chk(ok, max(max(abs(d - [-sin(q); cos(q)]))) < 1e-5, ...
         sprintf('the derivative closure is the derivative, to %.1e', max(max(abs(d - [-sin(q); cos(q)])))));
ok = chk(ok, seam.value < 1e-12 && seam.deriv < 1e-12, ...
         sprintf('seam continuous: value %.1e, derivative %.1e', seam.value, seam.deriv));

[~, seamNK] = periodic_pp(th, y, struct('scheme', 'notaknot'));
ok = chk(ok, seamNK.deriv > 1e-4, ...
         sprintf('an ordinary spline JUMPS at the seam: derivative mismatch %.1e', seamNK.deriv));

yOpen = y;  yOpen(end,:) = yOpen(end,:) + 1e-6;     % a table that does not close
[~, seamOpen] = periodic_pp(th, yOpen);
ok = chk(ok, abs(seamOpen.value - norm([1e-6 1e-6])) < 1e-9, ...
         sprintf('seam.value reports the data''s own closure (%.2e)', seamOpen.value));

ok = chk(ok, refuses(@() periodic_pp(th, y(1:end-1,:)), 'periodic_pp:size'), ...
         'a table whose length does not match the time vector is refused');
ok = chk(ok, refuses(@() periodic_pp(flipud(th), y), 'periodic_pp:tIncreasing'), ...
         'a non-increasing time vector is refused');
ok = chk(ok, refuses(@() periodic_pp(th(1:3), y(1:3,:)), 'periodic_pp:tooFew'), ...
         'fewer than four samples is refused');

if ok, fprintf('TEST_PERIODIC_PP: ALL PASS\n'); else, fprintf('TEST_PERIODIC_PP: FAIL\n'); end
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
