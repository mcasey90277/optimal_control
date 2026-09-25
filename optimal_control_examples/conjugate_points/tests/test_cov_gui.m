function ok = test_cov_gui()
%% Purpose:
%
%   Drive the explorer headless (Visible off) through its public handles
%   and check that what it SAYS matches the instruments:
%
%     1. every preset builds, solves and writes a verdict line;
%     2. catenary: two extremals; the shallow one reads "weak local
%        minimiser", the deep one "NOT a local minimiser" with the Morse
%        count agreeing;
%     3. oscillator past pi: the neighbours' first crossing is pi;
%     4. a malformed Lagrangian is reported as an ERROR in the readout
%        instead of throwing out of the callback.
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
app = conjugate_point_explorer('Visible', 'off');
closer = onCleanup(@() delete(app.fig));

P = cov_presets;
for k = 1:numel(P)
    app.setPreset(k);
    r = app.readout();
    ok = rep(ok, contains(r, 'VERDICT:') && contains(r, 'agree') && ~contains(r, 'DISAGREE'), ...
             sprintf('preset %d builds, verdict, Morse agrees', k), P(k).name);
end

kc = find(strcmp({P.name}, 'Minimal surface (two catenaries)'));
app.setPreset(kc);
s = app.state();
ok = rep(ok, numel(s.E) == 2, 'catenary: two extremals', sprintf('%d', numel(s.E)));
app.selectExtremal(1);
ok = rep(ok, contains(app.readout(), 'a weak local minimiser'), 'catenary shallow: minimiser', '');
app.selectExtremal(2);
ok = rep(ok, contains(app.readout(), 'NOT a local minimiser'), 'catenary deep: not a minimiser', '');

ko = find(strcmp({P.name}, 'Oscillator, b = 4 (past pi)'));
app.setPreset(ko);
app.setDelta(0.3);
s = app.state();
tc = arrayfun(@(f) firstOr(f.tCross), s.fan);
% numel first: all() of an EMPTY list is true, which once let this pass
ok = rep(ok, numel(tc) == 4 && all(abs(tc - pi) < 1e-6), ...
         'oscillator: all 4 neighbours cross at pi', mat2str(tc, 8));

% difference panel: raw zeros are the crossings; scaled by 1/delta it IS h
% on a linear problem, for every delta
app.setScaled(false);
s = app.state();
zOK = true;
for kk = 1:numel(s.fan)
    dz = interp1(s.diff.t, s.diff.D(kk,:), s.fan(kk).tCross(1));
    zOK = zOK && abs(dz) < 1e-6*max(abs(s.diff.D(kk,:)));
end
ok = rep(ok, zOK && size(s.diff.D, 1) == 4, 'difference: zero at each crossing', '');
app.setScaled(true);
s = app.state();
inAB = s.diff.t <= s.prob.b;
errOsc = max(abs(s.diff.D(:, inAB) - s.diff.h(inAB)), [], 'all');
ok = rep(ok, s.diff.scaled && errOsc < 1e-7, 'oscillator: (y_d - y0)/d = h exactly', ...
         sprintf('max err %.1e', errOsc));

% nonlinear: the scaled difference approaches h at first order in delta
kp = find(strcmp({P.name}, 'Pendulum (nonlinear)'));
app.setPreset(kp);  app.selectExtremal(2);  app.setScaled(true);
s = app.state();
inAB = s.diff.t <= s.prob.b;
errP = max(abs(s.diff.D(:, inAB) - s.diff.h(inAB)), [], 2).';
rat = errP(1:end-1)./errP(2:end);
ok = rep(ok, all(diff(errP) < 0) && abs(rat(end) - 2) < 0.2, ...
         'pendulum: (y_d - y0)/d -> h, first order', ...
         sprintf('errors %s, ratios %s', mat2str(errP, 3), mat2str(rat, 3)));
app.setScaled(false);
app.setPreset(ko);

f = findall(app.fig, 'Type', 'uieditfield', 'Value', P(ko).F);
f(1).Value = 'yp^2 - 4*y^2';                 % Jacobi field sin(2t)/2: t_c = pi/2
app.solve();
s = app.state();
ok = rep(ok, abs(s.S.tConj(1) - pi/2) < 1e-8 && contains(app.readout(), 'NOT a local minimiser'), ...
         'typed F = yp^2 - 4y^2: t_c = pi/2', sprintf('%.10f', s.S.tConj(1)));

f(1).Value = 'y^2 - yp^2';                   % P = -2: Legendre's necessary condition fails
app.solve();
r = app.readout();
ok = rep(ok, contains(r, 'necessary condition fails') && contains(r, 'not applicable') ...
             && ~contains(r, 'DISAGREE'), 'P < 0: Legendre verdict, Morse marked n/a', '');

f(1).Value = 'yp^2 - q';                     % q is not an allowed variable
app.solve();
ok = rep(ok, startsWith(app.readout(), 'ERROR'), 'bad Lagrangian reported, not thrown', '');
fprintf('test_cov_gui: %s\n', pf(ok));
end

% ---------------------------------------------------------------------------
function ok = rep(ok, c, name, msg)
% REP  Print one check and fold it into the verdict.
% INPUTS: ok running verdict, c this check, name, msg. OUTPUTS: ok.
fprintf('  %-44s %s   %s\n', name, pf(c), msg);
ok = ok && c;
end

function t = firstOr(v)
% FIRSTOR  First element, or NaN when empty. INPUTS: v. OUTPUTS: t scalar.
if isempty(v), t = NaN; else, t = v(1); end
end

function s = pf(ok)
% PF  PASS/FAIL text. INPUTS: ok logical. OUTPUTS: s char.
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
