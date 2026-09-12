function ok = test_validate_flight()
%% Purpose:
%
%   Tests validate_flight -- ONE flight validator shared by the certifier
%   and the study script, so "the flight is admissible" means the same thing
%   everywhere: it reached t_f, its data are finite, mass obeys the all-burn
%   law and stays positive, and nothing came closer to a primary than the
%   stated clearance. A returned array is not a completed flight, and finite
%   positive-mass samples can describe an integration that stopped early.
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
mu = 0.012150585609624;  lStar = 389703.264829278;
Tnd = 0.1756418;  cnd = 8.6737;  tf = 4.0;

% a synthetic all-burn flight far from both primaries
t = linspace(0, tf, 400).';
r = [0.5 + 0.1*t, 0.8*ones(size(t)), 0.05*t];   v = ones(400, 3);
m = 1 - (Tnd/cnd)*t;
Y = [r, v, m, zeros(400, 7)];

V = validate_flight(t, Y, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, V.ok, sprintf('a clean flight validates (%s)', V.reason));
ok = chk(ok, V.moonKm > 1e4 && V.earthKm > 1e4, ...
         sprintf('clearances reported: Moon %.0f km, Earth %.0f km', V.moonKm, V.earthKm));

% stopped early: t_end < t_f
V2 = validate_flight(t(1:300), Y(1:300, :), tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V2.ok && contains(V2.reason, 't_f'), sprintf('early stop refused: %s', V2.reason));

% mass law violated
Yb = Y;  Yb(:, 7) = 1 - 0.5*(Tnd/cnd)*t;
V3 = validate_flight(t, Yb, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V3.ok && contains(V3.reason, 'mass'), sprintf('wrong mass law refused: %s', V3.reason));

% a sample inside the Moon
Yc = Y;  Yc(200, 1:3) = [1 - mu, 0, 0] + 1000/lStar*[1 0 0];
V4 = validate_flight(t, Yc, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V4.ok && contains(lower(V4.reason), 'moon'), sprintf('lunar penetration refused: %s', V4.reason));

% non-finite data
Yd = Y;  Yd(50, 4) = NaN;
V5 = validate_flight(t, Yd, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V5.ok && contains(lower(V5.reason), 'finite'), 'non-finite data refused');

% (Astra review 2026-09-11) the mass law must hold at EVERY sample, not
% just the last one: an interior glitch with a correct endpoint
Ye = Y;  Ye(150:160, 7) = Ye(150:160, 7) + 1e-3;
V6 = validate_flight(t, Ye, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V6.ok && contains(V6.reason, 'mass law'), sprintf('interior mass glitch refused: %s', V6.reason));
ok = chk(ok, isfield(V, 'massLawErr') && V.massLawErr < 1e-12, ...
         sprintf('and the clean flight''s pointwise mass-law error is reported: %.1e', V.massLawErr));

% costates carried and non-finite
Yf = Y;  Yf(70, 12) = Inf;
V7 = validate_flight(t, Yf, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V7.ok && contains(V7.reason, 'costate'), sprintf('non-finite costate refused: %s', V7.reason));

% all-burn to a t_f that exhausts the mass: the EXPECTED final mass is
% non-positive, and no reported mass may rescue it
tfx = 1.05*cnd/Tnd;  tx = linspace(0, tfx, 400).';
Yx = [r, v, max(1 - (Tnd/cnd)*tx, 1e-3), zeros(400, 7)];
V8 = validate_flight(tx, Yx, tfx, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V8.ok && contains(V8.reason, 'exhausts'), sprintf('exhaustion refused: %s', V8.reason));

% non-monotone time
tb = t;  tb(200) = tb(199) - 1e-3;
V9 = validate_flight(tb, Y, tf, Tnd, cnd, mu, lStar);
ok = chk(ok, ~V9.ok && contains(V9.reason, 'times'), sprintf('non-monotone time refused: %s', V9.reason));

if ok, fprintf('TEST_VALIDATE_FLIGHT: ALL PASS\n'); else, fprintf('TEST_VALIDATE_FLIGHT: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
