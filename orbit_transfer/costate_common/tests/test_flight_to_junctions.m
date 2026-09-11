function ok = test_flight_to_junctions()
%% Purpose:
%
%   Tests flight_to_junctions -- a flown trajectory cut into the K+1
%   junction states a multiple-shooting seed needs. Six engines wrote the
%   same two lines inline (`interp1(tu/tu(end), yj(iu,1:14), sGrid,
%   'pchip')'` plus, for a rung change, the all-burn mass row), and
%   seed_from_z8 held a seventh copy that queried in ABSOLUTE time.
%
%   Checks:
%     1. it reproduces the samples at the junction times;
%     2. duplicate propagator times are dropped, and the result is the same
%        as if they had never been there;
%     3. the mass row, when asked for, is the all-burn identity
%        m = 1 - T t / c EXACTLY at every junction -- derived, never
%        rescaled from the old trajectory (a rescale carried a scaling
%        mistake once; see extend_thrust_ladder's comment);
%     4. without it the flight's own mass row survives;
%     5. the grid spans [0, tf], with tf overridable for a rung-scaled seed;
%     6. malformed input is refused by name;
%     7. EQUIVALENCE: on a real 70 mN flight it reproduces the inline
%        formula the engines used, BITWISE.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

% a synthetic "flight": 14 rows whose exact values are known
tj = linspace(0, 3, 301).';
yj = [sin(tj), cos(tj), tj, tj.^2, sin(2*tj), cos(2*tj), 1 - 0.1*tj, ...
      repmat(tj, 1, 7)];
K = 12;
[Y, tGrid] = flight_to_junctions(tj, yj, K);
ok = chk(ok, isequal(size(Y), [14 K+1]) && numel(tGrid) == K+1, ...
         sprintf('shape [14 x %d] and a matching grid', K+1));
ok = chk(ok, abs(tGrid(end) - 3) < 1e-15 && tGrid(1) == 0, 'the grid spans [0, tf]');
err = max(abs(Y(1,:) - sin(tGrid)));
ok = chk(ok, err < 1e-6, sprintf('it interpolates the flight (row 1 error %.1e)', err));

dup = [tj(1:100); tj(100); tj(101:end)];  ydup = [yj(1:100,:); yj(100,:); yj(101:end,:)];
Ydup = flight_to_junctions(dup, ydup, K);
ok = chk(ok, isequal(Ydup, Y), 'duplicate propagator times change nothing');

Tnd = 0.1756;  cnd = 3.3;  tf = 2.0;
Ym = flight_to_junctions(tj, yj, K, struct('tf', tf, 'massLaw', struct('Tnd', Tnd, 'cnd', cnd)));
mExact = 1 - Tnd*linspace(0, tf, K+1)/cnd;
ok = chk(ok, isequal(Ym(7,:), mExact), 'the mass row is the all-burn identity, exactly');
ok = chk(ok, isequal(Y(7,:), interp1(tj/tj(end), yj(:,7), linspace(0,1,K+1), 'pchip')), ...
         'without a mass law the flight''s own mass row survives');

ok = chk(ok, refuses(@() flight_to_junctions(tj, yj(1:end-1,:), K), 'flight_to_junctions:size'), ...
         'a state table that does not match the time vector is refused');
ok = chk(ok, refuses(@() flight_to_junctions(tj, yj, 0), 'flight_to_junctions:segments'), ...
         'K < 1 is refused');

% 7. the equivalence gate, on a real flight
lib = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', 'mintime_70mN_anchor.mat');
if isfile(lib)
    L = load(lib);  z8 = L.z(:);
    ndp = nd_propulsion(0.070, 900, 150);  mu = 0.012150585609624;
    B = arclength_arrival('setup');  rv0 = B.stateD(0);
    [tp, yp] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], ndp.Tnd, ndp.cnd, mu);
    [tu, iu] = unique(tp);
    sGrid = linspace(0, 1, 25);
    Yold = interp1(tu/tu(end), yp(iu,1:14), sGrid, 'pchip')';      % the engines' inline form
    Ynew = flight_to_junctions(tp, yp, 24);
    ok = chk(ok, isequal(Yold, Ynew), 'reproduces the engines'' inline formula BITWISE on a real flight');
else
    fprintf('  SKIP  no 70 mN anchor on disk for the equivalence gate\n');
end

if ok, fprintf('TEST_FLIGHT_TO_JUNCTIONS: ALL PASS\n'); else, fprintf('TEST_FLIGHT_TO_JUNCTIONS: FAIL\n'); end
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
