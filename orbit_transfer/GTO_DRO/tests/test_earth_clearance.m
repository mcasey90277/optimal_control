function ok = test_earth_clearance()
%% Purpose:
%
%   earth_clearance on three arcs: (1) a circular 6678.137 km Earth orbit
%   converted with pumpkyn fromPCI reads 300 km altitude to 1e-6 km and
%   PASSES a 200 km floor; (2) the same orbit with a 1000 km floor FAILS
%   (the gate can fail); (3) a straight radial arc that crosses the
%   Earth's radius reports a negative minimum altitude and FAILS; (4) the
%   Moon-centre mode on the tau = 1 DRO reproduces get_family_orbit's
%   perilune to 1 km.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
p = struct('muStar', muStar, 'lStar', lStar);

% (1) circular Earth orbit, radius 6678.137 km, 64 samples, converted at t = 0
rC   = 6678.137;
muE  = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
vC   = sqrt(muE/rC);
th   = linspace(0, 2*pi, 64).';
rvKm = [rC*cos(th), rC*sin(th), zeros(64,1), -vC*sin(th), vC*cos(th), zeros(64,1)];
rvND = zeros(64, 6);
for kk = 1:64
    rvND(kk,:) = pumpkyn.cr3bp.fromPCI(0, rvKm(kk,:), muStar, tStar, lStar, 1);
end
tau = linspace(0, 1, 64).';
[hMin, tAt, pass1] = earth_clearance(tau, rvND, 200, p);
ok = check(ok, abs(hMin - 300) < 1e-6, sprintf('circular 6678.137 km orbit reads %.9f km altitude (300 expected)', hMin));
ok = check(ok, pass1, '300 km orbit passes a 200 km floor');
ok = check(ok, tAt >= 0 && tAt <= 1, 'time of minimum lies in the arc');

% (2) the gate can fail
[~, ~, pass2] = earth_clearance(tau, rvND, 1000, p);
ok = check(ok, ~pass2, '300 km orbit FAILS a 1000 km floor');

% (3) a radial arc through the Earth's centre
xs  = linspace(0.05, -0.05, 21).';                       % ND, crosses x = -muStar
arc = [xs - muStar, zeros(21, 5)];
[hMin3, ~, pass3] = earth_clearance(linspace(0, 1, 21).', arc, 200, p);
ok = check(ok, hMin3 < 0 && ~pass3, sprintf('radial arc through the centre reads %.1f km and FAILS', hMin3));

% (4) Moon-centre mode on the tau = 1 DRO vs get_family_orbit's own geometry
[tA, rvA] = get_family_orbit('dro', struct('tau', 1));
pM = struct('muStar', muStar, 'lStar', lStar, 'center', [1 - muStar, 0, 0], 'radiusKm', 1737.4);
[hMoon, ~, passM, outM] = earth_clearance(tA, rvA, 500, pM);
dMoon = lStar*sqrt(sum((rvA(:,1:3) - [1 - muStar, 0, 0]).^2, 2));
ok = check(ok, abs(hMoon - (min(dMoon) - 1737.4)) < 1, sprintf('Moon mode perilune altitude %.1f km matches the direct computation', hMoon));
ok = check(ok, passM, 'tau = 1 DRO clears a 500 km lunar floor');
ok = check(ok, numel(outM.hKm) == size(rvA, 1), 'altitude history has one entry per sample');

fprintf('test_earth_clearance: %s\n', tern(ok, 'PASS', 'FAIL'));
end

function ok = check(ok, c, msg)
% CHECK  Print one line and fold the verdict.  INPUTS: ok; c; msg.  OUTPUTS: ok.
if c, s = 'ok  '; else, s = 'FAIL'; end
fprintf('  %s %s\n', s, msg);
ok = ok && c;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
