function ok = test_gto_dro_entry()
%% Purpose:
%
%   gto_dro_entry's PLANNING and REFUSALS, without a solve: (1) thrustN
%   outside [1, 15] is refused by name; (2) the ladder plan is
%   [15] for 15 N, [15 12 10 7] for 7 N, all nine rungs for 1 N and
%   [15 12 10 7 5 4] for 4 N (a non-standard rung is appended, never
%   cold-solved); (3) the GTO parameter passthrough and the DRO period
%   land in the engine options the driver builds (planOnly returns them in
%   entry.engineOpts); (3b) a GTO perigee altitude below floorEarthKm is
%   refused by name before any solve. The solve itself is exercised end to end by
%   test_gto_dro_study.
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
po = struct('planOnly', true);

% (1) refusals
for bad = [0.5, 16, -1]
    threw = false;
    try
        gto_dro_entry(0, 1, 0, 0, bad, po);
    catch err
        threw = strcmp(err.identifier, 'gto_dro_entry:thrust');
    end
    ok = check(ok, threw, sprintf('thrustN = %g refused with gto_dro_entry:thrust', bad));
end

% (2) ladder plans
plans = {15, [15];  7, [15 12 10 7];  1, [15 12 10 7 5 3 2 1.5 1];  4, [15 12 10 7 5 4]};
for kk = 1:size(plans, 1)
    e = gto_dro_entry(0, 1, 0, 0, plans{kk, 1}, po);
    ok = check(ok, isequal(e.rungsWalked, plans{kk, 2}), ...
               sprintf('thrustN = %g plans [%s]', plans{kk, 1}, num2str(e.rungsWalked)));
end

% (3) passthrough: GTO params, DRO period, phases, recipe reach the engine options
e = gto_dro_entry(90, 2, 0.25, 0.5, 10, struct('planOnly', true, ...
        'gtoParams', struct('sma_km', 25000, 'ecc', 0.7), 'tf0', 0.4, 'm0kg', 120, 'ispS', 2000));
eo = e.engineOpts;
ok = check(ok, strcmp(eo.depFamily, 'gto') && strcmp(eo.arrFamily, 'dro'), 'families gto -> dro');
ok = check(ok, eo.depParams.orientDeg == 90 && eo.depParams.sma_km == 25000 && abs(eo.depParams.ecc - 0.7) < 1e-15, 'GTO params reach depParams');
ok = check(ok, eo.arrParams.tau == 2, 'DRO period reaches arrParams.tau');
ok = check(ok, eo.sD0 == 0.25 && eo.sA0 == 0.5 && eo.nD == 1 && eo.nA == 1, 'phases reach a 1x1 grid');
ok = check(ok, eo.tf0 == 0.4 && eo.m0kg == 120 && eo.ispS == 2000 && eo.thrLock, 'recipe, mass and Isp reach the engine');
ok = check(ok, eo.tauDRO == 90, 'sheet key carries orientDeg (the GTO catalog rule)');
ok = check(ok, isequal(eo.rungs, [15 12 10]), 'planned rungs are the engine rungs');
ok = check(ok, eo.floorEarthKm == 200, 'Earth floor default 200 km reaches the engine');
eE = gto_dro_entry(0, 1, 0, 0, 15, struct('planOnly', true, 'floorEarthKm', 300));
ok = check(ok, eE.engineOpts.floorEarthKm == 300, 'Earth floor override 300 km reaches the engine');

% (3b) preflight: a GTO perigee below the Earth floor is refused by name, planOnly included
threw = false;
try
    gto_dro_entry(0, 1, 0, 0, 15, struct('planOnly', true, 'gtoParams', struct('sma_km', 6500, 'ecc', 0)));
catch err
    threw = strcmp(err.identifier, 'gto_dro_entry:earthFloor');
end
ok = check(ok, threw, 'GTO perigee 122 km < floorEarthKm 200 km refused with gto_dro_entry:earthFloor');

% (4) scratch name carries the phases (parallel cells must not collide); fractional orient is not mangled
e1 = gto_dro_entry(22.5, 1, 0.1, 0.2, 10, po);
e2 = gto_dro_entry(22.5, 1, 0.1, 0.3, 10, po);
[~, n1] = fileparts(e1.scratchMat);
ok = check(ok, ~isempty(strfind(n1, 'sD0.1000')) && ~isempty(strfind(n1, 'sA0.2000')) && ~isempty(strfind(n1, 'o22.5_')), ...
           sprintf('scratch name %s carries orient and both phases', n1));
ok = check(ok, ~strcmp(e1.scratchMat, e2.scratchMat), 'cells differing only in arrival phase get distinct scratch files');
ok = check(ok, strcmp([e1.scratchMat(1:end-4) '.log'], e1.engineOpts.logFile), 'log path is the scratch name with .log');

fprintf('test_gto_dro_entry: %s\n', tern(ok, 'PASS', 'FAIL'));
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
