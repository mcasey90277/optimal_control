function ok = test_ladder_ms_hardcap()
%% Purpose:
%
%   Tests opts.msHardCapSec of thrust_ladder_library -- the hard fence on
%   each ms_tfmin refinement attempt and tfMin acceptance (run_capped on a
%   parallel worker) -- and the mass-depletion source guard of
%   mintime_prop_seg. Engine checks run on the stage 1 GTO -> DRO demo cell
%   (orient 0, DRO tau 1, departure phase 0.10, arrival phase 0, Earth floor
%   200 km; the cell of test_ladder_seedtop):
%     (a) malformed options refuse before any solve: msHardCapSec -1 and
%         Inf, msFcn not a function handle
%     (a2) mintime_prop_seg refuses a segment whose exact end-mass lower
%          bound m0 - (Tmax/c) dt is <= 0 (mintime_prop_seg:massDepleted)
%     (b) THE CAP FIRES: through the test seam opts.msFcn the refinement
%         sleeps 1000 s; with msHardCapSec = 20 the engine returns, the log
%         carries 'ms capped after <s> s' with <s> in [20, 55] (cap +
%         run_capped's verified-cancel window), the rung is recorded as NOT
%         refined (OK false, Z8 / ACCDZ / RES NaN) while the direct-solve
%         fields are written, and with msSeg [12 24] the capped K = 12
%         attempt ends the msSeg loop (no K = 24 attempt)
%     (c) DEFAULT-OFF PIN + IDENTITY: rungs [15 12] with the option off
%         reproduce the PRE-CHANGE engine's measured t_f to 1e-12 ND
%         (0.304038463224718 at 15 N, 0.307000714257882 at 12 N; measured
%         2026-11-10 from HEAD's thrust_ladder_library, where all of
%         TF/FLYKM/ACCDZ/RES/OK/Z8 were isequaln), and a generous cap
%         (600 s) agrees with option off: same OK, |dTF| < 1e-10 ND,
%         |dZ8| < 1e-8 (measured differences are printed)
%     (d) A WORKER ERROR CARRIES ITS CAUSE: a seam that throws on the worker
%         is logged with its identifier and message, and the second
%         consecutive worker error throws thrust_ladder_library:fenceBroken
%   The session's pool state is left as found: a pool the engine opened is
%   deleted when the engine returns, an existing one is reused and kept.
%   Needs the Parallel Computing Toolbox and a free PCT licence. Without a
%   pool the seams refuse to run in-process, so the test errors rather than
%   hangs.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true when every check
%                                                   passes (errors otherwise)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
if isempty(which('oc.ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts(here))), 'oclib'));
end
if isempty(which('pumpkyn.cr3bp.tfMin'))
    od = cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(od);
end
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

outDir = fullfile(tempdir, 'test_ladder_ms_hardcap');
if isfolder(outDir), rmdir(outDir, 's'); end
mkdir(outDir);
cleanDir = onCleanup(@() rmdir_quiet(outDir));

hadPool = ~isempty(current_pool());

base = struct('rungs', [15 12], 'ispS', 1710, 'm0kg', 150, 'N', 400, ...
    'floorKm', 500, 'floorEarthKm', 200, 'maxIter', 6000, 'maxCpuSec', 600, ...
    'gateKm', 100, 'accTol', 1e-6, 'nD', 1, 'nA', 1, 'sD0', 0.10, 'sA0', 0, ...
    'depFamily', 'gto', 'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'dro', 'arrParams', struct('tau', 1), 'tauDRO', 0, ...
    'thrLock', true, 'tf0', 0.30, 'resume', false);

% (a) malformed option refuses before any solve
o = base;  o.msHardCapSec = -1;  o.logFile = fullfile(outDir, 'a.log');
try
    thrust_ladder_library(fullfile(outDir, 'a.mat'), o);
    error('test:noThrow', 'msHardCapSec = -1 did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:msHardCapSec'), ...
        'wrong error for a negative cap: %s', err.identifier);
end
o.msHardCapSec = Inf;
try
    thrust_ladder_library(fullfile(outDir, 'a.mat'), o);
    error('test:noThrow', 'msHardCapSec = Inf did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:msHardCapSec'), ...
        'wrong error for an infinite cap: %s', err.identifier);
end
o.msHardCapSec = 60;  o.msFcn = 'ms_tfmin';
try
    thrust_ladder_library(fullfile(outDir, 'a.mat'), o);
    error('test:noThrow', 'a non-handle msFcn did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:msFcn'), ...
        'wrong error for a non-handle msFcn: %s', err.identifier);
end
assert(~isfile(fullfile(outDir, 'a.mat')), 'a refused call must not write the sheet');

% (a2) the source guard: Tmax/c = 2 per ND time, m0 = 0.1 -> bound 0.1 - 2*0.06 < 0
y0 = [0.5; 0; 0; 0; 0.5; 0; 0.1; ones(7, 1)];
try
    mintime_prop_seg(0.06, y0, false, 2, 1, 0.0121505856);
    error('test:noThrow', 'a mass-depleting segment did not refuse');
catch err
    assert(strcmp(err.identifier, 'mintime_prop_seg:massDepleted'), ...
        'wrong error for a mass-depleting segment: %s', err.identifier);
end

% (b) the cap fires on a refinement that never returns
assert(~isempty(ver('parallel')), 'test_ladder_ms_hardcap needs the Parallel Computing Toolbox');
capS = 20;
o = base;  o.rungs = 15;  o.msSeg = [12 24];  o.msHardCapSec = capS;
% the seam sleeps ONLY on a worker; run in-process (no pool, engine fell back
% to unfenced) it throws instead, so this test can never hang the session
o.msFcn = @(varargin) deal(evalc(['if isempty(getCurrentTask()), ' ...
    'error(''test:unfenced'', ''refinement ran in-process''), end, pause(1000)']));
o.logFile = fullfile(outDir, 'b.log');
tB = tic;
Pb = thrust_ladder_library(fullfile(outDir, 'b.mat'), o);
wallB = toc(tB);
logB = fileread(o.logFile);
assert(contains(logB, 'ms hard cap 20 s ACTIVE'), ...
    'no pool: the fence did not open (log:\n%s)', logB);
tok = regexp(logB, 'K=12  ms capped after (\d+) s', 'tokens', 'once');
assert(~isempty(tok), 'no "ms capped after" line in the log:\n%s', logB);
cappedS = str2double(tok{1});
fprintf('  cap %d s: capped attempt took %d s, whole engine call %.1f s (direct solve + pool + fence)\n', ...
    capS, cappedS, wallB);
assert(cappedS >= capS && cappedS <= capS + 35, 'capped attempt took %d s for a %d s cap', cappedS, capS);
assert(~contains(logB, 'pool lost'), 'the cancelled worker did not stop within 30 s');
assert(~contains(logB, 'K=24'), 'a capped attempt must end the msSeg loop for the rung (K = 24 was tried)');
assert(wallB < 900, 'engine call took %.0f s: the sleeping refinement was not fenced', wallB);
assert(~Pb.OK(1, 1, 1), 'a capped refinement must not mark the rung OK');
assert(all(isnan(Pb.Z8(:, 1, 1, 1))) && isnan(Pb.ACCDZ(1, 1, 1)) && isnan(Pb.RES(1, 1, 1)), ...
    'a capped refinement must leave Z8 / ACCDZ / RES empty');
assert(isfinite(Pb.TF(1, 1, 1)) && isfinite(Pb.FLYKM(1, 1, 1)), ...
    'the direct-solve fields of the rung must still be recorded');
assert(~isempty(current_pool()) == hadPool, 'pool state changed by the engine call (b)');

% (c) identity: option off vs a generous cap
o = base;  o.logFile = fullfile(outDir, 'off.log');
Poff = thrust_ladder_library(fullfile(outDir, 'off.mat'), o);
o = base;  o.msHardCapSec = 600;  o.logFile = fullfile(outDir, 'on.log');
Pon = thrust_ladder_library(fullfile(outDir, 'on.mat'), o);
logOn = fileread(o.logFile);
assert(contains(logOn, 'ms hard cap 600 s ACTIVE'), 'the generous-cap call did not run fenced');
assert(~contains(logOn, 'capped after') && ~contains(logOn, 'worker error'), ...
    'the generous-cap call had a failed fenced attempt:\n%s', logOn);
dTF = max(abs(Pon.TF(:) - Poff.TF(:)));
dZ8 = max(abs(Pon.Z8(:) - Poff.Z8(:)));
dAcc = max(abs(Pon.ACCDZ(:) - Poff.ACCDZ(:)));
dFly = max(abs(Pon.FLYKM(:) - Poff.FLYKM(:)));
fprintf('  off vs cap 600 s: OK %s / %s, max|dTF| = %.3g ND, max|dZ8| = %.3g, max|dACCDZ| = %.3g, max|dFLYKM| = %.3g km\n', ...
    mat2str(squeeze(Poff.OK).'), mat2str(squeeze(Pon.OK).'), dTF, dZ8, dAcc, dFly);
assert(all(Poff.OK(:)), 'both rungs must be OK on the demo cell with the option off');
tfPre = [0.304038463224718 0.307000714257882];   % pre-change engine, measured 2026-11-10
dPre = max(abs(squeeze(Poff.TF(1, 1, :)).' - tfPre));
fprintf('  option off vs the pre-change engine: max|dTF| = %.3g ND; Z8(15 N) = %s\n', ...
    dPre, mat2str(Poff.Z8(:, 1, 1, 1).', 16));
assert(dPre < 1e-12, 'option-off t_f moved from the pre-change engine by %.3g ND', dPre);
assert(all(Poff.ACCDZ(:) < 1e-6) && all(Poff.FLYKM(:) < 100), 'option-off acceptance / flown bounds');
assert(isequal(Poff.OK, Pon.OK), 'OK differs between option off and a generous cap');
assert(dTF < 1e-10, 'TF differs by %.3g ND', dTF);
assert(dZ8 < 1e-8, 'Z8 differs by %.3g', dZ8);
assert(~isempty(current_pool()) == hadPool, 'pool state changed by the engine call (c)');

% (d) a worker error is logged with its cause; the second consecutive one is
% infrastructure and throws fenceBroken. In-process (no pool) the seam's own
% error would surface instead, failing the identifier check below.
o = base;  o.rungs = 15;  o.msSeg = [12 24];  o.msHardCapSec = 60;
o.msFcn = @(varargin) deal(evalc('error(''test:seamBoom'', ''injected worker failure'')'));
o.logFile = fullfile(outDir, 'd.log');
try
    thrust_ladder_library(fullfile(outDir, 'd.mat'), o);
    error('test:noThrow', 'two consecutive worker errors did not throw');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:fenceBroken'), ...
        'wrong error after two worker errors: %s (%s)', err.identifier, err.message);
end
logD = fileread(o.logFile);
fprintf('  worker-error log lines:\n%s', strjoin(regexp(logD, '[^\n]*worker error[^\n]*\n', 'match'), ''));
assert(contains(logD, 'K=12  ms worker error after') && contains(logD, 'test:seamBoom: injected worker failure'), ...
    'the first worker error must be logged with its cause:\n%s', logD);
assert(contains(err.message, 'injected worker failure'), 'fenceBroken must carry the last cause');
assert(~isempty(current_pool()) == hadPool, 'pool state changed by the engine call (d)');

ok = true;
fprintf('test_ladder_ms_hardcap: PASS\n');
end

% ------------------------------------------------------------------------
function rmdir_quiet(d)
% RMDIR_QUIET  Remove a temp folder if present.  INPUTS: d char.  OUTPUTS: none.
if isfolder(d)
    try rmdir(d, 's'); catch, end
end
end
