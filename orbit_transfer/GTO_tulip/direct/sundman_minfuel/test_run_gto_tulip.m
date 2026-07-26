% TEST_RUN_GTO_TULIP  No-solve guardrail checks for the front door.
%
% run_gto_tulip's job is as much about what it REFUSES as what it solves. Two
% of its parameters look adjustable but name open campaign problems, and a
% front door that accepted them would hand the user a solve that is expected to
% fail, with a confusing solver error instead of an explanation. This test pins
% those refusals, and the schedule the epsMin knob builds.
%
% Everything here runs in SECONDS -- every check trips a precondition that fires
% before the first solve. The solve path itself is exercised by an actual run
% (see the reference case in the commit that added this file).
%
% NOTE ON HOW THE SCHEDULE IS OBSERVED. run_gto_tulip is a script, so it shares
% this workspace. It builds `sched` BEFORE the energy-backbone existence check,
% so requesting a factor with no backbone leaves a fully-built `sched` behind
% when it errors. That ordering is what makes the schedule observable without a
% solve; if the checks are ever reordered, this test's schedule section will
% stop seeing `sched` and must be revisited rather than deleted.
%
% INPUTS:  none
% OUTPUTS: none (throws on any behaviour change)
%
% REFERENCES:
%   [1] run_gto_tulip.m (the front door under test).
%   [2] process/LADDER_PREP_PILOT_FINDINGS.md (the thrust topology wall).
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths;
cfg = minfuel_config();

fprintf('\n=== run_gto_tulip guardrail checks ===\n');

%% 1. Off-nominal thrust is refused, and the message explains WHY -----------
TULIP_OVERRIDES = struct('thrustN', 0.020); %#ok<NASGU>
try
    run_gto_tulip;
    error('test_run_gto_tulip:noThrow', 'off-nominal thrust was accepted');
catch ME
    assert(strcmp(ME.identifier, 'run_gto_tulip:offNominalThrust'), ...
        'wrong error for off-nominal thrust: %s', ME.identifier);
    % The point is not that it throws -- it is that the message routes the user
    % to the open problem instead of leaving them to guess.
    assert(contains(ME.message, 'pilot_rung_20mN') && contains(ME.message, 'topology wall'), ...
        ['the off-nominal-thrust refusal must name the pilot script and the ' ...
         'topology wall, else the user cannot act on it. Got:\n%s'], ME.message);
end
clear TULIP_OVERRIDES
fprintf('  thrust     : off-nominal refused, message names the open problem\n');

%% 2. A factor with no energy backbone is refused and lists what exists -----
TULIP_OVERRIDES = struct('factor', 1.05); %#ok<NASGU>
try
    run_gto_tulip;
    error('test_run_gto_tulip:noThrow', 'a factor with no backbone was accepted');
catch ME
    assert(strcmp(ME.identifier, 'run_gto_tulip:noBackbone'), ...
        'wrong error for a missing backbone: %s', ME.identifier);
    assert(contains(ME.message, '1.150') && contains(ME.message, 'gen_tulip_energy_2p'), ...
        ['the missing-backbone refusal must list the backbones on disk and name ' ...
         'the generator. Got:\n%s'], ME.message);
end
clear TULIP_OVERRIDES
fprintf('  factor     : missing backbone refused, lists what is available\n');

%% 3. epsMin builds a schedule that is a PREFIX of the certified one -------
% Truncation matters: it keeps every partial run on the same certified step
% sequence, so an eps-optimal solution is a genuine waypoint of the min-fuel
% run rather than a different homotopy.
for e = [0, 0.05, 0.5]
    clear sched
    TULIP_OVERRIDES = struct('factor', 1.05, 'epsMin', e); %#ok<NASGU>
    try, run_gto_tulip; catch, end     % expected to stop at the backbone check
    clear TULIP_OVERRIDES
    assert(exist('sched','var') == 1, ...
        ['sched was not built before the backbone check for epsMin=%g -- the ' ...
         'precondition order changed; see this file''s header note.'], e);
    assert(abs(sched(end) - e) < 1e-12, ...
        'epsMin=%g: schedule must END at the request, ends at %g', e, sched(end));
    assert(all(diff(sched) < 0), 'epsMin=%g: schedule must be strictly decreasing', e);
    kept = cfg.schedSharpen(cfg.schedSharpen > e);
    assert(isequal(sched(1:numel(kept)), kept), ...
        ['epsMin=%g: schedule is not a prefix of cfg.schedSharpen -- partial runs ' ...
         'must follow the same certified step sequence'], e);
    fprintf('  epsMin=%-4g : %2d steps, ends exactly at %g, prefix of schedSharpen\n', ...
            e, numel(sched), sched(end));
end

% epsMin = 1 is the min-ENERGY request: a single tight eps=1 solve, no sharpening.
clear sched
TULIP_OVERRIDES = struct('factor', 1.05, 'epsMin', 1); %#ok<NASGU>
try, run_gto_tulip; catch, end
clear TULIP_OVERRIDES
assert(isequal(sched, 1), 'epsMin=1 must give the single-step schedule [1], got %s', mat2str(sched));
fprintf('  epsMin=1    : single tight eps=1 step (min-energy, no sharpening)\n');

%% 4. epsMin outside [0,1] is rejected --------------------------------------
for bad = [-0.1, 1.5, NaN]
    TULIP_OVERRIDES = struct('factor', 1.05, 'epsMin', bad); %#ok<NASGU>
    threw = false;
    try, run_gto_tulip; catch ME, threw = contains(ME.identifier, 'run_gto_tulip'); end
    clear TULIP_OVERRIDES
    assert(threw, 'epsMin=%g should have been rejected by validateattributes', bad);
end
fprintf('  epsMin      : values outside [0,1] rejected\n');

fprintf('\ntest_run_gto_tulip: ALL PASS (no solves run)\n');
