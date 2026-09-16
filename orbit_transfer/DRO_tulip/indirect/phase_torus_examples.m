%% PHASE_TORUS_EXAMPLES  Three ways to call run_phase_torus, from one entry to a new problem.
%
%   Set `example` and `live`, run the script. Every case builds a spec in
%   the shape section 0 of run_phase_torus reads (grid, departure orbit,
%   arrival orbit, engine, anchors, folder, budgets), prints the plan, and
%   -- only with live = true -- launches the campaign. Plan mode costs
%   nothing and shows exactly what would run.
%
%     example = 1   ONE ENTRY. One departure phase, one arrival phase: the
%                   arcs walk from the anchor to that arrival phase, the
%                   sheet certifies the crossing, there is nothing to rib
%                   (one departure phase), the entry is packaged, audited
%                   and swept. For a single transfer at 70 mN the front door
%                   run_dro_tulip(sD, sA) is the cheaper tool; this shows the
%                   torus machinery at its smallest.
%     example = 2   THE 24 x 24 LIBRARY of FINDINGS 60-72 as one call: the
%                   shipped 70 mN DRO(tau 1) -> 7-petal torus, started from
%                   the five anchors the campaign found. Live, this is a
%                   day or two of machine time.
%     example = 3   ANOTHER PROBLEM: DRO of period 2.0 -> 8-petal tulip at
%                   the same engine. No catalog holds a root for it (the
%                   shipped 1710 s catalogs cover tau {0.5,1,2,3} x Np
%                   {5,7,9,12} at 15-0.5 N), so step 1 MAKES the first anchor
%                   by a direct solve warm-started from the 7-petal anchor
%                   (anchor_by_direct_solve, minutes), then a 6 x 6 torus
%                   runs from it.
%
%   Everything here is a parameter of the case; nothing is hard-coded in
%   the driver. FINDINGS 73; process/PHASE_TORUS_RUNBOOK.md section 1b.
%
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.

if ~exist('example', 'var'), example = 1; end    % 1 | 2 | 3   (set before running to choose)
if ~exist('live', 'var'),    live = false;  end   % false: plan only (prints, launches nothing); true: run it

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
res = fullfile(here, 'results');
sevenPetalAnchor = fullfile(res, 'mintime_70mN_anchor.mat');   % DRO tau 1 -> 7-petal, 70 mN, sA 0.0754, 17.80 d

switch example
%% ========================================================================
%  1. ONE ENTRY: (sD, sA) = (0, 0.4087)
%% ========================================================================
case 1
    spec = struct( ...
        'sD',        0, ...                                     % one departure phase: no ribs
        'sA',        0.0754 + 8/24, ...                         % one arrival phase (0.4087)
        'departure', struct('family', 'dro', 'tau', 1.0), ...   % the DRO, period 1.0 ND
        'arrival',   struct('family', 'tulip', 'Np', 7, 'pm', -1), ...
        'engine',    struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150), ...
        'anchors',   {{'anchor', sevenPetalAnchor, 0.0754, 'fast'}}, ...
        'outDir',    fullfile(res, 'torus_one_entry'), 'tag', 'one', ...
        'arc',       struct('nStep', 800, 'deadlineSec', 3600, 'span', 0.45), ...   % just far enough to cross 0.4087
        'nWorkers',  1, 'maxRounds', 2, 'discover', false);

%% ========================================================================
%  2. THE 24 x 24 LIBRARY (the shipped 70 mN torus, from its five anchors)
%% ========================================================================
case 2
    anchors = { ...                                                       % name       file                                                 sA       label
        'fast',     fullfile(res, 'mintime_70mN_anchor.mat'),           0.0754,          'fast';    ...
        'A2',       fullfile(res, 'mintime_70mN_certified.mat'),        0.0754 + 10/12,  'A2';      ...
        'fast2',    fullfile(res, 'mintime_70mN_anchor_fast2.mat'),     0.8671,          'fast2';   ...
        'direct18', fullfile(res, 'mintime_70mN_anchor_direct18.mat'),  0.0754 + 17/24,  'direct18'; ...
        'direct11', fullfile(res, 'mintime_70mN_anchor_direct11.mat'),  0.0754 + 10/24,  'direct11'};
    spec = struct( ...
        'sD',        (0:23)/24, ...
        'sA',        sort(mod(0.0754 + (0:23)/24, 1)), ...       % the 24-lattice from the first anchor's phase, wrapped
        'departure', struct('family', 'dro', 'tau', 1.0), ...
        'arrival',   struct('family', 'tulip', 'Np', 7, 'pm', -1), ...
        'engine',    struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150), ...
        'anchors',   {anchors}, ...
        'outDir',    fullfile(res, 'torus_70mN_24x24'), 'tag', '70mN', ...
        'nWorkers',  4, 'maxRounds', 6, 'improveDays', 2, 'slowDays', 2);

%% ========================================================================
%  3. ANOTHER PROBLEM: DRO tau 2.0 -> 8-petal tulip, 70 mN, 6 x 6 phases
%% ========================================================================
case 3
    departure = struct('family', 'dro', 'tau', 2.0);
    arrival   = struct('family', 'tulip', 'Np', 8, 'pm', -1);
    engine    = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150);
    sA0 = 0.0754;                                             % the first anchor's arrival phase
    anchorFile = fullfile(res, 'mintime_70mN_anchor_dro2_tulip8.mat');
    % -- step 1: the first anchor, by a direct solve from the 7-petal root
    if ~isfile(anchorFile)
        if live
            [C, info] = anchor_by_direct_solve(sevenPetalAnchor, 0, sA0, ...
                struct('tauDRO', departure.tau, 'NpTulip', arrival.Np, 'pmTulip', arrival.pm), engine, anchorFile);
            assert(C.ok, 'no anchor for DRO 2.0 -> 8-petal: %s (direct %.1f d). Try another sA0 or seed.', C.reason, info.tfDirectDays);
        else
            fprintf(['example 3, plan: the anchor %s does not exist yet; live = true makes it by a direct solve\n' ...
                     '  warm-started from %s (anchor_by_direct_solve), then runs the torus from it.\n'], anchorFile, sevenPetalAnchor);
        end
    end
    spec = struct( ...
        'sD',        (0:5)/6, ...
        'sA',        sort(mod(sA0 + (0:5)/6, 1)), ...
        'departure', departure, 'arrival', arrival, 'engine', engine, ...
        'anchors',   {{'a0', anchorFile, sA0, 'first'}}, ...
        'outDir',    fullfile(res, 'torus_dro2_tulip8_6x6'), 'tag', 'dro2t8', ...
        'nWorkers',  4, 'maxRounds', 4);
    if ~live && ~isfile(anchorFile)
        fprintf('example 3, plan: with the anchor in place the call would be run_phase_torus(spec) with\n');  disp(spec);  return
    end
otherwise
    error('phase_torus_examples: example must be 1, 2 or 3');
end

%% ========================================================================
%  4. RUN: plan first (always), then live if asked
%% ========================================================================
planSpec = spec;  planSpec.plan = true;
out = run_phase_torus(planSpec);                         % the plan, printed
fprintf('example %d: plan state "%s"\n', example, out.state);
if live
    out = run_phase_torus(spec);
    fprintf('example %d: %s (%s) -> %s\n', example, out.state, out.reason, out.final);
end
