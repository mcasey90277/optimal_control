% TEST_RUN_GTO_TULIP  No-solve guardrail checks for the campaign front door.
%
% run_gto_tulip is the single entry point for the direct GTO->tulip campaign.
% This pins the preconditions that fire BEFORE any solve, so a bad request
% fails with an explanation instead of a confusing solver error twenty minutes
% in -- and pins the layout invariant the 2026-07-26 flatten established.
%
% Everything here runs in seconds. The solve path is covered by the certified
% reference run and by test_minfuel_lib.
%
% HISTORY. These checks were written for a second, redundant front door added
% under sundman_minfuel/ on 2026-07-26, before it was noticed that run_psr (now
% run_gto_tulip) had always been the real one -- its stage 2 calls the same
% minfuel_at_tf driver. That duplicate is gone; its guards were folded into the
% survivor and this file follows them.
%
% INPUTS:  none
% OUTPUTS: none (throws on any behaviour change)
%
% REFERENCES:
%   [1] ../run_gto_tulip.m (the front door under test).
%   [2] ../README.md (folder layout after the flatten).
% run_gto_tulip is a SCRIPT and shares this workspace, so it overwrites any
% variable it also uses -- `here` and `cfg` among them. Capture everything this
% test needs up front, under names the front door does not touch.
tstDir_ = fileparts(mfilename('fullpath'));
addpath(fullfile(tstDir_, '..'));  setup_paths();
cfgT_   = minfuel_config();
certT_  = fullfile(tstDir_, '..', 'lib', 'sundman_minfuel_certified.mat');

fprintf('\n=== run_gto_tulip guardrail checks ===\n');

%% 1. A factor with no energy backbone is refused, and says what exists ------
% Below ~1.12 this is the OPEN near-min-time band, not a missing run. The
% message has to say so, or the user goes looking for a file to generate that
% cannot be generated.
TULIP_OVERRIDES = struct('factor', 1.05); %#ok<NASGU>
try
    run_gto_tulip;
    error('test_run_gto_tulip:noThrow', 'a factor with no backbone was accepted');
catch ME
    assert(strcmp(ME.identifier, 'run_gto_tulip:noBackbone'), ...
        'wrong error for a missing backbone: %s (%s)', ME.identifier, ME.message);
    assert(contains(ME.message, '1.150'), ...
        'the refusal must list the backbones on disk. Got:\n%s', ME.message);
    assert(contains(ME.message, 'OPEN near-min-time band'), ...
        ['the refusal must say that below ~1.12 this is an open problem, not a ' ...
         'missing run. Got:\n%s'], ME.message);
end
clear TULIP_OVERRIDES
fprintf('  factor    : missing backbone refused, lists what is available\n');

%% 2. The homotopy schedule ends EXACTLY at epsMin, as a prefix -------------
% Truncation (not a fresh schedule) keeps every partial run on the campaign's
% certified step sequence, so an eps-optimal solution is a genuine waypoint of
% the min-fuel run rather than a different homotopy.
for e = [0, 0.05, 0.5]
    base     = cfgT_.schedSharpen;
    effSched = [base(base > e), e];
    assert(abs(effSched(end) - e) < 1e-12, ...
        'epsMin=%g: schedule must END at the request, ends at %g', e, effSched(end));
    assert(all(diff(effSched) < 0), 'epsMin=%g: schedule must be strictly decreasing', e);
    kept = base(base > e);
    assert(isequal(effSched(1:numel(kept)), kept), ...
        'epsMin=%g: schedule is not a prefix of cfgT_.schedSharpen', e);
    fprintf('  epsMin=%-4g: %2d steps, ends exactly at %g, prefix of schedSharpen\n', ...
            e, numel(effSched), effSched(end));
end

%% 3. ONE copy of every shared name on the path -----------------------------
% The flatten dissolved PSR/lib, which had vendored 20 files and left two
% definitions of the solver reachable at once -- with the vendored copies
% silently shadowed, so PSR ran code its own README said it did not. This
% asserts the duplication has not crept back.
for nm = {'casadi_minfuel_sundman','minfuel_at_tf','refine_loop','minfuel_config', ...
          'insertion_states','verify_direct_pmp','warmstart_on_mesh'}
    hits = which(nm{1}, '-all');
    hits = hits(~contains(hits, 'attic'));      % the attic is off-path by design
    assert(numel(hits) == 1, ...
        ['%s resolves to %d copies -- vendoring has crept back:\n  %s\n' ...
         'The flatten deliberately left exactly one of each.'], ...
        nm{1}, numel(hits), strjoin(hits, sprintf('\n  ')));
end
fprintf('  layout    : one copy of each shared name on the path\n');

%% 4. Moved artifacts are where the code now expects them -------------------
assert(isfile(certT_), ...
    ['the certified reference is missing at %s -- it moved into lib/ during the ' ...
     'flatten and several files load it by that path'], certT_);
nBB = numel(dir(fullfile(cfgT_.dirs.energy, 'energy_f*.mat')));
assert(isfolder(cfgT_.dirs.energy) && nBB > 0, ...
    ['no energy backbones under %s -- minfuel_config''s results pointer and the ' ...
     'moved results tree have gone out of sync'], cfgT_.dirs.energy);
fprintf('  artifacts : certified reference in lib/, %d backbones under results/\n', nBB);

fprintf('\ntest_run_gto_tulip: ALL PASS (no solves run)\n');
