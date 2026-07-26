% TEST_CERTIFIED_GUARD  Pin every branch of the shared certified-quantity guard.
%
% certified_guard replaced four hand-written copies across four campaigns. Those
% copies stood between "verified the certified point" and "verified some other
% point and said so confidently", so the extraction is only safe if every
% refusal branch still fires on exactly the same inputs.
%
% This pins all three refusals, BOTH quantity directions, and -- the branch that
% matters most and is easiest to get backwards -- that a BETTER re-solve is
% accepted-and-flagged rather than refused.
%
% INPUTS:  none
% OUTPUTS: none (throws on any behaviour change)
%
% REFERENCES:
%   [1] verify_common/certified_guard.m (the function under test).
%   [2] orbit_transfer/CODE_STRUCTURE.md (Tier-1 extraction record).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

ok = @(varargin) struct('success', true, 'ipoptStatus', 'Solve_Succeeded', ...
                        'maxDefect', 1e-12, 'value', varargin{1});

% min-fuel: maximize m_f (higher is better) -- the tulip/ELFO/earth/CR3BP case
FUEL = struct('caller','run_foc_tulip','label','row10N','saved',1.0, ...
              'name','m_f','errName','mass','better','higher');
% min-time: minimize t_f (lower is better) -- the ELFO Route-B anchor case
TIME = struct('caller','run_foc_elfo','label','anchor','saved',6.0962, ...
              'name','t_f','errName','tf','better','lower','units',' ND');

fprintf('\n=== certified_guard branch gate ===\n');

expectThrow = @(fn, wantId, what) assertThrows(fn, wantId, what);

%% 1. success class ---------------------------------------------------------
% Both a false .success and a non-whitelisted status must refuse. The second is
% the case the Sundman solvers need: they set success=true on any non-throwing
% solve, so the status whitelist is the only thing catching a bad exit.
r = ok(1.0);  r.success = false;
expectThrow(@() certified_guard(r, FUEL), 'run_foc_tulip:resolveFailed', 'success=false');

r = ok(1.0);  r.ipoptStatus = 'Maximum_Iterations_Exceeded';
expectThrow(@() certified_guard(r, FUEL), 'run_foc_tulip:resolveFailed', 'bad status');

r = ok(1.0);  r.ipoptStatus = 'Infeasible_Problem_Detected';
expectThrow(@() certified_guard(r, FUEL), 'run_foc_tulip:resolveFailed', 'infeasible status');

% 'Solved_To_Acceptable_Level' is explicitly ALLOWED -- it is a real converged
% point and several certified rows exit this way.
r = ok(1.0);  r.ipoptStatus = 'Solved_To_Acceptable_Level';
info = certified_guard(r, FUEL);
assert(~info.improved, 'acceptable-level exit should pass unflagged');
fprintf('  success class : refuses false/bad/infeasible, admits acceptable-level\n');

%% 2. feasibility -----------------------------------------------------------
r = ok(1.0);  r.maxDefect = 1e-6;                       % > default feasTol 1e-8
expectThrow(@() certified_guard(r, FUEL), 'run_foc_tulip:notFeasible', 'loose defect');

r = ok(1.0);  r.maxDefect = 1e-8;                       % exactly AT tol -> allowed
info = certified_guard(r, FUEL);
assert(isequal(info.maxDefect, 1e-8), 'maxDefect not reported through');

% An explicit feasTol must override the default.
r = ok(1.0);  r.maxDefect = 1e-6;
s = FUEL;  s.feasTol = 1e-5;
info = certified_guard(r, s);
assert(info.maxDefect == 1e-6, 'explicit feasTol not honoured');
fprintf('  feasibility   : refuses loose defect, honours explicit feasTol\n');

%% 3. certified quantity, HIGHER-is-better ---------------------------------
% Degraded: mass dropped beyond tolerance -> refuse.
expectThrow(@() certified_guard(ok(1.0 - 1e-3), FUEL), ...
    'run_foc_tulip:massDegraded', 'mass dropped');

% Improved: mass rose -> ACCEPT and flag. Refusing here would reject the
% legitimate "warm re-solve found a nearby better optimum" case that the earth
% campaign documents at 5 N and 10 N.
info = certified_guard(ok(1.0 + 1e-3), FUEL);
assert(info.improved, 'a higher m_f must be flagged improved, not refused');
assert(abs(info.relChange - 1e-3) < 1e-12, 'relChange wrong: %g', info.relChange);

% Unchanged (within tol): pass, not flagged.
info = certified_guard(ok(1.0 + 1e-9), FUEL);
assert(~info.improved, 'a within-tolerance change must not be flagged improved');
fprintf('  m_f (higher)  : refuses drop, flags rise as improved, ignores noise\n');

%% 4. certified quantity, LOWER-is-better ----------------------------------
% The direction that inverts. t_f going UP is the degradation for min-time.
expectThrow(@() certified_guard(ok(6.0962 + 1e-3), TIME), ...
    'run_foc_elfo:tfDegraded', 't_f rose');

info = certified_guard(ok(6.0962 - 1e-3), TIME);
assert(info.improved, 'a lower t_f must be flagged improved, not refused');

info = certified_guard(ok(6.0962), TIME);
assert(~info.improved, 'identical t_f must not be flagged');
fprintf('  t_f (lower)   : refuses rise, flags drop as improved\n');

%% 5. direction must be declared, not guessed ------------------------------
% A silently-defaulted direction would invert the entire gate.
s = FUEL;  s.better = 'bigger';
expectThrow(@() certified_guard(ok(1.0), s), 'certified_guard:badSpec', 'bad direction');
fprintf('  spec          : rejects an unknown .better rather than defaulting\n');

%% 6. error ids carry the CALLER, so messages stay attributable ------------
s = FUEL;  s.caller = 'refresh_duals_mee';
expectThrow(@() certified_guard(ok(0.5), s), 'refresh_duals_mee:massDegraded', 'caller id');
fprintf('  error ids     : prefixed by spec.caller\n');

fprintf('\ntest_certified_guard: ALL PASS\n');

function assertThrows(fn, wantId, what)
% ASSERTTHROWS  Assert fn errors with exactly wantId.
% INPUTS:  fn - zero-arg handle; wantId - expected error identifier [char];
%          what - short description for the failure message [char]
% OUTPUTS: none
try
    fn();
catch ME
    assert(strcmp(ME.identifier, wantId), ...
        'test_certified_guard [%s]: expected id %s, got %s (%s)', ...
        what, wantId, ME.identifier, ME.message);
    return
end
error('test_certified_guard:noThrow', ...
    'test_certified_guard [%s]: expected %s, but nothing was thrown', what, wantId);
end
