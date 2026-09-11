function ok = test_run_dro_tulip()
%% Purpose:
%
%   Tests run_dro_tulip -- the single-transfer front door. One call with a
%   departure phase and an arrival phase must solve, certify and report ONE
%   DRO -> tulip minimum-time transfer, reproducing the certified library
%   values: (0, 0.0754) -> 17.7976 d and (11/12, 0.0754) -> 17.8775 d.
%   A phase pair reached only by a walk must either come back certified or
%   come back with a named reason -- never a bare number.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

T = run_dro_tulip(0, 0.0754, struct('quiet', true));
ok = chk(ok, T.ok, sprintf('the anchor pair solves and certifies (%s)', T.reason));
ok = chk(ok, abs(T.tfDays - 17.7976) < 2e-3, sprintf('t_f = %.4f d (library 17.7976)', T.tfDays));
ok = chk(ok, abs(T.dvKms - 0.7485) < 1e-3 && abs(T.propellantKg - 12.20) < 0.02, ...
         sprintf('dV = %.4f km/s, fuel = %.2f kg', T.dvKms, T.propellantKg));
ok = chk(ok, T.conj == 1 && T.g.dimS == 1 && T.dz <= 1e-6 && T.flyKm < 1, ...
         sprintf('gates: conj %d, dim S %d, |dz| %.1e, flown %.3f km', T.conj, T.g.dimS, T.dz, T.flyKm));
ok = chk(ok, strcmp(T.source, 'library'), sprintf('a library pair is served from the library (source: %s)', T.source));

T2 = run_dro_tulip(11/12, 0.0754, struct('quiet', true));
ok = chk(ok, T2.ok && abs(T2.tfDays - 17.8775) < 2e-3, ...
         sprintf('departure 11/12: t_f = %.4f d (library 17.8775)', T2.tfDays));

% a pair off the library must be WALKED to, and the walk reported
T3 = run_dro_tulip(0, 0.0754 + 1/24, struct('quiet', true, 'wallSec', 300));
ok = chk(ok, ~strcmp(T3.source, 'library') && ~isempty(T3.reason), ...
         sprintf('an off-library pair is walked: source %s, %s, t_f = %.4f d', T3.source, T3.reason, T3.tfDays));
ok = chk(ok, ~T3.ok || (T3.conj == 1 && T3.dz <= 1e-6), ...
         'a walked pair is reported certified only with its gates passed');

if ok, fprintf('TEST_RUN_DRO_TULIP: ALL PASS\n'); else, fprintf('TEST_RUN_DRO_TULIP: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
