function ok = test_phase_lists()
% TEST_PHASE_LISTS  The phase grid is two LISTS, not a lattice: the sheet
% indexes columns by an explicit arrival list, the catalog file takes an
% explicit departure list, the ribs walk explicit targets, and the rib
% validator checks points against the list.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

% ---- rib_targets: nearest first, unwrapped, the spine excluded ----------
sD = [0 0.1 0.35 0.6 0.95];
[t, p] = rib_targets(sD, 0, -1);
ok = chk(ok, isequal(round(t, 12), round([-0.05 -0.4 -0.65 -0.9], 12)) && isequal(round(p, 12), round([0.95 0.6 0.35 0.1], 12)), ...
         sprintf('targets downward from 0 over a non-uniform list: %s', mat2str(t, 3)));
[t, p] = rib_targets(sD, 0, +1);
ok = chk(ok, isequal(round(t, 12), round([0.1 0.35 0.6 0.95], 12)) && isequal(round(p, 12), round([0.1 0.35 0.6 0.95], 12)), ...
         'targets upward keep the list order');
[t, p] = rib_targets(sD, 0.35, -1);
ok = chk(ok, isequal(round(p, 12), round([0.1 0 0.95 0.6], 12)) && all(t < 0) && issorted(-t), ...
         'a spine inside the list: the other four phases, nearest first, wrapping past 0');
[t, ~] = rib_targets((0:23)/24, 0, -1);
ok = chk(ok, numel(t) == 23 && max(abs(t - (-(1:23)/24))) < 1e-12, 'the 24-lattice reproduces the old -k/24 targets');

% ---- sheet_from_arcs on an explicit, non-uniform arrival list -----------
sA = [0.05 0.2 0.5 0.9];
certFn = @(pp, s) struct('ok', true, 'reason', 'certified', 'z', [pp(1:7); pp(8)], 'tfDays', pp(8), 'sA', s);
cr = @(level, tf) struct('level', level, 'q', level, 'p', [zeros(7,1); tf], 'converged', true, 'normR', 1e-12, 'afterIndex', 1);
arc1 = struct('crossings', [cr(0.2, 20), cr(0.5, 21), cr(1.05, 19.5)]);      % 1.05 wraps onto 0.05
arc2 = struct('crossings', [cr(0.9, 22), cr(0.5, 20.5)]);
S = sheet_from_arcs({arc1, arc2}, struct('sA', sA, 'certFn', certFn, 'logFile', ''));
ok = chk(ok, isequal(S.sA, sA) && S.nA == 4, 'the sheet''s columns are the list, in its order');
ok = chk(ok, isequal(S.TF, [19.5 20 20.5 22]), sprintf('candidates land in their columns, fastest wins: %s', mat2str(S.TF)));
ok = chk(ok, numel(S.cand{3}) == 2, 'two distinct roots at 0.5 are both kept');
bad = struct('crossings', cr(0.3, 18));
ok = chk(ok, refuses(@() sheet_from_arcs({bad}, struct('sA', sA, 'certFn', certFn, 'logFile', ''))), ...
         'a crossing at a level off the list is refused');
S12 = sheet_from_arcs({}, struct('sA0', 0.0754, 'nA', 12, 'certFn', certFn, 'logFile', ''));
ok = chk(ok, S12.nA == 12 && abs(S12.sA(2) - (0.0754 + 1/12)) < 1e-12, 'without a list the lattice is built as before');

% ---- rib_validate against a departure list --------------------------------
tmp = tempname;  mkdir(tmp);  f = fullfile(tmp, 'rib.mat');
pts = struct('ok', {true, true, true}, 'sD', {0.95, 0.6, 0.35}, 'z', {ones(8,1), ones(8,1), ones(8,1)}, 'tfDays', {19, 19.5, 20});
R = struct('j', 2, 'sA', 0.2, 'pts', pts, 'stop', 'complete', 'nSolve', 3);
problem = struct('thrustN', 0.07, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1, 'sD', 0);
save(f, 'R', 'problem');
[g, why] = rib_validate(f, struct('nPts', 4, 'col', 2, 'sA', 0.2, 'sD', sD, 'problem', problem));
ok = chk(ok, g, sprintf('a rib on the list validates (%s)', why));
[g, why] = rib_validate(f, struct('nPts', 4, 'col', 2, 'sA', 0.2, 'sD', [0 0.1 0.5 0.9], 'problem', problem));
ok = chk(ok, ~g && contains(why, 'off the'), sprintf('a rib off the list is refused: %s', why));
rmdir(tmp, 's');

if ok, fprintf('TEST_PHASE_LISTS: ALL PASS\n'); else, fprintf('TEST_PHASE_LISTS: FAIL\n'); end
end

function r = refuses(fh)
% REFUSES  True when fh throws.  INPUTS: fh.  OUTPUTS: r.
try, fh();  r = false; catch, r = true; end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
