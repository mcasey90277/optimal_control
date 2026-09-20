function ok = test_family_costate_identity()
% TEST_FAMILY_COSTATE_IDENTITY  A root belongs to a family when the family's
% arc passes through it -- the same arrival phase, the same flight time AND
% THE SAME COSTATES. Flight time alone cannot say so where two families cross
% in the (phase, t_f) plane, which is exactly where a new family is found.
%
% Two synthetic families, A and C, whose arcs cross at sA = 0.30 with the SAME
% flight time, 19 d, and different initial costates (along e1 and along e2):
%   a root at (0.30, 19 d) with A's costates   -> family A
%   the same phase and time with C's costates  -> family C
%   the same phase and time, costates along e3 -> UNATTACHED: a third branch
%       through that point, which the driver must anchor and walk, and which
%       the flight-time rule filed under whichever arc it met first
%   no costates given (a legacy caller)        -> the flight-time rule, as before
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
arcDir = fullfile(tempdir, sprintf('famz_%s', char(java.util.UUID.randomUUID())));  mkdir(arcDir);
cleaner = onCleanup(@() rmdir(arcDir, 's'));
tStar = 86400*4;  n = 12;  nExtra = 1;  rho = 0.05;             % homogeneous chart: lam0_hom = rho * lam0
e = eye(7);
q = linspace(0.10, 0.50, 41);
writeArc(arcDir, 'arrival_arc_alpha_up_long.mat', q, (16 + 10*q)/4, rho, e(:, 1), n, nExtra);   % t_f = 16 + 10 sA
writeArc(arcDir, 'arrival_arc_gamma_up_long.mat', q, (22 - 10*q)/4, rho, e(:, 2), n, nExtra);   % t_f = 22 - 10 sA: crosses at 0.30, 19 d
sA = 0.05 + (0:9)/10;
S = struct('arcs', {{'arrival_arc_alpha_up_long.mat', 'arrival_arc_gamma_up_long.mat'}}, 'sA', sA, ...
           'TF', nan(1, 10), 'cand', {cell(1, 10)}, 'problem', struct('tStar', tStar));
F = family_map(S, struct('arcDir', arcDir));

z = @(dir) [dir; 19/4];
ok = chk(ok, F.attach(0.30, 19, z(3.7*e(:, 1))) == 1, 'A''s costates at the crossing -> family A (any positive scale of the costate)');
ok = chk(ok, F.attach(0.30, 19, z(0.2*e(:, 2))) == 2, 'C''s costates at the same phase and time -> family C');
[fam, gap, cosGap] = F.attach(0.30, 19, z(e(:, 3)));
ok = chk(ok, fam == 0 && gap < 1e-9 && cosGap > 0.5, sprintf('a THIRD root through the same (phase, t_f): unattached (t_f gap %.1e d, 1 - cos %.2f)', gap, cosGap));
ok = chk(ok, F.attach(0.30, 19) >= 1, 'no costates given: the flight-time rule, as before');
ok = chk(ok, F.attach(0.20, 18, z(e(:, 1))) == 1 && F.attach(0.20, 18, z(e(:, 2))) == 0, 'away from the crossing: A''s arc with A''s costates, and not with C''s');
lam = e(:, 1) + 5e-3*e(:, 4);                                    % 1 - cos ~ 1.2e-5: a re-polish of the same root
ok = chk(ok, F.attach(0.20, 18, z(lam)) == 1, 'costates within 1 - cos 1e-4 are the same root');
if ok, fprintf('test_family_costate_identity: ALL PASS\n'); else, fprintf('test_family_costate_identity: FAIL\n'); end
end

function writeArc(arcDir, name, q, tfNd, rho, lamDir, n, nExtra)
% WRITEARC  A synthetic arclength_arrival output whose points carry an initial
% costate (rows 1:7, homogeneous chart).  INPUTS: as named.  OUTPUTS: none.
p = cell(1, numel(q));
for k = 1:numel(q)
    v = zeros(n, 1);  v(1:7) = rho*lamDir;  v(n - nExtra) = tfNd(k);  v(n) = rho;
    p{k} = v;
end
A = struct('q', q, 'p', {p}, 'folds', struct([]), 'stop', 'nStep', 'crossings', struct([]), ...
           'anc', struct('n', n, 'nExtra', nExtra)); %#ok<NASGU>
save(fullfile(arcDir, name), 'A');
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
