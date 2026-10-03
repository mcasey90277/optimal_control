function ok = test_rib_as_catalog()
% TEST_RIB_AS_CATALOG  A departure rib (one arrival phase, many departure
% phases) wrapped as a one-COLUMN phase catalog: the spine crossing it was
% hung from plus every certified rib point become entries with the same
% eight numbers, sorted by departure phase; points off the departure lattice
% (bisection steps) are dropped when a lattice is given; a rib that is not
% one arrival phase, or not the sheet's column, is refused; the problem
% identity is the reference catalog's; the interpolator and the scorer's
% one-column rule read it.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  fn = fieldnames(L);  ref = L.(fn{1});

% a sheet at s_D = 0 with two columns; the rib hangs off column 2 (s_A 0.3)
S = struct('sA', [0.1 0.3], 'TF', [17 18], 'Z8', [repmat((1:7).', 1, 2); 4.0 4.2], 'problem', struct('sD', 0));
S.Z8(1, :) = [1 2];
% walked in the -1 direction at 1/8: 7/8, 6/8, a bisection point at 0.8125 off the lattice, 5/8
pt = @(sD, z1, tf) struct('sD', sD, 'sA', 0.3, 'z', [z1; (2:7).'; tf], 'tfDays', NaN, 'ok', true);
R = struct('j', 2, 'sA', 0.3, 'pts', [pt(-1/8, 2.2, 4.3), pt(-2/8, 2.4, 4.4), pt(-0.1875, 2.3, 4.35), pt(-3/8, 2.6, 4.5)], ...
           'stop', 'done', 'nSolve', 9);

c = rib_as_catalog(R, S, ref, struct('nD', 8));
sh = c.sheets(1);
ok = chk(ok, isscalar(sh.sA_frac) && sh.sA_frac == 0.3 && isequal(size(sh.has_solution), [4 1]) && c.n_entries == 4, ...
         'one column: the spine crossing + the three lattice points');
ok = chk(ok, max(abs(sh.sD_frac - [0 5/8 6/8 7/8])) < 1e-12 && issorted(sh.sD_frac), 'departure phases wrapped to [0,1) and sorted');
ok = chk(ok, isequal(sh.z8(:, sh.entry_index(1)), S.Z8(:, 2)) && sh.z8(1, sh.entry_index(4)) == 2.2 && sh.tf_nd(3) == 4.4, ...
         'the same eight numbers: the spine from the sheet, the rest from the rib; t_f in ND from z8(8)');
ok = chk(ok, c.rib.nOffLattice == 1 && c.rib.nPts == 4 && strcmp(c.rib.stop, 'done'), 'the bisection point is dropped, and counted');
c2 = rib_as_catalog(R, S, ref);
ok = chk(ok, c2.n_entries == 5 && any(abs(c2.sheets(1).sD_frac - 0.8125) < 1e-12), 'with no lattice given, every certified point is kept');
ok = chk(ok, all(sh.family_index == 1) && isequal(c.families.labels, {'rib'}), 'one family: a rib is one continuation from one crossing');
ok = chk(ok, isequal(c.thruster, ref.thruster) && isequal(c.constants, ref.constants) && sh.Np == ref.sheets(1).Np, ...
         'the problem identity is the reference catalog''s');

% refusals
Rb = R;  Rb.pts(2).sA = 0.31;
ok = chk(ok, refuses(@() rib_as_catalog(Rb, S, ref)), 'refuses a rib whose points are not all at one arrival phase');
Rc = R;  Rc.j = 1;
ok = chk(ok, refuses(@() rib_as_catalog(Rc, S, ref)), 'refuses a rib that is not the sheet column it names');
Rd = R;  Rd.pts(1).z(3) = NaN;
ok = chk(ok, refuses(@() rib_as_catalog(Rd, S, ref)), 'refuses a non-finite point');

% the readers: the interpolator blends along departure; the scorer pins arrival
[z, info] = phase_catalog_interp(sh, 0.8125, 0.3);
ok = chk(ok, strcmp(info.tier, 'linear') && abs(z(1) - 2.3) < 1e-12, 'the interpolator reads it (one arrival line: the blend is along departure phase)');
ok = chk(ok, isscalar(sh.sA_frac) && ~isscalar(sh.sD_frac), 'score_interpolator''s one-column rule applies (arrival pinned)');

if ok, fprintf('test_rib_as_catalog: ALL PASS\n'); else, fprintf('test_rib_as_catalog: FAIL\n'); end
end

function r = refuses(f)
% REFUSES  True if f() throws.  INPUTS: f (handle).  OUTPUTS: r.
try, f();  r = false; catch, r = true; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
