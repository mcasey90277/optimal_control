function ok = test_alternative_rows()
% TEST_ALTERNATIVE_ROWS  make_alternative builds one row from a certify_root
% result (stamped or legacy), refuses what is below the floor or overridden;
% dedup_alternatives removes a primary's twin and repeated roots; same_root
% is symmetric; the delegates agree with the shared units.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
z = [(1:7).'; 4.1];  Y = rand(14, 24);
C = struct('ok', false, 'z', z, 'Y', Y, 'flyKm', 0.2, 'flyVms', 0.01, 'reason', 'H6 not established: x', ...
           'stage', 6, 'conjFound', false, 'hypAfterConj', '', 'overridden', false, 'status', 3, ...
           'status_reason', 'necessary only: H6 not established: x', 'conj', 1, 'conjVerdict', 'PASS', ...
           'g', struct('minLamV', 0.1, 'minQmt', 0.2, 'dimS', 1), 'h6Margin', 0.9, 'liftMargin', 30);
A = make_alternative(C, 0.25, 0.5, 7, 13, 'test');
ok = chk(ok, isstruct(A) && A.status == 3 && isequal(A.z8, z) && isequal(A.junctions, Y) && A.tf_nd == 4.1 && ~A.inferred ...
         && A.iD == 7 && A.dimS == 1 && strcmp(A.source, 'test'), 'a stamped result becomes one row with its numbers');
L = rmfield(C, {'stage', 'conjFound', 'hypAfterConj', 'status', 'status_reason', 'overridden', 'conjVerdict'});
L.reason = 'conjugate test verdict 0';
A2 = make_alternative(L, 0.25, 0.5, 7, 13, 'legacy');
ok = chk(ok, A2.status == 3 && A2.inferred, 'a legacy result is classified on the way in (inferred)');
Cb = C;  Cb.flyKm = 500;
ok = chk(ok, isempty(make_alternative(Cb, 0.25, 0.5, 7, 13, 'x')), 'below the floor -> no row');
Co = C;  Co.overridden = true;
ok = chk(ok, isempty(make_alternative(Co, 0.25, 0.5, 7, 13, 'x')), 'an overridden (test-seam) result -> no row');
Cc = C;  Cc.ok = true;  Cc.status = 4;
A4 = make_alternative(Cc, 0.25, 0.5, 7, 13, 'x');
ok = chk(ok, A4.status == 4, 'a certified slower root is a status-4 row');

% dedup: a primary's twin, a repeated root, a distinct root kept
sheet = struct('sD_frac', [0 0.25], 'sA_frac', [0.5 0.75], 'has_solution', [false false; true false], ...
               'entry_index', [0 0; 1 0], 'z8', z);
Ad = [make_alternative(C, 0.25, 0.5, 2, 1, 'twin of primary'), ...
      make_alternative(setz(C, z + [1; zeros(7, 1)]), 0.25, 0.5, 2, 1, 'distinct'), ...
      make_alternative(setz(C, z + [1; zeros(7, 1)]), 0.25, 0.5, 2, 1, 'repeat')];
Ad = dedup_alternatives(Ad, sheet);
ok = chk(ok, numel(Ad) == 1 && strcmp(Ad.source, 'distinct'), 'dedup drops the primary''s twin and the repeat');

% same_root and the delegates
ok = chk(ok, same_root(z, z + 1e-9) && same_root(z + 1e-9, z) && ~same_root(z, z + [0.01; zeros(7, 1)]), 'same_root: symmetric, 1e-6');
H = run_phase_torus('localfunctions');
ok = chk(ok, H.sameRoot(z, z + 1e-9) == same_root(z, z + 1e-9), 'run_phase_torus''s sameRoot delegates');
Lc = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = Lc.(char(fieldnames(Lc)));
Hf = fill_holes_direct('localfunctions');
ok = chk(ok, isequal(Hf.physicsFromCatalog(c, c.sheets(1), 0), catalog_setup_request(c, 0)), 'fill_holes_direct''s physicsFromCatalog delegates');

if ok, fprintf('test_alternative_rows: ALL PASS\n'); else, fprintf('test_alternative_rows: FAIL\n'); end
end

function C = setz(C, z)
% SETZ  Replace a result's z.  INPUTS: C; z.  OUTPUTS: C.
C.z = z;
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
