function ok = test_polish_costate_guess()
% TEST_POLISH_COSTATE_GUESS  Polishing a costate guess into a root, on one
% entry of the library of record: from the root itself (nothing to do), from a perturbed
% root (must come back to the SAME root), and with too few iterations (must
% say it did not converge and return no root -- never a best iterate).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
% a certified entry of the library of record, with its EXACT phases (the anchor in
% dro_tulip_library carries a rounded arrival phase, 49 km off its own endpoint)
catMat = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
assert(isfile(catMat), 'the library of record is missing: %s', catMat);
L = load(catMat);  fn = fieldnames(L);  sh = L.(fn{1}).sheets(1);
[iD, iA] = find(sh.has_solution(:, :, 1), 1);
E = struct('z', sh.z8(:, sh.entry_index(iD, iA)), 'sD', sh.sD_frac(iD), 'sA', sh.sA_frac(iA));
[B, ~] = arclength_arrival('setup', struct('physicsOnly', true));
rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'muStar', B.mu);

[z, info] = polish_costate_guess(E.z(:), rv0(1:6), rvf(1:6), phys);
ok = chk(ok, info.converged && info.iters <= 3 && max(abs(z - E.z(:))) < 1e-8, sprintf('from the root: converged in %d iteration(s), unmoved', info.iters));
ok = chk(ok, info.guessMissKm < 1, sprintf('the guess''s own flown miss is reported (%.3f km)', info.guessMissKm));

zG = E.z(:);  zG(1:7) = zG(1:7)*(1 + 2e-4);
[z, info] = polish_costate_guess(zG, rv0(1:6), rvf(1:6), phys);
ok = chk(ok, info.converged && max(abs(z - E.z(:))) < 1e-7, sprintf('from a 2e-4 perturbation: back to the same root in %d iterations (guess flew %.0f km off)', info.iters, info.guessMissKm));

zBad = E.z(:);  zBad(1:7) = zBad(1:7)*(1 + 1e-2);
[z, info] = polish_costate_guess(zBad, rv0(1:6), rvf(1:6), phys, struct('maxIter', 1, 'polishMax', 0));
ok = chk(ok, ~info.converged && all(isnan(z)) && isfinite(info.normR), sprintf('a 1%% perturbation, ONE search iteration and no Newton steps after it: not converged (%d), no root returned (%d), the residual reported (%.1e)', info.converged, all(isnan(z)), info.normR));

ok = chk(ok, throws(@() polish_costate_guess([zG(1:7); NaN], rv0(1:6), rvf(1:6), phys)), 'a non-finite guess is refused before anything is flown');
if ok, fprintf('test_polish_costate_guess: ALL PASS\n'); else, fprintf('test_polish_costate_guess: FAIL\n'); end
end

function tf = throws(f)
% THROWS  Does calling f throw?  INPUTS: f (handle).  OUTPUTS: tf.
tf = false;
try, f(); catch, tf = true; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
