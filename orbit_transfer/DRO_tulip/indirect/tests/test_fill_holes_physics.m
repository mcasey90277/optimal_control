function ok = test_fill_holes_physics()
% TEST_FILL_HOLES_PHYSICS  The filler reconstructs the physics from the
% catalog WITHOUT polishing an anchor (FINDINGS 78): its setup request asks
% for the closures only, at the catalog's spine departure phase, with the
% catalog's engine and orbits.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
catMat = fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});  sh = cat_.sheets(1);
H = fill_holes_direct('localfunctions');
so = H.physicsFromCatalog(cat_, sh, sh.sD_frac(1));
ok = chk(ok, isfield(so, 'physicsOnly') && so.physicsOnly, 'closures only: no anchor is polished');
ok = chk(ok, ~isfield(so, 'anchorMat'), 'no anchor file is named');
ok = chk(ok, so.sD == sh.sD_frac(1) && so.thrustN == cat_.rungs_N(1) && so.NpTulip == sh.Np && so.pmTulip == sh.pm, ...
         'the catalog''s spine phase, engine and orbits');
% ---- the rib builder: the same rule, from the SHEET's identity ------------
Hr = build_ribs('localfunctions');
P = struct('thrustN', 0.05, 'ispS', 1200, 'm0kg', 300, 'tauDRO', 2, 'NpTulip', 8, 'pmTulip', +1, 'sD', 0.2);
so = Hr.setupFromSheet(struct('problem', P), struct());
ok = chk(ok, isfield(so, 'physicsOnly') && so.physicsOnly && ~isfield(so, 'anchorMat'), 'build_ribs: closures only, no anchor');
ok = chk(ok, so.sD == 0.2 && so.thrustN == 0.05 && so.NpTulip == 8 && so.pmTulip == 1 && so.m0kg == 300, ...
         'build_ribs: the sheet''s engine, orbits, branch and departure phase');
if ok, fprintf('test_fill_holes_physics: ALL PASS\n'); else, fprintf('test_fill_holes_physics: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
