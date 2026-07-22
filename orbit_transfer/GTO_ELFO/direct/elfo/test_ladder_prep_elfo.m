% TEST_LADDER_PREP_ELFO  chain helper, seed fp filter, cBox rule (no solves).
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
cfg = minfuel_config();  p20 = cr3bp_lt_params(0.020, cfg.m0kg, cfg.ispS);
% (a) chain helper pass-through + fp
S = load(fullfile(here,'results','energy_elfo_f1200.mat'));
[S2, fp] = chain_rung_seed_elfo(S, p20, struct('note','test'));
assert(isequal(S2.X, S.X) && isequal(S2.U, S.U), 'pass-through must not touch X/U');
assert(fp.thrustN==0.020 && isfield(fp,'chainedFrom'), 'fp');
ok=false; p25 = cr3bp_lt_params(0.025, cfg.m0kg, cfg.ispS);
try, chain_rung_seed_elfo(S, p25, struct());
catch err, ok = strcmp(err.identifier,'chain_rung_seed_elfo:sameThrust'); end
assert(ok, 'same-thrust refusal');
% (b) seed fp filter: legacy seeds eligible under a 25 mN fp, skipped under 20 mN
fp25 = cr3bp_fingerprint(p25); fp20 = cr3bp_fingerprint(p20);
w = warning('off','all');
[sfA,~,~] = elfo_find_energy_seed(fullfile(here,'results'), S.X(8,end), 0.02, fp25);
[sfB,~,~] = elfo_find_energy_seed(fullfile(here,'results'), S.X(8,end), 0.02, fp20);
warning(w);
assert(~isempty(sfA), 'legacy seed eligible under nominal fp');
assert(isempty(sfB), 'legacy (nominal) seed must NOT satisfy a 20 mN request');
fprintf('test_ladder_prep_elfo: ALL PASS\n');
