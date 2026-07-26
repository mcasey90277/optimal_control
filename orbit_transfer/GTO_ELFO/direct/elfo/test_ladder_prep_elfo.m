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
% The 20 mN request must never be served by a NOMINAL seed. This used to assert
% isempty(sfB) -- "nothing is returned" as a proxy for "nominal seeds are
% rejected". That proxy expired when the 20 mN pilot succeeded and banked
% energy_elfo_f1200_T20mN.mat: the filter now correctly RETURNS a seed, and the
% old assertion failed on the campaign's own progress rather than on a defect.
% Assert the property that was actually meant.
if ~isempty(sfB)
    B = load(sfB);
    assert(isfield(B,'fp') && abs(B.fp.thrustN - 0.020) < 1e-12, ...
        ['a 20 mN request was served by a seed whose fingerprint is not 20 mN ' ...
         '(%s) -- the thrust filter is not doing its job'], sfB);
    assert(~strcmp(sfB, fullfile(here,'results','energy_elfo_f1200.mat')), ...
        'the nominal legacy seed was returned for a 20 mN request');
end
fprintf('test_ladder_prep_elfo: ALL PASS\n');
