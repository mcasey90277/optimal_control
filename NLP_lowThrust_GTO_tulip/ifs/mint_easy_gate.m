function outFile = mint_easy_gate(outFile)
% MINT_EASY_GATE  Build an ifs_seed-compatible few-switch rendezvous seed.
%
% The genuine 3-switch min-fuel local optimum (burn-then-coast, tf=1.15x) is
% stored only in the OLD 7-state fmincon layout (minfuel_from_energy_seed.mat:
% nlp.X[7x4001], no time state, no costates, no factor), which ifs_seed cannot
% read. This mints the modern layout: map the 7-state solution into Sundman
% 8-state (sundman_seed_map) and re-solve eps=0, warmTight
% (casadi_minfuel_sundman) to regenerate out.X[8]/out.lamDef and stamp factor,
% preserving the low-switch structure. Output layout matches what ifs_seed /
% sms_seed_duals require (out, factor, tauf0, sigma, rv0, rvf).
%
% INPUTS:
%   outFile - destination .mat [char, default ifs/seed_3sw_1p15.mat]
% OUTPUTS:
%   outFile - the written path [char]
%
% REFERENCES:
%   [1] sundman_minfuel/refine/prep_refine_seed.m (the dual-regen pattern).
%   [2] ifs/PLAN_OF_ATTACK.md (Rung 0 minted easy gate).

here = fileparts(mfilename('fullpath'));
addpath(here);  setup_paths();
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));      % CasADi (defensive)
if nargin < 1 || isempty(outFile)
    outFile = fullfile(here, 'seed_3sw_1p15.mat');
end

S     = load(fullfile(here, '..', 'sundman_minfuel', 'minfuel_from_energy_seed.mat'));
p     = cr3bp_lt_params(0.025, 15, 2100);
pSund = 1.5;
tf    = S.tf;                                           % this seed's own tf (1.15x)
rv0   = S.rv0;  rvf = S.rvf;

fprintf('mint_easy_gate: mapping 7-state 3-switch seed -> Sundman 8-state (tf=%.4f)\n', tf);
[sigma, X0, U0, tauf0] = sundman_seed_map(S.nlp.X, S.nlp.U, tf, S.sigma, ...
                                          pSund, p.muStar, rv0, rvf);

fprintf('mint_easy_gate: eps=0 warmTight re-solve to regenerate duals...\n');
out = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                             X0, U0, tauf0, pSund, 3000, 0, true);
assert(out.success && out.maxDefect < 1e-6 && ~isempty(out.lamDef), ...
       'mint re-solve failed: success=%d defect=%.2e', out.success, out.maxDefect);

factor = round(tf/6.290694, 2);                         % tfMin campaign constant (ND)
save(outFile, 'out', 'factor', 'tauf0', 'sigma', 'rv0', 'rvf');
fprintf('mint_easy_gate: wrote %s (factor=%.2f, switches=%d, defect=%.2e, mf=%.6f)\n', ...
        outFile, factor, out.switches, out.maxDefect, out.mf);
end
