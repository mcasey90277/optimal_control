function outFile = prep_refine_seed(seedFile, outFile)
% PREP_REFINE_SEED  Ensure a direct solution carries duals for refinement.
%
% The certified 1.15x .mat was saved before dual extraction existed and
% lacks out.lamDef and a factor field. This loads such a file, re-solves
% eps=0 warmTight from its own (X,U) to regenerate the KKT-dual costates
% (out.lamDef), stamps factor = round(tf/tfMin, 2), and writes the layout
% pmp_refine_indicator / sms_seed_duals require. A file that already carries
% out.lamDef is passed through with fields normalized (no re-solve).
%
% INPUTS:
%   seedFile - source solution .mat (certified layout: out, sigma, tauf0,
%              rv0, rvf, tf; out may lack lamDef/factor)
%   outFile  - destination .mat path
%
% OUTPUTS:
%   outFile - the written path [char] (echoes the input)
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md
%   [2] sundman_minfuel/casadi_minfuel_sundman.m (solver + dual extraction)

S   = load(seedFile);
p   = cr3bp_lt_params(0.025, 15, 2100);
out = S.out;
tf  = S.out.X(8, end);                       % carried terminal time = tf
tfMin = 6.290694;                            % campaign constant (ND)
factor = round(tf/tfMin, 2);

if ~isfield(out, 'lamDef') || isempty(out.lamDef)
    fprintf('prep_refine_seed: regenerating duals (eps=0 warmTight re-solve)...\n');
    out = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, ...
              p.muStar, S.out.X, S.out.U, S.tauf0, 1.5, 3000, 0, true);
    assert(out.success && out.maxDefect < 1e-6 && ~isempty(out.lamDef), ...
           'seed re-solve failed: success=%d defect=%.2e', out.success, out.maxDefect);
end

sigma = S.sigma;  tauf0 = S.tauf0;  rv0 = S.rv0;  rvf = S.rvf; %#ok<NASGU>
save(outFile, 'out', 'factor', 'tauf0', 'sigma', 'rv0', 'rvf');
fprintf('prep_refine_seed: wrote %s (factor=%.2f, switches=%d)\n', ...
        outFile, factor, out.switches);
end
