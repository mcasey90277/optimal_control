function res = run_ifs_1p12()
% RUN_IFS_1P12  Rung-2 / make-or-break gate: full 1.12x min-fuel IFS solve.
%
% INPUTS:  none
% OUTPUTS: res - struct with out, cert, meta
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
seed = fullfile(here,'..','sundman_minfuel','results','minfuel','legacy_ms_f1120.mat');
[Z0, prob, meta] = ifs_seed(seed, struct('mode','full'));
fprintf('RUNG 2 (GATE) full 1.12x: k=%d seedRes=%.3e\n', prob.k, meta.seedResNorm);
out  = ifs_solve(Z0, prob, struct('tolR',1e-8,'maxIter',400));
save(fullfile(here,'ifs_1p12_results.mat'), 'out');          % persist the solve first
try
    cert = ifs_certify(out.Z, prob, meta);
    res  = struct('out',out,'cert',cert,'meta',meta);
catch ME
    fprintf('WARNING: ifs_certify errored (%s); saving solve without certificate.\n', ME.message);
    cert = [];
    res  = struct('out',out,'cert',[],'meta',meta);
end
save(fullfile(here,'ifs_1p12_results.mat'), 'res');
fprintf('\n=== RUNG 2 SUMMARY (DONE) ===\nk=%d seed||R||=%.3e -> ||R||=%.3e success=%d iters=%d', ...
        prob.k, out.seedResNorm, out.resNorm, out.success, out.iterations);
if ~isempty(cert)
    fprintf(' maxSwitchMove=%.3e cert.ok=%d\n%s\n', cert.switchMoveFromSeed, cert.ok, cert.text);
else
    fprintf(' (no certificate)\n');
end
end
