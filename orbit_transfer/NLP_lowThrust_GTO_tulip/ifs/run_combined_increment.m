function res = run_combined_increment(runs)
% RUN_COMBINED_INCREMENT  Rung 0 + Rung 1: scaled truncated-SVD Newton (ifs_solve2)
% on the IFS rendezvous gates, in both the sigmoid and direct-tau
% parameterizations, compared against the old lsqnonlin crawl (1.96 -> 0.023,
% stalled).
%
% NB: the planned minted "3-switch easy gate" is NOT cleanly available -- the
% 3-switch LOCAL optimum is not a fixed point of the Sundman eps=0 re-solve (it
% slides toward the many-switch global basin), and every compatible legacy
% solution is many-switch (1.12x=10sw is in fact the SMALLEST). So the gates are
% the real cases: 1.12x (hard terminal cluster, k=10) is both the smallest
% system and the actual open target; 1.25x (k=50) is the GPT-recommended clean,
% well-separated benchmark.
%
% INPUTS:
%   runs - optional Nx4 cell {label, seedFile, tauParam, maxIter}; default runs
%          the 1.12x gate in both parameterizations
% OUTPUTS: res - struct array of per-run results (label, k, seedRes, resNorm,
%          success, iters, flag, rankEq, condEq, cert)
%
% REFERENCES: ifs/PLAN_OF_ATTACK.md (section 3b).
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
minf = @(f) fullfile(here, '..', 'sundman_minfuel', 'results', 'minfuel', f);

if nargin < 1 || isempty(runs)
    gate12 = minf('legacy_ms_f1120.mat');
    runs = { '1.12x sigmoid', gate12, 'sigmoid', 200; ...
             '1.12x direct',  gate12, 'direct',  200 };
end

res = struct('label',{},'k',{},'seedRes',{},'resNorm',{},'success',{}, ...
             'iters',{},'flag',{},'rankEq',{},'condEq',{},'cert',{});
for i = 1:size(runs,1)
    label = runs{i,1};  seed = runs{i,2};  tp = runs{i,3};  mi = runs{i,4};
    fprintf('\n================ %s ================\n', label);
    [Z0, prob, meta] = ifs_seed(seed, struct('mode','full','tauParam',tp));
    out = ifs_solve2(Z0, prob, struct('tolR',1e-8,'maxIter',mi,'relTrunc',1e-10));
    cert = [];
    try
        cert = ifs_certify(out.Z, prob, meta);
    catch ME
        fprintf('  certify errored: %s\n', ME.message);
    end
    res(end+1) = struct('label',label,'k',prob.k,'seedRes',out.seedResNorm, ...
        'resNorm',out.resNorm,'success',out.success,'iters',out.iterations, ...
        'flag',out.flag,'rankEq',out.rankEq,'condEq',out.condEq, ...
        'cert',cert); %#ok<AGROW>
end

fprintf('\n\n=================== SUMMARY ===================\n');
fprintf('%-20s k  seed||R||  final||R||  succ iters flag rankEq  condEq   certOK\n', 'label');
for i = 1:numel(res)
    r = res(i);
    certok = ~isempty(r.cert) && isfield(r.cert,'ok') && r.cert.ok;
    fprintf('%-20s %2d  %.2e  %.2e   %d  %4d  %3d  %4d  %.2e   %d\n', ...
        r.label, r.k, r.seedRes, r.resNorm, r.success, r.iters, r.flag, ...
        r.rankEq, r.condEq, certok);
end
save(fullfile(here, 'combined_increment_results.mat'), 'res');
end
