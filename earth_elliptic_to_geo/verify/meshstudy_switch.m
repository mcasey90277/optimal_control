function rows = meshstudy_switch(srcMat, densities, opts)
% MESHSTUDY_SWITCH  Primal mesh-convergence study of a certified low-thrust
% rung's bang-bang switch structure. Warm-refines a converged solution to a
% sequence of finer UNIFORM node densities (chained, each refining the previous)
% and re-solves to eps=0 at each, then reports the switch count / duty / mass so
% their convergence (or non-convergence) can be read directly. This is the
% load-bearing test of P0 (process/P0_SWITCH_MESH_CONVERGENCE.md): it is purely
% primal, so it sidesteps the raw-dual anomaly that makes the PMP
% switching-function cross-check inconclusive at high eccentricity.
%
% METHOD (per density): interp_warmstart the previous converged eps=0 solution
% onto the finer uniform sigma grid, then run a SHORT warm eps-tail
% (default [0.01 0.003 0.001 0]) via homotopy_mee -- per-eps-step cached, so a
% MEX crash mid-density resumes near where it died -- with the deep-rung levers
% (liftDL, generous maxIter). Same-rung refinement keeps DeltaL/rev-count fixed,
% so the sigma-interp of beta does NOT phase-alias (unlike a cross-thrust chain,
% where warmstart_phase_beta is required); interp_warmstart is correct here.
%
% RESUME-SAFE across process death: each density caches to
% opts.resDir/meshconv_<npr>pr.mat; a density whose cache exists is loaded and
% skipped, so re-invoking after any crash continues the study.
%
% INPUTS:
%   srcMat    - path to a certified res-struct .mat (fields res.fuel.{X,U,dL,
%               m_f_kg,maxDefect}, res.sigma, res.tf, res.fp.{thrustN,m0kg,
%               ispS}) -- e.g. 'results/MEE_M2_0p2N.mat' [char]
%   densities - target nodes/rev to refine THROUGH, ascending [1xK]
%               (e.g. [16 24 40]); the source density is row 1 automatically
%   opts      - struct (optional): .resDir [tempdir/meshconv_<tag>],
%               .maxIter [4000], .tail [ [0.01 0.003 0.001 0] ],
%               .liftDL [true]
%
% OUTPUTS:
%   rows - 1x(K+1) struct array (source density first), each with fields
%          .npr .nodes .revs .nSw .duty .m_f_kg .defect .certified .swPhase
%          Also written to opts.resDir/meshconv_summary.mat.
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/process/P0_SWITCH_MESH_CONVERGENCE.md (finding).
%   [2] earth_elliptic_to_geo/verify/switch_structure.m (the primal metric).
%   [3] earth_elliptic_to_geo/core/{interp_warmstart,homotopy_mee}.m (reused).
if nargin < 3, opts = struct(); end
d = @(f,v) optdef(opts, f, v);

S0  = load(srcMat);
assert(isfield(S0,'res'), 'meshstudy_switch: %s must hold a res-struct', srcMat);
r0  = S0.res;
par = kepler_lt_params(r0.fp.thrustN, r0.fp.m0kg, r0.fp.ispS);
tfT = r0.tf;

tag     = sprintf('%gN', r0.fp.thrustN);
resDir  = d('resDir', fullfile(tempdir, ['meshconv_' tag]));
maxIter = d('maxIter', 4000);
tail    = d('tail', [0.01 0.003 0.001 0]);
liftDL  = d('liftDL', true);
if ~exist(resDir,'dir'), mkdir(resDir); end

% row 0: the source (certified) density
base = struct('X',r0.fuel.X,'U',r0.fuel.U,'dL',r0.fuel.dL,'sigma',r0.sigma(:), ...
    'npr',round((numel(r0.sigma)-1)/(r0.fuel.dL/(2*pi))), ...
    'm_f_kg',r0.fuel.m_f_kg,'maxDefect',r0.fuel.maxDefect);
ss0  = switch_structure(base.X, base.U, base.dL, base.sigma);
rows = struct('npr',base.npr,'nodes',numel(base.sigma),'revs',ss0.revs, ...
    'nSw',ss0.nSwitch,'duty',ss0.duty,'m_f_kg',base.m_f_kg, ...
    'defect',base.maxDefect,'certified',1,'swPhase',ss0.swPhase);

prev = base;
for npr = densities(:).'
    dcache = fullfile(resDir, sprintf('meshconv_%dpr.mat', npr));
    if exist(dcache,'file')
        L = load(dcache);  rows(end+1) = L.row;  prev = L.prev; %#ok<AGROW>
        fprintf('[cached] %d/rev: nodes=%d nSw=%d m_f=%.3f defect=%.2e cert=%d\n', ...
            npr, L.row.nodes, L.row.nSw, L.row.m_f_kg, L.row.defect, L.row.certified);
        continue;
    end
    N        = round(npr * prev.dL/(2*pi));
    sigmaDst = linspace(0,1,N+1).';
    fprintf('\n>>> refining %d/rev -> %d/rev : N=%d nodes\n', prev.npr, npr, N);
    W  = interp_warmstart(prev.X, prev.U, prev.dL, prev.sigma, sigmaDst);
    ho = struct('par',par,'x0',W.X(:,1),'tfTarget',tfT,'maxIter',maxIter, ...
        'resDir', fullfile(resDir, sprintf('steps_%dpr', npr)), ...
        'tag', sprintf('meshconv_%s_%dpr', tag, npr), 'printLevel', 0, ...
        'fp', struct('thrustN',r0.fp.thrustN,'npr',npr), 'xf',[1;0;0;0;0], ...
        'sched', tail, 'liftDL', liftDL, 'scaleNLP', false);
    [best,~] = homotopy_mee(sigmaDst, W.X, W.U, W.dL, ho);

    ss  = switch_structure(best.X, best.U, best.dL, sigmaDst);
    row = struct('npr',npr,'nodes',N+1,'revs',ss.revs,'nSw',ss.nSwitch, ...
        'duty',ss.duty,'m_f_kg',best.m_f_kg,'defect',best.maxDefect, ...
        'certified',best.certified,'swPhase',ss.swPhase);
    prev = struct('X',best.X,'U',best.U,'dL',best.dL,'sigma',sigmaDst,'npr',npr);
    save(dcache, 'row', 'prev');
    rows(end+1) = row; %#ok<AGROW>
    fprintf('[done]  %d/rev: nodes=%d nSw=%d sw/rev=%.3f m_f=%.3f defect=%.2e cert=%d\n', ...
        npr, row.nodes, row.nSw, row.nSw/row.revs, row.m_f_kg, row.defect, row.certified);
end

fprintf('\n===== %s N MESH-CONVERGENCE (primal) =====\n', tag);
fprintf('%6s %8s %8s %6s %8s %11s %10s %5s\n', ...
    'npr','nodes','revs','nSw','sw/rev','m_f[kg]','defect','cert');
for i = 1:numel(rows)
    R = rows(i);
    fprintf('%6d %8d %8.2f %6d %8.3f %11.3f %10.1e %5d\n', ...
        R.npr, R.nodes, R.revs, R.nSw, R.nSw/R.revs, R.m_f_kg, R.defect, R.certified);
end
save(fullfile(resDir,'meshconv_summary.mat'), 'rows');
end
