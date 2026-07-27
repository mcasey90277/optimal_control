function m = mesh_diagnostic(solFile, outPng)
% MESH_DIAGNOSTIC  Show what the collocation grid actually is, and test it.
%
% The mesh in this campaign is INHERITED from the energy backbone, not chosen
% (see ../README.md and doc/gto_tulip_guide.pdf sec 3.4). Nothing in the
% pipeline reports what it is, so "is the grid adequate?" has historically been
% unanswerable by inspection. This makes it answerable.
%
% Three panels, each testing a different property:
%
%   (a) TIME STEP vs TIME.  dt_k across the transfer, log scale. Shows the
%       perigee/apogee contrast directly -- the whole point of the Sundman
%       clock. A flat curve would mean the regularization is doing nothing.
%
%   (b) KAPPA vs EARTH DISTANCE, against the predicted power law. Note this
%       plots dt/dsigma, NOT dt: the sigma grid in this campaign is strongly
%       NON-uniform (measured ratio 1.6e8 between the largest and smallest
%       dsigma, identical in every backbone), so plotting dt alone conflates the
%       Sundman law with the mesh spacing and yields a meaningless slope. The
%       first version of this function made exactly that mistake and reported a
%       spurious "mesh is not following r1^p". Dividing by dsigma isolates
%       kappa, and THEN the slope must be p (nominally 1.5).
%
%   (c) THE SIGMA GRID ITSELF, as dsigma vs node index. This is the panel that
%       tells you the grid was never designed: it is inherited unchanged from
%       the no-resample map of the original time-mesh seed, and is bit-identical
%       across every energy backbone regardless of t_f. Nodes-per-revolution is
%       reported numerically alongside.
%
% WHAT THIS DOES NOT DO. It characterises ONE mesh; it cannot tell you whether
% that mesh is converged. Convergence needs the same problem solved on two or
% three densities and compared (observed order of accuracy; do m_f and the
% switch times stop moving). That study is specified in
% docs/superpowers/plans/2026-07-25-mesh-convergence-study.md and has not been
% run. Do not read a healthy plot here as evidence of a converged solution.
%
% INPUTS:
%   solFile - solution .mat path or struct: needs out.X ([r;v;m;t] in rows 1-8)
%             and sigma; a 9th cScale row is ignored [char|struct]
%   outPng  - (optional) write the figure here [char]; omit for on-screen only
%
% OUTPUTS:
%   m - struct of measurements:
%       .N .revs .nodesPerRev_median .nodesPerRev_min
%       .dt_min .dt_max .dt_ratio      time-step range and contrast
%       .pFit                          fitted exponent from panel (b)
%       .pNominal .pResidual           nominal p and |pFit - pNominal|
%
% REFERENCES:
%   [1] ../README.md "The mesh is inherited, not chosen".
%   [2] ../../doc/gto_tulip_guide.tex sec 3.4 (the mathematics).

if nargin < 2, outPng = ''; end
if isstruct(solFile), S = solFile; else, S = load(solFile); end
if isfield(S, 'out'), X = S.out.X; else, X = S.X; end
sigma = S.sigma(:);
p     = cr3bp_lt_params(0.025, 15, 2100);          % constants only: muStar, tStar
pNom  = 1.5;                                        % cfg.pSund

t   = X(8,:);                                       % carried time state
dt  = diff(t);                                      % per-interval physical step
tm  = 0.5*(t(1:end-1) + t(2:end));                  % interval midpoints
r1  = vecnorm(X(1:3,:) - [-p.muStar;0;0], 2, 1);    % Earth distance
r1m = 0.5*(r1(1:end-1) + r1(2:end));

% revolutions: unwrapped polar angle about the Earth
th   = unwrap(atan2(X(2,:), X(1,:) + p.muStar));
revs = (max(th) - min(th)) / (2*pi);

% nodes per revolution, per revolution
edges = linspace(min(th), max(th), max(floor(revs),1) + 1);
nPer  = histcounts(th, edges);

% panel (b): the Sundman law lives in dt/dsigma, not dt. dsigma is strongly
% non-uniform here, so it must be divided out before fitting the exponent.
dsg   = diff(sigma(:)).';
kap   = dt ./ dsg;                                  % proportional to kappa
good  = kap > 0 & r1m > 0 & isfinite(kap);
cf    = polyfit(log(r1m(good)), log(kap(good)), 1);
pFit  = cf(1);
dsRatio = max(dsg)/min(dsg);

m = struct('N', numel(sigma), 'revs', revs, ...
           'nodesPerRev_median', median(nPer), 'nodesPerRev_min', min(nPer), ...
           'dt_min', min(dt), 'dt_max', max(dt), 'dt_ratio', max(dt)/min(dt), ...
           'pFit', pFit, 'pNominal', pNom, 'pResidual', abs(pFit - pNom), ...
           'dsigma_ratio', dsRatio, 'dsigma_cv', std(dsg)/mean(dsg));

fprintf('\n=== mesh diagnostic ===\n');
fprintf('  nodes                 : %d\n', m.N);
fprintf('  revolutions           : %.1f\n', m.revs);
fprintf('  nodes/rev  median/min : %.0f / %.0f\n', m.nodesPerRev_median, m.nodesPerRev_min);
fprintf('  dt  min/max (ND)      : %.3e / %.3e   (ratio %.1f)\n', m.dt_min, m.dt_max, m.dt_ratio);
fprintf('  dsigma  ratio / CV    : %.3g / %.3f   (uniform grid => 1 / 0)\n', ...
        m.dsigma_ratio, m.dsigma_cv);
fprintf('  Sundman exponent p    : fitted %.3f vs nominal %.3f  (residual %.3f)\n', ...
        m.pFit, m.pNominal, m.pResidual);
if m.pResidual > 0.15
    fprintf('  ^^ fitted p is far from nominal -- the mesh is NOT following r1^p.\n');
end
fprintf('  NOTE: this characterises ONE mesh. It is not a convergence test.\n');

fig = figure('Color','w','Position',[80 80 1150 380]);

subplot(1,3,1);
semilogy(tm*p.tStar/86400, dt, '-', 'LineWidth', 0.9); grid on;
xlabel('time [days]'); ylabel('\Deltat per interval [ND]');
title(sprintf('(a) step vs time  (ratio %.0f\\times)', m.dt_ratio));

subplot(1,3,2);
loglog(r1m, kap, '.', 'MarkerSize', 3); hold on; grid on;
rr = linspace(min(r1m), max(r1m), 50);
loglog(rr, exp(cf(2))*rr.^pNom, 'r--', 'LineWidth', 1.4);
xlabel('r_1  (Earth distance, ND)'); ylabel('\Deltat/\Delta\sigma  (\propto \kappa)');
title(sprintf('(b) fitted p = %.3f  (nominal %.2f)', m.pFit, pNom));
legend({'nodes', sprintf('r_1^{%.1f}', pNom)}, 'Location','northwest');

subplot(1,3,3);
semilogy(dsg, '-', 'LineWidth', 0.8); grid on;
xlabel('node index'); ylabel('\Delta\sigma');
title(sprintf('(c) the grid itself: \\Delta\\sigma ratio %.1e\nnodes/rev median %.0f, min %.0f', ...
      m.dsigma_ratio, m.nodesPerRev_median, m.nodesPerRev_min));

if ~isempty(outPng)
    exportgraphics(fig, outPng, 'Resolution', 150);
    fprintf('  figure -> %s\n', outPng);
end
end
