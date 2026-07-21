function chk = hamiltonian_const_check(matPath, plotPath)
% HAMILTONIAN_CONST_CHECK  Independent first-order optimality test: verify that
% the (time-domain) Hamiltonian is CONSTANT along a min-fuel MEE solution, via
% the time-costate.
%
% THE PHYSICS. The problem is solved in the L-domain (independent variable =
% true longitude L) with TIME carried as a state (X = [P;ex;ey;hx;hy;m;t]). The
% dynamics and cost are AUTONOMOUS in t -- t appears in no right-hand side and
% in no cost term. Pontryagin's costate ODE then gives dlambda_t/dL = -dH/dt = 0,
% so lambda_t (the costate of the time state) is a FIRST INTEGRAL: constant along
% the optimum. Its constant value is (minus) the time-domain Hamiltonian. Hence
% "the Hamiltonian is constant along the optimal trajectory" is EXACTLY
% "the recovered time-costate lambda_t is flat".
%
% THE DISCRETE FORM (sharpest, checked here). In the collocation NLP the row-7
% (time) defect couples only adjacent time nodes and depends on no state's t.
% KKT stationarity w.r.t. each interior time node t_k therefore forces the
% row-7 interval defect duals to be EQUAL across all intervals:
%     lamDef(7,k-1) = lamDef(7,k)   for every interior node.
% So lamDef(7,:) should be constant to solver (KKT) tolerance. A non-conserved
% costate (mass lambda_m, or lambda_P) is reported alongside as an internal
% control -- it should be visibly NON-constant, showing the test discriminates.
%
% INPUTS:
%   matPath  - path to a run_transfer/homotopy results .mat holding `res`
%              (uses res.fuel.lamDef [7xN], res.sigma [(N+1)x1], res.fuel.dL),
%              OR the res struct itself [char | struct]
%   plotPath - optional PNG path for the lambda_t-vs-revs diagnostic plot
%              (skipped if omitted/empty) [char]
%
% OUTPUTS:
%   chk - struct with fields:
%           .lamt_mean .lamt_std        mean/std of the row-7 interval duals
%           .lamt_cov                   coefficient of variation std/|mean|
%           .lamt_maxreldev             max_k |lamDef(7,k)-mean|/|mean|
%           .Hconst                     the constant time-Hamiltonian value = -mean(lambda_t)
%           .cov_mass .cov_P            control CoVs (non-conserved costates)
%           .discrim                    cov_mass / lamt_cov (>>1 => test discriminates)
%           .verdict                    'CONSTANT' | 'APPROX-CONSTANT' | 'NOT-CONSTANT'
%           .nNodes                     N+1
%
% REFERENCES:
%   [1] Pontryagin et al., Mathematical Theory of Optimal Processes (H is a
%       first integral for autonomous problems; constant, generally nonzero
%       for fixed final time).
%   [2] verify/mee_dual_to_costate.m (interval-dual -> nodal-costate recovery).
%   [3] earth_elliptic_to_geo/core/lt_mee_rhs.m (state ordering; t is X(7),
%       autonomous in t).
if ischar(matPath) || isstring(matPath)
    S = load(char(matPath));  res = S.res;
else
    res = matPath;
end
assert(isfield(res,'fuel') && isfield(res.fuel,'lamDef') && ~isempty(res.fuel.lamDef), ...
    'hamiltonian_const_check: res.fuel.lamDef is required (defect duals)');
LamDef = res.fuel.lamDef;                 % [7 x N] interval defect duals
sigma  = res.sigma(:);                     % [(N+1) x 1]
dL     = res.fuel.dL;
N      = size(LamDef, 2);

% --- row-7 (time-state) interval duals: should be constant across intervals ---
lamt   = LamDef(7, :);
mn     = mean(lamt);   sd = std(lamt);
cov    = sd / max(abs(mn), realmin);
maxrel = max(abs(lamt - mn)) / max(abs(mn), realmin);
Hconst = -mn;                              % time-Hamiltonian = -lambda_t (constant)

% --- controls: non-conserved costates (mass row 6, semilatus row 1) ---
lamm   = LamDef(6, :);   covM = std(lamm) / max(abs(mean(lamm)), realmin);
lamP   = LamDef(1, :);   covP = std(lamP) / max(abs(mean(lamP)), realmin);
discrim = covM / max(cov, realmin);

% --- verdict (KKT tol ~1e-6..1e-8; be generous, the control must dwarf it) ---
if     cov < 1e-4 && discrim > 50,  verdict = 'CONSTANT';
elseif cov < 1e-2 && discrim > 10,  verdict = 'APPROX-CONSTANT';
else,                               verdict = 'NOT-CONSTANT';
end

chk = struct('lamt_mean', mn, 'lamt_std', sd, 'lamt_cov', cov, ...
    'lamt_maxreldev', maxrel, 'Hconst', Hconst, 'cov_mass', covM, ...
    'cov_P', covP, 'discrim', discrim, 'verdict', verdict, 'nNodes', N+1);

fprintf(['[H-const] N+1=%d | lambda_t: mean=%.6g std=%.3g CoV=%.3e maxRelDev=%.3e\n' ...
         '          H(time)=-lambda_t=%.6g (const) | control CoV: mass=%.3e P=%.3e | ' ...
         'discriminate x%.0f | VERDICT: %s\n'], ...
        N+1, mn, sd, cov, maxrel, Hconst, covM, covP, discrim, verdict);

% --- optional diagnostic plot: lambda_t (flat) vs a non-conserved costate ---
if nargin >= 2 && ~isempty(plotPath)
    Lmid = pi + 0.5*(sigma(1:end-1)+sigma(2:end))*dL;    % interval-midpoint L
    revs = (Lmid - pi) / (2*pi);
    fig = figure('Visible','off','Color','w','Position',[100 100 1000 640]);
    tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; grid on;
    plot(revs, lamt, '-', 'Color',[0.10 0.45 0.85], 'LineWidth',1.2);
    yline(mn, '--', 'Color',[0.85 0.30 0.20], 'LineWidth',1.0);
    ylabel('\lambda_t  (time costate)');
    title(sprintf('Hamiltonian constancy:  \\lambda_t  flat  (CoV = %.2e)  \\Rightarrow  H = %.5g = const', cov, Hconst));
    % zoom y to show it really is flat at fractional scale
    yr = max(3*sd, abs(mn)*1e-6);
    ylim([mn-8*yr, mn+8*yr]);

    nexttile; hold on; grid on;
    plot(revs, lamm/max(abs(mean(lamm)),realmin), '-', 'Color',[0.20 0.60 0.25], 'LineWidth',1.0);
    xlabel('revolutions'); ylabel('\lambda_m / |mean|  (control)');
    title(sprintf('Internal control: mass costate \\lambda_m is NOT constant (CoV = %.2e, %.0f\\times larger)', covM, discrim));

    exportgraphics(fig, plotPath, 'Resolution', 150);
    close(fig);
    fprintf('          plot -> %s\n', plotPath);
end
end
