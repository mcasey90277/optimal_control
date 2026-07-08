function results = tf_continuation_minfuel(N, factors, saveFile)
% TF_CONTINUATION_MINFUEL  Grow the min-fuel switch structure by tf-homotopy.
%
% Energy->fuel continuation, extended in TIME. Starting from the min-time
% solution (always-burn, tf = tfMin) with the FIXED tulip target, step the
% transfer time up. With a fixed endpoint and growing time budget the
% min-fuel optimum spends the slack as per-revolution coasts, so the number
% of bang-bang switches grows from 0 (min-time) toward the research-grade
% ~80-switch regime. Each step is warm-started from the previous converged
% solution (consecutive tf are close, so the switch structure grows one
% relaxation at a time -- the only way a single-shot solve can't reach it).
%
% Robust for long unattended runs: per-step fmincon iteration cap, a defect
% guard (a step that fails to converge does NOT poison the next warm start),
% and an incremental save after every step.
%
% INPUTS:
%   N        - trapezoidal segments [default 2500]
%   factors  - increasing tf multipliers of tfMin [default see below]
%   saveFile - .mat path for incremental results [default alongside this file]
%
% OUTPUTS:
%   results - struct array, one per factor: .factor .tf .switches .mProp_kg
%             .maxDefect .flag .burnFrac .converged
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.  [2] Bertrand & Epenoy, OCAM 23(4), 2002.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(N), N = 2500; end
if nargin < 2 || isempty(factors)
    factors = [1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.50 1.60 1.75 1.90 2.00];
end
if nargin < 3 || isempty(saveFile)
    saveFile = fullfile(here, 'tf_continuation_results.mat');
end

run(fullfile(here, 'setup_paths.m'));
addpath(fullfile(here, '..', 'lowThrust_GTO_tulip'));

muStar = 0.012150585609624; lStar = 389703.264829278; tStar = 382981.289129055;
m0kg = 15; g0 = 9.80665*tStar^2/(1000*lStar);
Tmax = (0.025/m0kg)*tStar^2/(lStar*1000); c = (2100/tStar)*g0;

% --- endpoints: GTO start, FIXED tulip target (the true rendezvous point) ---
muEarth = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2; ecc = (35786 - 350)/(2*sma);
[r0d, v0d] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0d, v0d], muStar, tStar, lStar, 1);
[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, rvTgt]   = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, muStar);
[~, idxF]    = max(rvTgt(:,5));
rvf = rvTgt(idxF, :);

% --- min-time arc: density-matched mesh + fine grid for the warm start ------
tfMin = 6.2906939607;
zMinTime = [190.4760481; -79.7060409; -0.4298691037; 0.3011592775; ...
            0.5866700046; -0.007117348902; 4.329378839];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[tauMT, ~] = ode113(@lt_pmp_eom, [0 tfMin], [rv0(:); 1; zMinTime], opts, Tmax, c, muStar);
tauFine = linspace(0, tfMin, 20000).';
[tauFine, yFine] = ode113(@lt_pmp_eom, tauFine, [rv0(:); 1; zMinTime], opts, Tmax, c, muStar);
[tauFine, keep] = unique(tauFine, 'stable'); yFine = yFine(keep, :);

sigAd = unique(tauMT)/tfMin;
sigma = interp1(linspace(0,1,numel(sigAd)).', sigAd, linspace(0,1,N+1).');
sigma(1) = 0; sigma(end) = 1;

% warm start at tf = tfMin: min-time states, throttle 1 (always burn), primer dir
tMesh0 = sigma*tfMin;
Xg   = interp1(tauFine, yFine(:,1:7), tMesh0, 'pchip').';
lamV = interp1(tauFine, yFine(:,11:13), tMesh0, 'pchip').';
alph = -lamV./sqrt(sum(lamV.^2, 1));
Xg(1:6,1) = rv0(:); Xg(7,1) = 1; Xg(1:6,end) = rvf(:);
s0 = ones(1, N+1);
Z  = [Xg(:); reshape([alph.*s0; s0], [], 1)];

% --- continuation loop -----------------------------------------------------
results = struct('factor',{},'tf',{},'switches',{},'mProp_kg',{}, ...
                 'maxDefect',{},'flag',{},'burnFrac',{},'converged',{});
Zgood = Z;                         % last warm start that produced a good solve
for f = factors
    tf = f*tfMin;
    ok = false;
    try
        [Znew, nlp] = solveMinfuelStep(Zgood, sigma, tf, rv0, 1, rvf, Tmax, c, muStar);
        s  = nlp.U(4,:);
        sw = sum(abs(diff(s > 0.5)));
        conv = nlp.maxDefect < 5e-3;          % coherent solution at the mesh floor
        results(end+1) = struct('factor',f,'tf',tf,'switches',sw, ...
            'mProp_kg',m0kg*(1-nlp.mf),'maxDefect',nlp.maxDefect,'flag',nlp.exitflag, ...
            'burnFrac',mean(s>0.5),'converged',conv); %#ok<AGROW>
        fprintf('f=%.2f tf=%.3f ND (%.2f d): switches=%3d  prop=%.4f kg  defect=%.2g  flag=%d  %s\n', ...
                f, tf, tf*tStar/86400, sw, m0kg*(1-nlp.mf), nlp.maxDefect, nlp.exitflag, ...
                ternary(conv,'OK','(loose)'));
        if conv, Zgood = Znew; ok = true; end   % only advance the warm start on a good solve
    catch meErr
        fprintf('f=%.2f FAILED: %s\n', f, meErr.message);
    end
    save(saveFile, 'results', 'sigma', 'Zgood', 'factors', 'N');
    if ~ok
        fprintf('  (kept previous warm start; step %.2f did not converge cleanly)\n', f);
    end
end

fprintf('\n=== tf-continuation summary (switches vs tf) ===\n');
for k = 1:numel(results)
    fprintf('  f=%.2f  switches=%3d  prop=%.4f kg  %s\n', results(k).factor, ...
            results(k).switches, results(k).mProp_kg, ternary(results(k).converged,'','LOOSE'));
end
end

% ---------------------------------------------------------------------------
function [Z, out] = solveMinfuelStep(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar)
% One min-fuel NLP solve with a per-step iteration cap (good warm start ->
% few iterations; caps a stalling step so the continuation keeps moving).
sigma = sigma(:); N = numel(sigma)-1; nN = N+1; nZ = 11*nN;
lb = -inf(nZ,1); ub = inf(nZ,1);
for k = 1:nN
    xi=(k-1)*7+(1:7); lb(xi)=[-3;-3;-3;-12;-12;-12;0.3]; ub(xi)=[3;3;3;12;12;12;1.0];
    ui=7*nN+(k-1)*4+(1:4); lb(ui)=[-1.1;-1.1;-1.1;0]; ub(ui)=[1.1;1.1;1.1;1];
end
lb(1:7)=[rv0(:);m0]; ub(1:7)=[rv0(:);m0];
xf=(nN-1)*7+(1:6); lb(xf)=rvf(:); ub(xf)=rvf(:);
idxMf = 7*nN; gJ = sparse(nZ,1); gJ(idxMf) = -1;
objFun = @(Z) deal(-Z(idxMf), gJ);
conFun = @(Z) nlp_constraints_minfuel(Z, sigma, tf, Tmax, c, muStar);
o = optimoptions('fmincon','Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
    'InitBarrierParam',1e-4,'HessianApproximation','lbfgs', ...
    'MaxIterations',1500,'MaxFunctionEvaluations',1e5, ...
    'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-7, ...
    'StepTolerance',1e-12,'Display','off');
[Z,~,ef,op] = fmincon(objFun, Z0, [],[],[],[], lb, ub, conFun, o); %#ok<ASGLU>
X = reshape(Z(1:7*nN),7,nN); U = reshape(Z(7*nN+(1:4*nN)),4,nN);
[~,ceq] = nlp_constraints_minfuel(Z,sigma,tf,Tmax,c,muStar);
% maxDefect = worst violation over ALL constraints (dynamics defects AND the
% throttle cone). The earlier version used only ceq(1:7*N) (defects), which
% let the guard accept cone-loose "solutions" (cone off by ~2e-2) -- fixed.
out = struct('X',X,'U',U,'mf',X(7,end),'exitflag',ef, ...
             'maxDefect',max(abs(ceq)), ...
             'maxDynDefect',max(abs(ceq(1:7*N))), 'maxCone',max(abs(ceq(7*N+1:end))));
end

function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end
