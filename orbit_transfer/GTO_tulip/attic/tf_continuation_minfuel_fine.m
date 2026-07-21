function results = tf_continuation_minfuel_fine(N, factors, saveFile)
% TF_CONTINUATION_MINFUEL_FINE  Refined tf-homotopy of the min-fuel switch
% structure -- finer mesh, higher per-step iteration cap, and DENSE tf steps
% through the 1.10->1.20 breakdown where the coarse run (N=2500) lost tight
% convergence (6 clean switches -> 35 loose). Goal: grow the switch structure
% gradually enough to stay machine-tight into the many-switch regime.
%
% Identical method to TF_CONTINUATION_MINFUEL; only the resolution knobs
% change (N=6000 default, MaxIterations 4000/step, dense factors). Writes to
% a SEPARATE results file so it can run alongside the coarse pass.
%
% INPUTS/OUTPUTS: see TF_CONTINUATION_MINFUEL.
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.  [2] Bertrand & Epenoy, OCAM 23(4), 2002.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(N), N = 6000; end
if nargin < 2 || isempty(factors)
    % dense through the 1.09-1.19 breakdown, coarser outside
    factors = [1.00 1.04 1.07 1.09 1.10 1.115 1.13 1.145 1.16 1.175 1.19 1.21 1.25 1.30];
end
if nargin < 3 || isempty(saveFile)
    saveFile = fullfile(here, 'tf_continuation_fine_results.mat');
end

run(fullfile(here, 'setup_paths.m'));
addpath(fullfile(here, '..', 'indirect', 'lowThrust_GTO_tulip'));

muStar = 0.012150585609624; lStar = 389703.264829278; tStar = 382981.289129055;
m0kg = 15; g0 = 9.80665*tStar^2/(1000*lStar);
Tmax = (0.025/m0kg)*tStar^2/(lStar*1000); c = (2100/tStar)*g0;

muEarth = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2; ecc = (35786 - 350)/(2*sma);
[r0d, v0d] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0d, v0d], muStar, tStar, lStar, 1);
[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, rvTgt]   = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, muStar);
[~, idxF]    = max(rvTgt(:,5));
rvf = rvTgt(idxF, :);

tfMin = 6.2906939607;
zMinTime = [190.4760481; -79.7060409; -0.4298691037; 0.3011592775; ...
            0.5866700046; -0.007117348902; 4.329378839];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[tauMT, ~] = ode113(@lt_pmp_eom, [0 tfMin], [rv0(:); 1; zMinTime], opts, Tmax, c, muStar);
tauFine = linspace(0, tfMin, 30000).';
[tauFine, yFine] = ode113(@lt_pmp_eom, tauFine, [rv0(:); 1; zMinTime], opts, Tmax, c, muStar);
[tauFine, keep] = unique(tauFine, 'stable'); yFine = yFine(keep, :);

sigAd = unique(tauMT)/tfMin;
sigma = interp1(linspace(0,1,numel(sigAd)).', sigAd, linspace(0,1,N+1).');
sigma(1) = 0; sigma(end) = 1;

tMesh0 = sigma*tfMin;
Xg   = interp1(tauFine, yFine(:,1:7), tMesh0, 'pchip').';
lamV = interp1(tauFine, yFine(:,11:13), tMesh0, 'pchip').';
alph = -lamV./sqrt(sum(lamV.^2, 1));
Xg(1:6,1) = rv0(:); Xg(7,1) = 1; Xg(1:6,end) = rvf(:);
s0 = ones(1, N+1);
Z  = [Xg(:); reshape([alph.*s0; s0], [], 1)];

results = struct('factor',{},'tf',{},'switches',{},'mProp_kg',{}, ...
                 'maxDefect',{},'flag',{},'burnFrac',{},'converged',{});
Zgood = Z;
for f = factors
    tf = f*tfMin;
    ok = false;
    try
        [Znew, nlp] = solveMinfuelStepFine(Zgood, sigma, tf, rv0, 1, rvf, Tmax, c, muStar);
        s  = nlp.U(4,:);
        sw = sum(abs(diff(s > 0.5)));
        conv = nlp.maxDefect < 5e-3;
        results(end+1) = struct('factor',f,'tf',tf,'switches',sw, ...
            'mProp_kg',m0kg*(1-nlp.mf),'maxDefect',nlp.maxDefect,'flag',nlp.exitflag, ...
            'burnFrac',mean(s>0.5),'converged',conv); %#ok<AGROW>
        fprintf('f=%.3f tf=%.3f ND (%.2f d): switches=%3d  prop=%.4f kg  defect=%.2g  flag=%d  %s\n', ...
                f, tf, tf*tStar/86400, sw, m0kg*(1-nlp.mf), nlp.maxDefect, nlp.exitflag, ...
                ternary(conv,'OK','(loose)'));
        if conv, Zgood = Znew; ok = true; end
    catch meErr
        fprintf('f=%.3f FAILED: %s\n', f, meErr.message);
    end
    save(saveFile, 'results', 'sigma', 'Zgood', 'factors', 'N');
    if ~ok
        fprintf('  (kept previous warm start; step %.3f did not converge cleanly)\n', f);
    end
end

fprintf('\n=== fine tf-continuation summary (switches vs tf) ===\n');
for k = 1:numel(results)
    fprintf('  f=%.3f  switches=%3d  prop=%.4f kg  defect=%.2g  %s\n', results(k).factor, ...
            results(k).switches, results(k).mProp_kg, results(k).maxDefect, ...
            ternary(results(k).converged,'','LOOSE'));
end
end

% ---------------------------------------------------------------------------
function [Z, out] = solveMinfuelStepFine(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar)
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
    'MaxIterations',4000,'MaxFunctionEvaluations',3e5, ...
    'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-7, ...
    'StepTolerance',1e-12,'Display','off');
[Z,~,ef,op] = fmincon(objFun, Z0, [],[],[],[], lb, ub, conFun, o); %#ok<ASGLU>
X = reshape(Z(1:7*nN),7,nN); U = reshape(Z(7*nN+(1:4*nN)),4,nN);
[~,ceq] = nlp_constraints_minfuel(Z,sigma,tf,Tmax,c,muStar);
% maxDefect now spans ALL constraints (dynamics defects AND the throttle
% cone). The old ceq(1:7*N)-only version let cone-loose solutions (cone off
% by ~2e-2) pass the guard -- fixed.
out = struct('X',X,'U',U,'mf',X(7,end),'exitflag',ef, ...
             'maxDefect',max(abs(ceq)), ...
             'maxDynDefect',max(abs(ceq(1:7*N))), 'maxCone',max(abs(ceq(7*N+1:end))));
end

function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end
