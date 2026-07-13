% DIRECT_MINTIME_ELFO  Direct (collocation) min-time GTO->ELFO by homotoping the
% target from the tulip. Uses the existing fmincon direct min-time NLP
% (attic/solve_tfmin_nlp: minimize tf, always-burn, trapezoid collocation).
%
% Why direct: the campaign's core lesson is that the DIRECT NLP converges where
% indirect shooting does not -- and retargeting to the ELFO is exactly where
% indirect (even multiple-shooting) min-time fought the shooting sensitivity.
% Direct collocation has all node states as variables (no STM products), so a
% target homotopy is robust; and min-time lets tf FLOAT, removing the fixed-tf
% "can't-reach-terminal" wall that stalled the direct ENERGY homotopy.
%
% Base vars: Nnodes [300], tgtName ['nearest'] ('apolune'|'nearest'|'both').

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','attic'));      % build_guess, solve_tfmin_nlp, nlp_constraints, lt_dynamics
resDir = fullfile(here,'results');

if evalin('base','exist(''Nnodes'',''var'')'), N = evalin('base','Nnodes'); else, N = 300; end
if evalin('base','exist(''tgtName'',''var'')'), tgtSel = evalin('base','tgtName'); else, tgtSel = 'nearest'; end

[rv0, rvf_tul, P] = mintime_params();  Tmax = P.Tmax25;  c = P.c;  mu = P.muStar;  lStar = P.lStar;

% --- tulip direct min-time (validate the solver + get the warm start) -------
fprintf('=== tulip DIRECT min-time (N=%d) ===\n', N);
[Z0, sigma] = build_guess('indirect', N, rv0, rvf_tul, Tmax, c, mu);
[Ztul, out] = solve_tfmin_nlp(Z0, sigma, rv0, rvf_tul, Tmax, c, mu);
fprintf('  tf=%.5f ND (%.2f d)  maxDefect=%.2e  exitflag=%d\n', ...
        out.tf, out.tf*P.tStar/86400, out.maxDefect, out.exitflag);
save(fullfile(resDir,'direct_mintime_tulip.mat'), 'Ztul','sigma','out','rv0','rvf_tul','P');
if out.maxDefect > 1e-6
    fprintf('  tulip direct min-time did NOT converge cleanly; stopping.\n');  return
end

% --- ELFO targets -----------------------------------------------------------
Mmass = 5.9736e24 + 7.35e22;  nu0 = pumpkyn.util.mean2True(0, 0.69);
oev = [12000, 0.69, 56.5*pi/180, 90*pi/180, 0, nu0];
x0  = pumpkyn.cr3bp.fromOrb(0, oev, Mmass, mu, P.tStar, lStar, 2);
tau = linspace(0,0.5,8000)';  [~, xs] = pumpkyn.cr3bp.prop(tau, x0, mu);
rM = [1-mu 0 0];  dM = sqrt(sum((xs(:,1:3)-rM).^2,2));
[~,iApo]=max(dM);  rvApo = xs(iApo,1:6);
d6 = sqrt(sum((xs(:,1:6)-rvf_tul(:).').^2,2)); [~,iNear]=min(d6); rvNear = xs(iNear,1:6);
allT = {'apolune', rvApo; 'nearest', rvNear};
if ~strcmp(tgtSel,'both'), allT = allT(strcmp(allT(:,1),tgtSel), :); end

% --- homotope the target: tulip -> ELFO (direct re-solves, tf floats) -------
for ti = 1:size(allT,1)
    nm = allT{ti,1};  rvf_elfo = allT{ti,2};
    fprintf('\n=== DIRECT min-time homotopy tulip -> ELFO %s (dMoon=%.0f km, ||drvf||=%.3f) ===\n', ...
            nm, norm(rvf_elfo(1:3)-rM)*lStar, norm(rvf_elfo(:)-rvf_tul(:)));
    Zc = Ztul;  s = 0;  step = 0.1;
    while s < 1 - 1e-9
        sTry = min(s+step, 1);
        rvf_s = (1-sTry)*rvf_tul(:).' + sTry*rvf_elfo(:).';
        [Z2, o2] = solve_tfmin_nlp(Zc, sigma, rv0, rvf_s, Tmax, c, mu);
        if o2.exitflag > 0 && o2.maxDefect < 1e-6
            Zc = Z2;  s = sTry;
            fprintf('  s=%.3f  tf=%.4f (%.2f d)  maxDefect=%.2e  exit=%d\n', ...
                    s, o2.tf, o2.tf*P.tStar/86400, o2.maxDefect, o2.exitflag);
            if step < 0.1, step = min(1.5*step, 0.1); end
        else
            step = step/2;
            fprintf('  s=%.3f STALL (defect=%.2e exit=%d) -> step=%.3f\n', sTry, o2.maxDefect, o2.exitflag, step);
            if step < 0.01, fprintf('  -> direct min-time homotopy stuck for %s\n', nm); break; end
        end
    end
    if s >= 1-1e-9
        save(fullfile(resDir,sprintf('direct_mintime_elfo_%s.mat',nm)), 'Zc','sigma','rv0','rvf_elfo','P');
        fprintf('  ELFO %s REACHED (direct min-time). saved.\n', nm);
    end
end
