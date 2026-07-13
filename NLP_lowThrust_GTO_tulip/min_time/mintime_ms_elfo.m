% MINTIME_MS_ELFO  Homotope the min-time MS target from the tulip to the ELFO.
%
% Starts from the converged tulip min-time MS solution (mintime_tulip_ms.mat,
% ||R||~4e-9) and slides rvf tulip->ELFO in adaptive steps, re-solving the MS
% each step (warm-started). Min-time's tf FLOATS, so every intermediate target
% is reachable -- this is the robustness the fixed-tf energy homotopy lacked.
% Reaches both apolune (south-pole dwell point) and nearest-insertion; saves
% each converged min-time ELFO solution to seed a min-energy solve (-> PSR).
%
% Base vars: tgtName ['apolune'] ('apolune'|'nearest'|'both'), step0 [0.1].

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','ztl'));
warning('off','MATLAB:ode45:IntegrationTolNotMet');
resDir = fullfile(here,'results');

if evalin('base','exist(''tgtName'',''var'')'), tgtSel = evalin('base','tgtName'); else, tgtSel = 'apolune'; end
if evalin('base','exist(''step0'',''var'')'),   step0  = evalin('base','step0');   else, step0  = 0.1;       end

S = load(fullfile(resDir,'mintime_tulip_ms.mat'));
z_tul = S.z;  prob0 = S.prob;  prob0.resFun = @mintime_ms_residual;
rvf_tul = prob0.rvf(:).';
[~, ~, P] = mintime_params();  mu = prob0.muStar;  lStar = P.lStar;

% --- ELFO targets (apolune + nearest-to-tulip) ------------------------------
Mmass = 5.9736e24 + 7.35e22;  nu0 = pumpkyn.util.mean2True(0, 0.69);
oev = [12000, 0.69, 56.5*pi/180, 90*pi/180, 0, nu0];
x0  = pumpkyn.cr3bp.fromOrb(0, oev, Mmass, mu, P.tStar, lStar, 2);
tau = linspace(0,0.5,8000)';  [~, xs] = pumpkyn.cr3bp.prop(tau, x0, mu);
rM = [1-mu 0 0];  dM = sqrt(sum((xs(:,1:3)-rM).^2,2));
[~,iApo]=max(dM);  rvApo = xs(iApo,1:6);
d6 = sqrt(sum((xs(:,1:6)-rvf_tul).^2,2)); [~,iNear]=min(d6); rvNear = xs(iNear,1:6);

allT = {'apolune', rvApo; 'nearest', rvNear};
if ~strcmp(tgtSel,'both'), allT = allT(strcmp(allT(:,1),tgtSel), :); end

for ti = 1:size(allT,1)
    nm = allT{ti,1};  rvf_elfo = allT{ti,2};
    fprintf('\n=== MIN-TIME MS homotopy tulip -> ELFO %s (dMoon=%.0f km, ||drvf||=%.3f) ===\n', ...
            nm, norm(rvf_elfo(1:3)-rM)*lStar, norm(rvf_elfo(:)-rvf_tul(:)));
    % PREDICTOR-CORRECTOR (tangent) continuation. Min-time is sensitive: a small
    % target move shifts the solution far in z-space, so a plain warm start lands
    % in the slow globalization phase. Instead, use the tangent dz/ds (from the
    % Jacobian at the current solution) to PREDICT the next iterate near the new
    % branch, then correct -- a few iterations each.
    drvf = rvf_elfo(:).' - rvf_tul;            % d(rvf_s)/ds (constant)
    rendRows = 14*(prob0.M-1) + (1:6);
    zc = z_tul;  prob = prob0;  s = 0;  step = step0;
    while s < 1 - 1e-9
        sTry = min(s+step, 1);  ds = sTry - s;
        % tangent predictor at the current converged (zc, rvf_s)
        prob.rvf = (1-s)*rvf_tul + s*rvf_elfo(:).';
        [~, J, ~] = mintime_ms_residual(zc, prob, true);
        rhs = zeros(size(J,1),1);  rhs(rendRows) = drvf(:);   % J dz/ds = drvf (rend rows)
        dzds = J \ rhs;
        z_pred = zc + ds*dzds;
        % corrector at the new target
        prob.rvf = (1-sTry)*rvf_tul + sTry*rvf_elfo(:).';
        [z2, out] = ztl_ms_solve_tr(z_pred, prob, struct('tolR',1e-8,'maxIter',80,'verbose',false));
        if out.resNorm < 1e-6
            zc = z2;  s = sTry;
            fprintf('  s=%.3f  ||R||=%.3e  tf=%.4f (%.2f d)  iters=%d\n', ...
                    s, out.resNorm, zc(end), zc(end)*P.tStar/86400, out.iters);
            if step < step0, step = min(1.5*step, step0); end
        else
            step = step/2;
            fprintf('  s=%.3f STALL (||R||=%.2e) -> step=%.3f\n', sTry, out.resNorm, step);
            if step < 0.005, fprintf('  -> min-time homotopy stuck for %s\n', nm); break; end
        end
    end
    if s >= 1-1e-9
        probF = prob;  zF = zc;
        save(fullfile(resDir,sprintf('mintime_elfo_%s_ms.mat',nm)), 'zF','probF','P');
        fprintf('  ELFO %s REACHED (min-time MS): tf=%.5f ND (%.2f d). saved.\n', ...
                nm, zc(end), zc(end)*P.tStar/86400);
    end
end
