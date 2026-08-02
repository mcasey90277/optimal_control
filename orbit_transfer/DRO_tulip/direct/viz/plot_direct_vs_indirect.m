function F = plot_direct_vs_indirect(o, tauInd, zInd, p, Tmax, c, outStem)
% PLOT_DIRECT_VS_INDIRECT  How closely the direct solve matches the indirect one.
%
% TWO FIGURES, because there are two separate claims and conflating them is how
% this campaign previously overstated its results.
%
%   FIGURE 1, STATE. Do the two methods fly the same trajectory? Note that
%   agreement at the endpoints is guaranteed by the boundary conditions and is
%   therefore worth nothing as evidence; what matters is the INTERIOR, where
%   nothing forces the two to agree.
%
%   FIGURE 2, COSTATE. Do they find the same extremal? This is the stronger
%   claim: the costates are never constrained to match, never compared during
%   either solve, and the direct method does not even represent them -- they are
%   recovered from NLP multipliers after the fact.
%
% A NOTE ON WHAT THE TIME AXIS HIDES. On this transfer the closest lunar
% approach sits at the very END of the trajectory, so any quantity that degrades
% near periselene also appears to degrade "at the boundary". Every panel that
% shows a rise at late time therefore also carries the lunar altitude, so the
% two can be told apart by eye. Failing to separate them is exactly the error
% that produced a wrong diagnosis here once already.
%
% INPUTS:
%   o       - direct solution from casadi_mintime_dro (returnModel = true for
%             the costate figure)
%   tauInd  - indirect time history [Mx1]
%   zInd    - indirect state+costate history [Mx14] from tfMinProp:
%             columns 1:7 = [r;v;m], 8:14 = [lambda_r; lambda_v; lambda_m]
%   p       - params (.muStar .lStar .tStar)
%   Tmax, c - ND thrust acceleration and exhaust velocity [scalars]
%   outStem - path stem; '_state.png' and '_costate.png' are appended [char]
%
% OUTPUTS:
%   F - struct: .posErrKm .velErrMs [1x(N+1)] interior state differences,
%       .posErrMaxKm .velErrMaxMs, .files {2x1}
%
% REFERENCES:
%   [1] certify/costate_compare.m -- the numbers these figures visualize.

mu = p.muStar;  lStar = p.lStar;  tStar = p.tStar;
rMoonKm = 1737.4;
t   = o.s(:).' * o.tf;
nN  = numel(t);

% resample the indirect solution onto the direct time grid
zI  = interp1(tauInd(:), zInd, t(:), 'spline').';        % 14 x nN
altKm = (vecnorm(o.X(1:3,:) - [1-mu;0;0], 2, 1))*lStar - rMoonKm;

dPos = vecnorm(o.X(1:3,:) - zI(1:3,:), 2, 1) * lStar;              % km

% THE SAME COMPARISON WITH PHASING REMOVED. The two solutions have slightly
% different final times (here 7.6e-06 ND = 2.9 s), and comparing them at matched
% ABSOLUTE time therefore charges that offset as a position error: at ~1.1 km/s
% it is worth up to 3.3 km of pure along-track displacement, which is most of
% what the raw difference shows. Re-sampling the indirect solution at matched
% FRACTIONAL time s = t/t_f removes it and leaves the difference in trajectory
% SHAPE, which is the quantity of interest.
zF   = interp1(tauInd(:)/tauInd(end), zInd, o.s(:), 'spline').';
dPosF= vecnorm(o.X(1:3,:) - zF(1:3,:), 2, 1) * lStar;
F_dtfSec = (o.tf - tauInd(end)) * tStar;
dVel = vecnorm(o.X(4:6,:) - zI(4:6,:), 2, 1) * lStar/tStar * 1000; % m/s
dMas = abs(o.X(7,:) - zI(7,:));

F = struct('posErrKm',dPos,'velErrMs',dVel,'massErr',dMas, ...
           'posErrMaxKm',max(dPos),'velErrMaxMs',max(dVel), ...
           'posErrPhasedKm',dPosF,'posErrPhasedMaxKm',max(dPosF), ...
           'dtfSec',F_dtfSec,'files',{{}});

% =========================== FIGURE 1: STATE ==============================
f1 = figure('Color','w','Position',[40 40 1500 900]);

subplot(2,3,1);
plot((o.X(1,:)-(1-mu))*lStar, o.X(2,:)*lStar, '-', 'LineWidth',1.6, ...
     'Color',[0.15 0.35 0.75]); hold on; grid on; axis equal;
plot((zI(1,:)-(1-mu))*lStar, zI(2,:)*lStar, '--', 'LineWidth',1.0, ...
     'Color',[0.90 0.35 0.15]);
th = linspace(0,2*pi,200);
fill(rMoonKm*cos(th), rMoonKm*sin(th), [0.6 0.6 0.6], 'EdgeColor','none');
xlabel('x - x_{Moon} [km]'); ylabel('y [km]');
legend('direct','indirect','Moon','Location','best');
title('(a) trajectory overlay (they superpose)');

subplot(2,3,2);
semilogy(t*tStar/86400, max(dPos,1e-12), '-', 'LineWidth',1.4); hold on;
semilogy(t*tStar/86400, max(dPosF,1e-12), '-', 'LineWidth',1.4);
grid on; xlabel('time [days]'); ylabel('position difference [km]');
legend('matched absolute time','matched fractional time','Location','southeast');
title(sprintf(['(b) POSITION difference: %.3f km raw, %.3f km once the\n' ...
    '%.1f s difference in t_f is removed -- it is PHASING, not shape'], ...
    max(dPos), max(dPosF), F_dtfSec));

subplot(2,3,3);
yyaxis left;  semilogy(t*tStar/86400, max(dVel,1e-12), '-', 'LineWidth',1.3);
ylabel('|v_{direct} - v_{indirect}| [m/s]');
yyaxis right; semilogy(t*tStar/86400, altKm, '-', 'LineWidth',1.0);
ylabel('lunar altitude [km]');
grid on; xlabel('time [days]');
title(sprintf('(c) VELOCITY difference, max %.4f m/s', max(dVel)));

subplot(2,3,4);
plot(t*tStar/86400, o.X(7,:), '-', 'LineWidth',1.6); hold on; grid on;
plot(t*tStar/86400, zI(7,:), '--', 'LineWidth',1.0);
xlabel('time [days]'); ylabel('mass fraction');
legend('direct','indirect','Location','best');
title(sprintf('(d) mass, max difference %.2e', max(dMas)));

subplot(2,3,5);
aD = o.U(1:3,:) ./ max(vecnorm(o.U(1:3,:),2,1), eps);
lv = zI(11:13,:);  aI = -lv ./ max(vecnorm(lv,2,1), eps);
angT = real(acosd(max(-1,min(1, sum(aD.*aI,1)))));
semilogy(t*tStar/86400, max(angT,1e-10), '-', 'LineWidth',1.3); grid on;
xlabel('time [days]'); ylabel('angle [deg]');
title(sprintf('(e) THRUST DIRECTION difference\nmedian %.2e deg, max %.2e deg', ...
    median(angT), max(angT)));

subplot(2,3,6);
loglog(altKm, max(dPos,1e-12), '.', 'MarkerSize',5); grid on;
xlabel('lunar altitude [km]'); ylabel('position difference [km]');
title('(f) difference vs altitude, not vs time');

sgtitle(sprintf(['DRO \\rightarrow tulip: DIRECT vs INDIRECT, STATE  |  N = %d, %s  |  ' ...
    't_f  %.7f vs %.7f'], nN-1, o.scheme, o.tf, tauInd(end)), 'FontWeight','bold');
file1 = [outStem '_state.png'];
exportgraphics(f1, file1, 'Resolution', 140);  close(f1);
F.files{end+1} = file1;

% ========================== FIGURE 2: COSTATE =============================
if isfield(o,'lamDef') && ~isempty(o.lamDef)
    N    = size(o.lamDef,2);
    tMid = 0.5*(t(1:end-1) + t(2:end));
    lamI = interp1(tauInd(:), zInd(:,8:14), tMid(:), 'spline').';
    % sign from the primal control, exactly as costate_compare does
    if isfield(o,'Um') && ~isempty(o.Um), aS = o.Um(1:3,:);
    else, aS = 0.5*(o.U(1:3,1:end-1)+o.U(1:3,2:end)); end
    aS = aS ./ max(vecnorm(aS,2,1),eps);
    fit = zeros(1,2);  sgnTry = [1 -1];
    for kk = 1:2
        a = -sgnTry(kk)*o.lamDef(4:6,:) ./ max(vecnorm(o.lamDef(4:6,:),2,1),eps);
        fit(kk) = mean(sum(a.*aS,1));
    end
    [~,ib] = max(fit);  lamD = sgnTry(ib)*o.lamDef;
    altM = interp1(t, altKm, tMid);

    f2 = figure('Color','w','Position',[40 40 1500 900]);
    lbl = {'\lambda_{r}','\lambda_{v}'};
    for blk = 1:2
        subplot(2,3,blk);
        rows = (1:3) + 3*(blk-1);
        co = lines(3);
        for j = 1:3
            plot(tMid*tStar/86400, lamD(rows(j),:), '-', 'LineWidth',1.6, 'Color',co(j,:)); hold on;
            plot(tMid*tStar/86400, lamI(rows(j),:), '--', 'LineWidth',1.0, 'Color',co(j,:)*0.6);
        end
        grid on; xlabel('time [days]'); ylabel(lbl{blk});
        title(sprintf('(%c) %s components: solid = direct duals,\ndashed = indirect costates', ...
            'a'+blk-1, lbl{blk}));
    end

    subplot(2,3,3);
    plot(tMid*tStar/86400, lamD(7,:), '-', 'LineWidth',1.6); hold on; grid on;
    plot(tMid*tStar/86400, lamI(7,:), '--', 'LineWidth',1.0);
    yline(0,'k:');
    xlabel('time [days]'); ylabel('\lambda_m');
    legend('direct','indirect','Location','best');
    title(sprintf('(c) \\lambda_m: transversality wants 0 at t_f\ndirect end %.2e', lamD(7,end)));

    subplot(2,3,4);
    uD = lamD(4:6,:) ./ max(vecnorm(lamD(4:6,:),2,1),eps);
    uI = lamI(4:6,:) ./ max(vecnorm(lamI(4:6,:),2,1),eps);
    angv = real(acosd(max(-1,min(1,sum(uD.*uI,1)))));
    uDr = lamD(1:3,:) ./ max(vecnorm(lamD(1:3,:),2,1),eps);
    uIr = lamI(1:3,:) ./ max(vecnorm(lamI(1:3,:),2,1),eps);
    angr = real(acosd(max(-1,min(1,sum(uDr.*uIr,1)))));
    semilogy(tMid*tStar/86400, max(angv,1e-10), '-', 'LineWidth',1.3); hold on;
    semilogy(tMid*tStar/86400, max(angr,1e-10), '-', 'LineWidth',1.3);
    grid on; legend('\lambda_v','\lambda_r','Location','best');
    xlabel('time [days]'); ylabel('angle to indirect [deg]');
    title(sprintf('(d) costate DIRECTION agreement\nmedian %.1e / %.1e deg', ...
        median(angv), median(angr)));

    subplot(2,3,5);
    ratio = vecnorm(lamI(4:6,:),2,1) ./ max(vecnorm(lamD(4:6,:),2,1),eps);
    plot(tMid*tStar/86400, ratio, '-', 'LineWidth',1.3); grid on;
    yline(1,'k--');
    xlabel('time [days]'); ylabel('||\lambda_v^{ind}|| / ||\lambda_v^{dir}||');
    title(sprintf('(e) SCALE factor: median %.6f\nCoV %.2e -- 1.0 means same normalization', ...
        median(ratio), std(ratio)/abs(mean(ratio))));

    subplot(2,3,6);
    if isfield(o,'Xm') && ~isempty(o.Xm), Xmid = o.Xm;
    else, Xmid = 0.5*(o.X(:,1:end-1)+o.X(:,2:end)); end
    if isfield(o,'Um') && ~isempty(o.Um), thr = o.Um(4,:); else, thr = ones(1,N); end
    H = zeros(1,N);
    for k = 1:N
        H(k) = lamD(:,k).' * local_f(Xmid(:,k), [aS(:,k); thr(k)], mu, Tmax, c);
    end
    yyaxis left;  semilogy(tMid*tStar/86400, max(abs(H+1),1e-14), '-', 'LineWidth',1.3);
    ylabel('|\lambda^Tf + 1|   (PMP wants 0)');
    yyaxis right; semilogy(tMid*tStar/86400, altM, '-', 'LineWidth',1.0);
    ylabel('lunar altitude [km]');
    grid on; xlabel('time [days]');
    title(sprintf('(f) HAMILTONIAN condition, median %.1e\nthe late rise tracks ALTITUDE, not the boundary', ...
        median(abs(H+1))));

    sgtitle('DRO \rightarrow tulip: DIRECT duals vs INDIRECT costates', 'FontWeight','bold');
    file2 = [outStem '_costate.png'];
    exportgraphics(f2, file2, 'Resolution', 140);  close(f2);
    F.files{end+1} = file2;
end

for k = 1:numel(F.files), fprintf('  figure -> %s\n', F.files{k}); end
end

% ---------------------------------------------------------------------------
function dz = local_f(z, u, mu, Tmax, c)
% LOCAL_F  CR3BP with thrust; mirrors certify/dro_residual.
% INPUTS: z [7x1]; u [4x1]; mu; Tmax; c   OUTPUTS: dz [7x1]
r = z(1:3);  v = z(4:6);  m = z(7);
al = u(1:3);  al = al/max(norm(al),eps);  th = min(max(u(4),0),1);
dd = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
gr = [r(1) - (1-mu)*(r(1)+mu)/dd^3 - mu*(r(1)-1+mu)/rr^3;
      r(2) - (1-mu)*r(2)/dd^3      - mu*r(2)/rr^3;
           - (1-mu)*r(3)/dd^3      - mu*r(3)/rr^3];
dz = [v; gr + [2*v(2); -2*v(1); 0] + (th*Tmax/m)*al; -(Tmax/c)*th];
end
