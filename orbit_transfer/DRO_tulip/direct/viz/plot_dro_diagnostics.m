function D = plot_dro_diagnostics(o, p, Tmax, c, outPng)
% PLOT_DRO_DIAGNOSTICS  Six-panel diagnostic for a direct DRO->tulip solution.
%
% The panels answer distinct questions, and the ordering is by what a reader
% should check first:
%
%   (a) TRAJECTORY, Moon-centred, with the lunar disc and the altitude floor
%       drawn to scale. Whether the constraint bites is visible, not inferred.
%   (b) LUNAR ALTITUDE vs time, against the floor. The close approach is the
%       whole story of this problem, so it gets its own panel.
%   (c) TIME-GRID SPACING. The mesh is uniform in NORMALIZED time, so physical
%       spacing is uniform too -- which is exactly the weakness at periselene,
%       where the dynamics are fastest and the nodes are no denser.
%   (d) CONTINUOUS-TIME RESIDUAL. Each interval is re-integrated at high
%       accuracy from its own left node and compared at its right node. This is
%       the panel that says whether the collocation is honest between nodes:
%       the NLP drives the DEFECT to ~1e-14, but that is a statement about the
%       trapezoid rule, not about the trajectory.
%   (e) THROTTLE. Min-time should saturate at 1 throughout; anything else is a
%       finding.
%   (f) RESIDUAL vs LUNAR ALTITUDE. If corner-cutting is happening it shows up
%       here as the residual rising where the altitude falls.
%
% INPUTS:
%   o      - solution struct from casadi_mintime_dro
%   p      - params struct from dro_tulip_endpoints (.muStar .lStar .tStar)
%   Tmax,c - ND thrust acceleration and exhaust velocity [scalars]
%   outPng - (optional) write the figure here [char]
%
% OUTPUTS:
%   D - struct: .tPhys .altKm .hPhys .Rx (per-interval continuous residual)
%       .RxMax .RxMed .altAtWorst
%
% REFERENCES:
%   [1] verify_common/pmp/pmp_residual.m -- the same continuous-residual idea,
%       applied here inline because this problem has no switches to split on.

if nargin < 5, outPng = ''; end
mu = p.muStar;  lStar = p.lStar;  tStar = p.tStar;
rMoonKm = 1737.4;
N1 = size(o.X,2);  N = N1-1;
tPhys = o.s(:).' * o.tf;                       % ND time at the nodes
r2 = vecnorm(o.X(1:3,:) - [1-mu;0;0], 2, 1);
altKm = r2*lStar - rMoonKm;

% --- continuous residual: re-integrate each interval at high accuracy -------
odeo = odeset('RelTol',1e-11,'AbsTol',1e-13);
Rx = nan(1,N);
for k = 1:N
    a = tPhys(k);  b = tPhys(k+1);
    uOf = @(t) o.U(:,k) + (t-a)/max(b-a,eps)*(o.U(:,k+1)-o.U(:,k));
    [~,Z] = ode113(@(t,z) local_rhs(z, uOf(t), mu, Tmax, c), [a b], o.X(:,k), odeo);
    Rx(k) = norm(Z(end,:).' - o.X(:,k+1));
end
hPhys = diff(tPhys);
[RxMax, kw] = max(Rx);
D = struct('tPhys',tPhys,'altKm',altKm,'hPhys',hPhys,'Rx',Rx, ...
           'RxMax',RxMax,'RxMed',median(Rx),'altAtWorst',altKm(kw));

fig = figure('Color','w','Position',[60 60 1500 900]);

% (a) trajectory, Moon-centred
subplot(2,3,1);
xk = (o.X(1,:)-(1-mu))*lStar;  yk = o.X(2,:)*lStar;
plot(xk, yk, '-', 'LineWidth',1.1, 'Color',[0.15 0.35 0.75]); hold on; grid on; axis equal;
th = linspace(0,2*pi,200);
fill(rMoonKm*cos(th), rMoonKm*sin(th), [0.6 0.6 0.6], 'EdgeColor','none');
if isfinite(o.minAltKm)
    plot((rMoonKm+o.minAltKm)*cos(th), (rMoonKm+o.minAltKm)*sin(th), 'r--', 'LineWidth',1.0);
end
plot(xk(1),yk(1),'go','MarkerFaceColor','g','MarkerSize',6);
plot(xk(end),yk(end),'ro','MarkerFaceColor','r','MarkerSize',6);
xlabel('x - x_{Moon} [km]'); ylabel('y [km]');
title('(a) trajectory, Moon-centred');

% (b) altitude
subplot(2,3,2);
plot(tPhys*tStar/86400, altKm, '-', 'LineWidth',1.2); hold on; grid on;
if isfinite(o.minAltKm)
    yline(o.minAltKm,'r--','LineWidth',1.2);
end
set(gca,'YScale','log');
xlabel('time [days]'); ylabel('lunar altitude [km]');
title(sprintf('(b) altitude, min %.0f km', min(altKm)));

% (c) grid spacing -- in time, and in the angle the mesh actually has to resolve
subplot(2,3,3);
tmid = 0.5*(tPhys(1:end-1)+tPhys(2:end))*tStar/86400;
thMoon = unwrap(atan2(o.X(2,:), o.X(1,:)-(1-mu)));
dTheta = abs(diff(thMoon))*180/pi;
yyaxis left
plot(tmid, hPhys*tStar/3600, '-', 'LineWidth',1.4);
ylabel('\Deltat per interval [hr]');  ylim([0 1.4*max(hPhys)*tStar/3600]);
yyaxis right
semilogy(tmid, dTheta, '-', 'LineWidth',1.1);
ylabel('lunar angle swept [deg/interval]');
grid on; xlabel('time [days]');
title(sprintf(['(c) grid spacing: UNIFORM in time (ratio %.3f)\n' ...
               'but %.0f\\times non-uniform in angle swept'], ...
    max(hPhys)/min(hPhys), max(dTheta)/median(dTheta)));

% (d) continuous residual
subplot(2,3,4);
semilogy(0.5*(tPhys(1:end-1)+tPhys(2:end))*tStar/86400, Rx, '.', 'MarkerSize',6); hold on; grid on;
yline(o.maxDefect,'r--','LineWidth',1.2);
xlabel('time [days]'); ylabel('||x_{prop} - x_{node}||');
title(sprintf('(d) CONTINUOUS residual, med %.1e\nvs NLP defect %.1e (dashed)', median(Rx), o.maxDefect));

% (e) throttle
subplot(2,3,5);
stairs(tPhys*tStar/86400, o.U(4,:), '-', 'LineWidth',1.1); grid on;
ylim([-0.05 1.05]); xlabel('time [days]'); ylabel('throttle');
title(sprintf('(e) throttle, min %.6f', o.thrMin));

% (f) residual vs altitude
subplot(2,3,6);
altc = 0.5*(altKm(1:end-1)+altKm(2:end));
loglog(altc, Rx, '.', 'MarkerSize',6); grid on;
xlabel('lunar altitude [km]'); ylabel('continuous residual');
title(sprintf('(f) residual vs altitude\nworst residual at %.0f km', D.altAtWorst));

sgtitle(sprintf('DRO \\rightarrow tulip min-time  |  N=%d  t_f=%.6f  minAlt=%s km', ...
    N, o.tf, local_num(o.minAltKm)), 'FontWeight','bold');

if ~isempty(outPng)
    exportgraphics(fig, outPng, 'Resolution', 140);
    fprintf('  figure -> %s\n', outPng);
end
end

% ---------------------------------------------------------------------------
function dz = local_rhs(z, u, mu, Tmax, c)
% LOCAL_RHS  CR3BP with thrust, matching casadi_mintime_dro / tfMinEoM.
% INPUTS: z [7x1]; u [4x1]; mu; Tmax; c   OUTPUTS: dz [7x1]
r = z(1:3);  v = z(4:6);  m = z(7);
al = u(1:3);  al = al/max(norm(al),eps);  th = min(max(u(4),0),1);
dd = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
gr = [r(1) - (1-mu)*(r(1)+mu)/dd^3 - mu*(r(1)-1+mu)/rr^3;
      r(2) - (1-mu)*r(2)/dd^3      - mu*r(2)/rr^3;
           - (1-mu)*r(3)/dd^3      - mu*r(3)/rr^3];
hv = [2*v(2); -2*v(1); 0];
dz = [v; gr + hv + (th*Tmax/m)*al; -(Tmax/c)*th];
end

% ---------------------------------------------------------------------------
function s = local_num(x)
% LOCAL_NUM  'none' for NaN, else a compact number. INPUTS: x  OUTPUTS: s
if isnan(x), s = 'none'; else, s = sprintf('%.0f', x); end
end
