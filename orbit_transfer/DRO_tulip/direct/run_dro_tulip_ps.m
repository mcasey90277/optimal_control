function run_dro_tulip_ps()
%% Purpose:
%
%  This routine solves the minimum-time low-thrust transfer from a DRO orbit
%  to a Tulip orbit in the Earth-Moon CR3BP by DIRECT COLLOCATION
%  (Hermite-Simpson, CasADi + IPOPT), and compares the answer against the
%  indirect (Pontryagin shooting) solution from pumpkyn.cr3bp.tfMin.
%
%  It is the pumpkyn-style companion to run_dro_tulip.m: one straight-line
%  script, the same problem and constants as demos/lowThrustDRO2Tulip.m, and
%  the same figures at the end -- so the two methods can be read side by side.
%  The heavy machinery (mesh sweeps, certification gates, movies) lives in
%  run_dro_tulip.m; nothing here writes into pumpkyn or pumpkynPie.
%
%  References:
%  pumpkynPie/demos/lowThrustDRO2Tulip.m  (the indirect original)
%  ../doc/dro_tulip_mintime.tex           (formulation + lessons learned)

%% Constants:
  muStar = 0.012150585609624;          % Mass ratio
   lStar = 389703.264829278;           % Characteristic length (km)
   tStar = 382981.289129055;           % Characteristic time (s)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

%% Initial Orbit Conditions for DRO:
             tauND0 = 1.0;
          [~,rvND0] = pumpkynPie.cr3bp.getDRO(tauND0);
              rvND0 = pumpkyn.cr3bp.cont_np(rvND0,tauND0,muStar,1e-12);
     [tauNDv,rvNDv] = pumpkyn.cr3bp.prop(tauND0,rvND0,muStar);

%% Final Orbit Conditions:
                  Np = 7;
              tauNDf = 5*2*pi/6;
           [~,rvNDF] = pumpkyn.cr3bp.getTulip(tauNDf,Np,-1);
               rvNDF = pumpkyn.cr3bp.cont_np(rvNDF,tauNDf,muStar,1e-12);
        [tauND,rvND] = pumpkyn.cr3bp.prop(tauNDf,rvNDF,muStar);
             dvTheta = pumpkyn.util.bsxAng(rvND(:,4:6),rvNDv(1,4:6),2);
             [~,idx] = max(dvTheta);
               rvNDF = rvND(idx,:);

%% Thruster:
 Tmax = 0.07;                                    % Max thrust (N)
   m0 = 150;                                     % Initial mass (kg)
   g0 = 9.80665 * tStar^2 / (1000 * lStar);      % Gravity in ND units
 Tmax = (Tmax / m0) * tStar^2 / (lStar * 1000);  % ND thrust accel
  Isp = 900 / tStar;                             % ND specific impulse
    c = Isp * g0;                                % ND exhaust velocity

%% Solve Minimum-Time Transfer (INDIRECT), for the reference:
% The converged costates below are the demo's own 0.07 N solution.
    lambda0_tf = [ 13.9579470518969
          6.26766054139842
          7.02087277890196
        -0.969930036390389
         -5.04449989501292
         -3.11628049459521
          5.50150618633038
          4.01524259262941];

sol_lambda0_tf = pumpkyn.cr3bp.tfMin(rvND0(1:6), rvNDF(1:6), lambda0_tf, ...
                     Tmax, c, muStar);
     [tau, rv] = pumpkyn.cr3bp.tfMinProp(sol_lambda0_tf(8), ...
                     [rvND0(1:6), 1, sol_lambda0_tf(1:7)'], Tmax, c, muStar);

%% Solve Minimum-Time Transfer (DIRECT), the point of this script:
% Hermite-Simpson collocation, N = 1600 intervals, warm-started from the
% indirect trajectory. See ../doc for why fourth order is required here:
% trapezoidal collocation converges to machine-tight defects on trajectories
% that are wrong by thousands of km at periselene.
    N = 1600;
   sN = linspace(0, tau(end), N+1);
   X0 = interp1(tau, rv(:,1:7), sN, 'spline').';
   LV = interp1(tau, rv(:,11:13), sN, 'spline');
   U0 = [(-LV ./ max(vecnorm(LV,2,2),eps)).'; ones(1,N+1)];

  out = casadi_mintime_dro(rvND0(1:6), rvNDF(1:6), Tmax, c, muStar, N, ...
            X0, U0, sol_lambda0_tf(8), ...
            struct('maxIter',5000,'scheme','hermite-simpson'));

%% Compare the Two Methods:
   tfInd = sol_lambda0_tf(8);
   tfDir = out.tf;
    rvI  = interp1(tau/tau(end), rv(:,1:6), out.s(:), 'spline');   % matched s
  posErr = vecnorm(out.X(1:3,:).' - rvI(:,1:3), 2, 2) .* lStar;    % km
  velErr = vecnorm(out.X(4:6,:).' - rvI(:,4:6), 2, 2) .* lStar ./ tStar .* 1000;

fprintf('\n Minimum-Time DRO -> Tulip, both methods \n');
fprintf('--------------------------------------------------\n');
fprintf(' t_f indirect  : %.7f ND = %.3f days\n', tfInd, tfInd*tStar/86400);
fprintf(' t_f direct    : %.7f ND = %.3f days\n', tfDir, tfDir*tStar/86400);
fprintf(' difference    : %.2e ND = %.2f s\n', abs(tfDir-tfInd), abs(tfDir-tfInd)*tStar);
fprintf(' position match: max %.3f km, median %.3f km (phasing removed)\n', ...
    max(posErr), median(posErr));
fprintf(' velocity match: max %.4f m/s\n', max(velErr));
fprintf(' final mass    : %.6f (indirect %.6f)\n', out.mf, rv(end,7));
fprintf(' NLP defect    : %.1e   min throttle: %.6f\n', out.maxDefect, out.thrMin);

%% Propellant Usage and dV Analysis:
          mTot = out.X(7,:) .* m0;                       % Mass vs. time (kg)
         dVtot = (c.*lStar./tStar).*log(m0 ./ mTot);     % dV (km/s)
fprintf(' propellant    : %.3f kg    dV: %.4f km/s\n\n', m0-mTot(end), dVtot(end));

%% Show 3D Trajectory in CR3BP:
hIn = figure('color', [0 0 0]);
pumpkyn.cr3bp.showMoon(lStar, muStar, hIn);
plot3(out.X(1,:), out.X(2,:), out.X(3,:), 'w');
plot3(rv(:,1), rv(:,2), rv(:,3), ':', 'color', [1 0.5 0]);
plot3(rvND0(1), rvND0(2), rvND0(3), '.g', 'markersize', 12);
plot3(rvNDF(1), rvNDF(2), rvNDF(3), '.r', 'markersize', 12);
legend({'','direct (white)','indirect (orange)'}, 'TextColor', 'w');
axis equal off;
set(gca, 'color', 'k', 'clipping', 'off');

%% Show the Agreement:
figure('Color', [1 1 1]);
yyaxis left
semilogy(out.s.*tfDir.*tStar./3600, max(posErr,1e-9), 'k', 'linewidth', 1.5);
grid on; ylabel('position difference [km]');
set(gca, 'YColor', 'k')
yyaxis right
semilogy(out.s.*tfDir.*tStar./3600, ...
    (vecnorm(out.X(1:3,:) - [1-muStar;0;0], 2, 1).*lStar - 1737.4), ...
    'linewidth', 1.5);
ylabel('lunar altitude [km]'); xlabel('Time [Hrs]');
title(sprintf('Direct vs indirect: max %.2f km over %.1f days', ...
    max(posErr), tfDir*tStar/86400));

end
