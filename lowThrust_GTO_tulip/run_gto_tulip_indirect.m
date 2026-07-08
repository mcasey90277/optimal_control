function out = run_gto_tulip_indirect(makePlot)
% RUN_GTO_TULIP_INDIRECT  Min-time GTO -> tulip transfer, indirect method.
%
% Reproduces the pumpkynPie demo lowThrust_GTO_Tulip.m end-to-end using the
% local shooting machinery (LT_PMP_EOM + SHOOT_RESIDUAL_TF +
% SOLVE_TFMIN_INDIRECT). Endpoint construction (GTO state, tulip seed)
% leverages pumpkyn; the PMP solver is the tutorial's own.
%
% INPUTS:
%   makePlot - (optional) true to draw the two-panel trajectory figure
%              [scalar logical, default false]
%
% OUTPUTS:
%   out - results struct:
%           .zSol      converged [lambda0(7); tf] [8x1]
%           .resNorm   terminal residual 2-norm [scalar]
%           .tf_days   transfer time (days) [scalar]
%           .mProp_kg  propellant used (kg) [scalar]
%           .dV_kms    total delta-V (km/s) [scalar]
%           .tau, .rv  trajectory time + augmented states [Nx1], [Nx14]
%
% REFERENCES:
%   [1] pumpkynPie Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m

if nargin < 1, makePlot = false; end

% --- Earth-Moon CR3BP characteristic quantities -------------------------
M      = 5.9736E24 + 7.35E22;      % combined primary mass (kg)
G      = 6.67384e-20;              % grav constant (km^3/kg/s^2)
muStar = 0.012150585609624;        % mass ratio
lStar  = 389703.264829278;         % characteristic length (km)
tStar  = 382981.289129055;         % characteristic time (s)

% --- GTO departure state (dimensional -> CR3BP ND) -----------------------
muEarth = G*(1 - muStar)*M;
rEarth  = 6378;                                  % km
rPer    = rEarth + 350;                          % perigee radius (km)
rApo    = rEarth + 35786;                        % apogee radius (km)
sma     = (rApo + rPer)/2;
ecc     = (rApo - rPer)/(rApo + rPer);
oev0    = [sma, ecc, 0, -25*pi/180, 0, 0];       % [a e i argp raan nu]
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, oev0, 2);
rv0ND   = pumpkyn.cr3bp.fromPCI(0, [r0, v0], muStar, tStar, lStar, 1);

% --- Tulip arrival state -------------------------------------------------
Np   = 7;                    % petals
pm   = -1;                   % southern hemisphere
tau0 = (5/6)*2*pi;           % tulip period (ND)
[~, x0Tulip] = pumpkyn.cr3bp.getTulip(tau0, Np, pm, 1e-12);

[~, rvTgt] = pumpkyn.cr3bp.prop(tau0, x0Tulip, muStar);
[~, idx_f] = max(rvTgt(:,5));                    % arrival: max y-velocity
rv0 = rv0ND;          % departure = GTO state at epoch (demo's idx_0 = 1)
rvf = rvTgt(idx_f, :);

% --- Spacecraft (15 kg, 25 mN, Isp 2100 s) in ND units -------------------
m0     = 15;                                     % kg
g0     = 9.80665*tStar^2/(1000*lStar);
Tmax_N = 0.025;                                  % N
Tmax   = (Tmax_N/m0)*tStar^2/(lStar*1000);       % ND accel at m = 1
Isp    = 2100/tStar;                             % ND
c      = Isp*g0;                                 % ND exhaust velocity

% --- Converged reference guess (from the pumpkynPie demo) ----------------
zGuess = [ 190.476497248065
           -79.7064866984696
            -0.430399154713168
             0.301159446575878
             0.586671892449694
            -0.00711582435720301
             4.32931089137559
             6.29081541876621];

[zSol, resNorm] = solve_tfmin_indirect(rv0, rvf, zGuess, Tmax, c, muStar);

% --- Propagate the converged trajectory ----------------------------------
y0   = [rv0(:); 1; zSol(1:7)];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[tau, rv] = ode113(@lt_pmp_eom, [0 zSol(8)], y0, opts, Tmax, c, muStar);

mProp = m0*(1 - rv(end,7));
dVtot = c*log(1/rv(end,7))*lStar/tStar;

out = struct('zSol', zSol, 'resNorm', resNorm, ...
             'tf_days', zSol(8)*tStar/86400, ...
             'mProp_kg', mProp, 'dV_kms', dVtot, ...
             'tau', tau, 'rv', rv);

fprintf('tf       = %.6f ND = %.4f days\n', zSol(8), out.tf_days);
fprintf('||R||    = %.3g\n', resNorm);
fprintf('propellant = %.4f kg of %g kg, dV = %.4f km/s\n', mProp, m0, dVtot);

% --- Two-panel visualization (top-down + edge-on), demo style ------------
if makePlot
    [tauCoast, rvCoast] = pumpkyn.cr3bp.prop(tau0, rv(end,1:6), muStar); %#ok<ASGLU>
    rvAll = [rv(:,1:3); rvCoast(2:end,1:3)];
    Lpts  = pumpkyn.cr3bp.lagrangePts(muStar);

    fig = figure('Color', [1 1 1]);
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    viewAngles = {[0 90], [90 0]};
    for panel = 1:2
        ax = nexttile(tl, panel);
        plot3(ax, Lpts(1,1), Lpts(1,2), Lpts(1,3), '.k', 'MarkerSize', 15);
        hold(ax, 'on');
        text(ax, Lpts(1,1), Lpts(1,2), Lpts(1,3), '\leftarrow L_1', ...
             'Color', 'k', 'Rotation', 45);
        plot3(ax, rvAll(:,1), rvAll(:,2), rvAll(:,3), 'k', 'LineWidth', 1);
        axis(ax, 'equal'); grid(ax, 'on'); ax.Clipping = 'off';
        xlabel(ax, 'X [ND]'); ylabel(ax, 'Y [ND]'); zlabel(ax, 'Z [ND]');
        view(ax, viewAngles{panel}(1), viewAngles{panel}(2));
        if panel == 2, zlim(ax, [-0.15 0.15]); end
    end
end
end
