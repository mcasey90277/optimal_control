function out = NLP_lowThrust_GTO_Tulip(N, guessMode, makePlot)
% NLP_LOWTHRUST_GTO_TULIP  Direct-NLP twin of the pumpkynPie indirect demo.
%
% Solves the same minimum-time low-thrust GTO -> tulip transfer as
% pumpkynPie's lowThrust_GTO_Tulip.m, but by DIRECT TRANSCRIPTION: states
% and controls on an N-segment trapezoidal mesh, minimum-tf objective, and
% fmincon (interior-point with analytic sparse gradients) instead of
% costate shooting. Produces the same outputs: transfer time, propellant
% and delta-V accounting, and the two-panel trajectory figure.
%
% INPUTS:
%   N         - (optional) trapezoidal segments [scalar, default 12000;
%               mesh-refinement convergence of tf vs the indirect 6.290694:
%               N=3000 -> ~6.2658 (stall scatter ~1e-3), 6000 ->
%               6.281802, 12000 -> 6.288574]
%   guessMode - (optional) 'indirect' (default) or 'tangential', see
%               BUILD_GUESS
%   makePlot  - (optional) true (default) to draw the two-panel figure
%
% OUTPUTS:
%   out - results struct:
%           .tf_ND, .tf_days   transfer time
%           .mProp_kg, .dV_kms propellant used, total delta-V
%           .maxDefect         inf-norm of mesh defects at the solution
%           .devPos_km, .devVel_kms  max node-wise deviation from the
%                              indirect reference arc (open-loop replay is
%                              NOT a usable check at this sensitivity)
%           .X, .W, .tauMesh   mesh solution
%           .exitflag          fmincon exit flag
%
% REFERENCES:
%   [1] pumpkynPie Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m (indirect
%       counterpart whose outputs this reproduces).
%   [2] Betts, "Practical Methods for Optimal Control and Estimation Using
%       Nonlinear Programming," 2nd ed., SIAM, 2010.

if nargin < 1 || isempty(N),         N = 12000;             end
if nargin < 2 || isempty(guessMode), guessMode = 'indirect'; end
if nargin < 3 || isempty(makePlot),  makePlot = true;        end

setup_paths();

% --- Earth-Moon CR3BP characteristic quantities ---------------------------
M      = 5.9736E24 + 7.35E22;      % combined primary mass (kg)
G      = 6.67384e-20;              % grav constant (km^3/kg/s^2)
muStar = 0.012150585609624;        % mass ratio
lStar  = 389703.264829278;         % characteristic length (km)
tStar  = 382981.289129055;         % characteristic time (s)

% --- GTO departure state (dimensional -> CR3BP ND) -------------------------
muEarth = G*(1 - muStar)*M;
rEarth  = 6378;
rPer    = rEarth + 350;
rApo    = rEarth + 35786;
sma     = (rApo + rPer)/2;
ecc     = (rApo - rPer)/(rApo + rPer);
oev0    = [sma, ecc, 0, -25*pi/180, 0, 0];
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, oev0, 2);
rv0ND   = pumpkyn.cr3bp.fromPCI(0, [r0, v0], muStar, tStar, lStar, 1);

% --- Tulip arrival state ---------------------------------------------------
Np   = 7;  pm = -1;  tau0 = (5/6)*2*pi;
[~, x0Tulip] = pumpkyn.cr3bp.getTulip(tau0, Np, pm, 1e-12);
[~, rvTgt] = pumpkyn.cr3bp.prop(tau0, x0Tulip, muStar);
[~, idx_f] = max(rvTgt(:,5));
rv0 = rv0ND;          % departure = GTO state at epoch (demo's idx_0 = 1)
rvf = rvTgt(idx_f, :);

% --- Spacecraft (15 kg, 25 mN, Isp 2100 s) in ND units ---------------------
m0     = 15;
g0     = 9.80665*tStar^2/(1000*lStar);
Tmax   = (0.025/m0)*tStar^2/(lStar*1000);
c      = (2100/tStar)*g0;

% --- Build guess and solve -------------------------------------------------
fprintf('Building %s guess on N = %d segments (density-matched mesh)...\n', ...
        guessMode, N);
[Z0, sigma] = build_guess(guessMode, N, rv0, rvf, Tmax, c, muStar);

[~, ceq0] = nlp_constraints(Z0, sigma, Tmax, c, muStar);
fprintf('warm-start max defect = %.3g\n', max(abs(ceq0(1:7*N))));

tSolve = tic;
[~, nlp] = solve_tfmin_nlp(Z0, sigma, rv0, rvf, Tmax, c, muStar);
tSolve = toc(tSolve);

tf      = nlp.tf;
tauMesh = sigma.'.*tf;
mf      = nlp.X(7, end);
mProp   = m0*(1 - mf);
dVtot   = c*log(1/mf)*lStar/tStar;

% --- Validation: node-wise deviation from the indirect reference arc ------
% (An open-loop replay is NOT a usable check here: ~40 perigee passes
% amplify control-interpolation error by ~1e6, so replay diverges even for
% a perfect control -- the same sensitivity that makes shooting hard.)
[X0g, ~, ~] = unpack_z(Z0, N);
devPos = max(sqrt(sum((nlp.X(1:3,:) - X0g(1:3,:)).^2, 1)))*lStar;        % km
devVel = max(sqrt(sum((nlp.X(4:6,:) - X0g(4:6,:)).^2, 1)))*lStar/tStar;  % km/s

fprintf('\n=== NLP solution (fmincon flag %d, %.1f min) ===\n', ...
        nlp.exitflag, tSolve/60);
fprintf('tf         = %.6f ND = %.4f days\n', tf, tf*tStar/86400);
fprintf('max defect = %.3g\n', nlp.maxDefect);
fprintf('propellant = %.4f kg of %g kg, dV = %.4f km/s\n', mProp, m0, dVtot);
fprintf('max node deviation from indirect arc: %.3g km, %.3g km/s\n', ...
        devPos, devVel);

out = struct('tf_ND', tf, 'tf_days', tf*tStar/86400, ...
             'mProp_kg', mProp, 'dV_kms', dVtot, ...
             'maxDefect', nlp.maxDefect, ...
             'devPos_km', devPos, 'devVel_kms', devVel, ...
             'X', nlp.X, 'W', nlp.W, 'tauMesh', tauMesh, ...
             'exitflag', nlp.exitflag);

% --- Two-panel visualization (top-down + edge-on), demo style --------------
if makePlot
    [~, rvCoast] = pumpkyn.cr3bp.prop(tau0, nlp.X(1:6, end).', muStar);
    rvAll = [nlp.X(1:3, :).'; rvCoast(2:end, 1:3)];
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
