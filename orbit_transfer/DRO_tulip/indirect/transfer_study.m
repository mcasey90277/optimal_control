%% TRANSFER_STUDY  One DRO -> tulip minimum-time transfer, end to end.
%
%   Edit the PARAMETERS block, press Run. Five sections:
%     1  parameters                 everything you choose, in one place
%     2  solve                      costates, trajectory, mass, Delta-V
%     3  NECESSARY conditions       Pontryagin, first order
%     4  SUFFICIENCY hypotheses     Bonnard-Caillau-Trelat
%     5  interactive 3D plot        rotate it with the mouse
%
%   Sections 3 and 4 are separate on purpose. The first-order conditions are
%   HYPOTHESES of the sufficiency theorem, not consequences of it, and the
%   conjugate test is computed ALONG the extremal -- off one its determinant
%   means nothing. On this very problem 12 of 14 candidates satisfied every
%   first-order condition and were then refuted by the second-order test.
%
%   Everything here runs through the same gate stack as the catalog, so a
%   number is never reported without its verdict.
%
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.

%% ------------------------------------------------------------------------
%  1. PARAMETERS -- the only block you need to edit
%% ------------------------------------------------------------------------
sD        = 0;          % departure phase, fraction of the DRO period   [0, 1)
sA        = 0.0754;     % arrival phase,  fraction of the tulip period  [0, 1)

thrustN   = 0.070;      % thrust [N]
ispS      = 900;        % specific impulse [s]
m0kg      = 150;        % initial mass [kg]

tauDRO    = 1.0;        % departure DRO period [ND]
NpTulip   = 7;          % tulip petal count

wallSec   = 1800;       % budget for a walk, if the pair is not in the library
outPng    = '';         % e.g. 'results/my_transfer.png' to save the figure

%% ------------------------------------------------------------------------
%  2. SOLVE
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));

popt = struct('thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
              'tauDRO', tauDRO, 'NpTulip', NpTulip, 'wallSec', wallSec);
T = run_dro_tulip(sD, sA, popt);

if ~T.ok
    fprintf(2, '\nNo certified transfer at (%.4f, %.4f): %s\n\n', sD, sA, T.reason);
    return
end

% the setup, for the flight and the plot below
[B, ~] = arclength_arrival('setup', setfield(popt, 'sD', sD)); %#ok<SFLD>
z8 = T.z;                      % [lambda_0(7); t_f]  -- the costates
Yj = T.Y;                      % [14 x K] multiple-shooting junctions
[tu, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [B.stateD(sD); 1; z8(1:7)], B.Tnd, B.cnd, B.mu);

fprintf('\nSolved from: %s\n', T.source);
fprintf('  t_f          %.6f ND   (%.4f days)\n', z8(8), T.tfDays);
fprintf('  Delta-V      %.6f km/s\n', T.dvKms);
fprintf('  propellant   %.3f kg of %g\n', T.mfKg, m0kg);
fprintf('  lambda_0     [%s]\n', strjoin(compose('%+.6g', z8(1:7)'), ' '));
fprintf('  trajectory   %d propagator samples, %d junctions\n', numel(tu), size(Yj, 2));

%% ------------------------------------------------------------------------
%  3-4. NECESSARY, then SUFFICIENCY
%% ------------------------------------------------------------------------
R = report_optimality(T);

%% ------------------------------------------------------------------------
%  5. INTERACTIVE 3D PLOT  (drag to rotate)
%% ------------------------------------------------------------------------
P = plot_transfer_3d(T, B, struct('outPng', outPng));
fprintf('Figure %d is rotatable: drag to spin, scroll to zoom.\n', P.fig.Number);
