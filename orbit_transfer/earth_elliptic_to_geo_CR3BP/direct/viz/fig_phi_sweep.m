function fig_phi_sweep()
% FIG_PHI_SWEEP  Moon effect vs lunar phase at 10 N: the quadrupole figure.
%
% Loads the four phi0-swept certified front-door products and the same-chain
% gain=0 control, computes dmf(phi0) = m_f(CR3BP,phi0) - m_f(2-body control),
% decomposes into harmonics (mean + cos/sin(phi0) dipole + cos/sin(2 phi0)
% quadrupole -- exactly identifiable from 4 quarter-phase samples, higher
% harmonics alias; caveat printed), and renders a two-panel figure:
%   (a) dmf vs phi0 [g]: data points (green=helps, red=hurts), the fitted
%       mean+dipole+quadrupole curve, zero line;
%   (b) polar view of the same fit -- the pi-periodic tidal-quadrupole
%       signature at a glance (radius = dmf offset to keep it positive).
% Writes doc/figs/phi_sweep_dmf.png (150 dpi) for the campaign note.
%
% INPUTS:  none (paths hardcoded to the Phase-1 artifacts)
% OUTPUTS: none (figure file written; fit coefficients printed)
%
% REFERENCES:
%   [1] doc/cr3bp_geo_phase1_note.tex (lunar-phase section; consumes this fig).
%   [2] results/compare_vs_2body.md (the same numbers in table form).
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, '..', 'results');
figDir = fullfile(here, '..', '..', 'doc', 'figs');
if ~exist(figDir, 'dir'), mkdir(figDir); end

C0 = load(fullfile(resDir, 'cr3bp_T10N_2body_control.mat'), 'm_f_kg');
runs = { 0,      'cr3bp_T10N_phi0_fuel.mat';
         pi/2,   'cr3bp_T10N_phiPi2_fuel.mat';
         pi,     'cr3bp_T10N_phiPi_fuel.mat';
         3*pi/2, 'cr3bp_T10N_phi3Pi2_fuel.mat' };
nR  = size(runs, 1);
phi = zeros(nR,1);  dmfG = zeros(nR,1);
for k = 1:nR
    S = load(fullfile(resDir, runs{k,2}), 'm_f_kg');
    phi(k)  = runs{k,1};
    dmfG(k) = 1000*(S.m_f_kg - C0.m_f_kg);       % grams
end

% Harmonic decomposition on the 4 quarter-phase samples (exact for k<=2;
% harmonics >=3 ALIAS onto these -- honest caveat, needs a finer grid).
c0 = mean(dmfG);
a1 = ( dmfG(1) - dmfG(3) ) / 2;                  % cos(phi) dipole
b1 = ( dmfG(2) - dmfG(4) ) / 2;                  % sin(phi) dipole
a2 = ( dmfG(1) - dmfG(2) + dmfG(3) - dmfG(4) ) / 4;   % cos(2phi) quadrupole
fprintf('harmonics [g]: mean=%+.1f  dipole=(%+.1f cos, %+.1f sin)  quadrupole=%+.1f cos2phi\n', ...
        c0, a1, b1, a2);
phid = linspace(0, 2*pi, 361).';
fit  = c0 + a1*cos(phid) + b1*sin(phid) + a2*cos(2*phid);

fig = figure('Visible','off','Position',[60 60 1240 480]);
% (a) linear panel
ax1 = subplot(1,2,1);  hold(ax1,'on');  grid(ax1,'on');
plot(ax1, phid/pi, fit, '-', 'Color',[0.3 0.3 0.9], 'LineWidth', 1.4);
yline(ax1, 0, 'k-');
helps = dmfG >= 0;
plot(ax1, phi(helps)/pi,  dmfG(helps),  'o', 'MarkerSize', 9, ...
     'MarkerFaceColor',[0.15 0.6 0.2], 'Color',[0.1 0.4 0.15]);
plot(ax1, phi(~helps)/pi, dmfG(~helps), 'o', 'MarkerSize', 9, ...
     'MarkerFaceColor',[0.85 0.25 0.2], 'Color',[0.55 0.1 0.1]);
for k = 1:nR
    text(ax1, phi(k)/pi, dmfG(k) + sign(dmfG(k))*9, sprintf('%+.1f g', dmfG(k)), ...
         'HorizontalAlignment','center', 'FontSize', 9);
end
xlabel(ax1, '\phi_0 / \pi'); ylabel(ax1, '\Delta m_f  [g]');
xlim(ax1, [-0.08 2.08]); ylim(ax1, [-85 95]);
title(ax1, {'Moon effect vs lunar phase (10 N, c_{tf}=1.5, same-chain baseline)', ...
    sprintf('fit: %+.1f %+.1f cos\\phi_0 %+.1f sin\\phi_0 %+.1f cos2\\phi_0  [g]', c0, a1, b1, a2)});
legend(ax1, {'mean+dipole+quadrupole fit','','Moon helps','Moon hurts'}, 'Location','southeast');
% (b) polar panel (offset radius so the curve stays positive)
off = 100;
ax2 = subplot(1,2,2, polaraxes);  hold(ax2,'on');
polarplot(ax2, phid, off + fit, '-', 'Color',[0.3 0.3 0.9], 'LineWidth', 1.4);
polarplot(ax2, phid, off + 0*phid, 'k--');
polarplot(ax2, phi(helps),  off + dmfG(helps),  'o', 'MarkerSize', 8, ...
     'MarkerFaceColor',[0.15 0.6 0.2], 'Color',[0.1 0.4 0.15]);
polarplot(ax2, phi(~helps), off + dmfG(~helps), 'o', 'MarkerSize', 8, ...
     'MarkerFaceColor',[0.85 0.25 0.2], 'Color',[0.55 0.1 0.1]);
ax2.ThetaZeroLocation = 'right';
title(ax2, {'polar view: the \pi-periodic quadrupole signature', ...
            '(radius = 100 g + \Delta m_f; dashed = zero effect; \phi_0 = Moon direction at t_0)'});
outPng = fullfile(figDir, 'phi_sweep_dmf.png');
exportgraphics(fig, outPng, 'Resolution', 150);
close(fig);
fprintf('WROTE %s\n', outPng);
end
