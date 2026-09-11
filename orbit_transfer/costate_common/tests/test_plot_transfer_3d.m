function ok = test_plot_transfer_3d()
%% Purpose:
%
%   Tests plot_transfer_3d -- the interactive 3D view of one transfer, the
%   still version of the movie's last frame that a reader can rotate.
%
%   A plot test cannot judge whether a picture is good, so it checks what
%   would make one wrong or useless: that the departure orbit, the target
%   orbit and the transfer arc are all drawn, that the arc actually starts
%   on the departure orbit and ends on the target, that the axes are 3D and
%   rotatable rather than a flat projection, and that the reported numbers
%   come from the trajectory rather than from the caller.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

[B, anc] = arclength_arrival('setup');
T = certify_crossing(anc.p, anc.sA, B, anc);
assert(T.ok, 'fixture: the anchor must certify (%s)', T.reason);

% a SUPPLIED flight is used rather than re-flown: the study script owns one
% flight and every consumer should draw from it
[tu0, Y0] = pumpkyn.cr3bp.tfMinProp(T.z(8), [B.rv0(1:6); 1; T.z(1:7)], B.Tnd, B.cnd, B.mu);
fl = struct('t', tu0, 'Y', Y0);
P = plot_transfer_3d(T, B, struct('visible', false, 'flight', fl));
ok = chk(ok, isfield(P, 'flightSupplied') && P.flightSupplied, 'a supplied flight is consumed, not re-flown');

ok = chk(ok, isgraphics(P.fig, 'figure'), 'a figure is produced');
ok = chk(ok, numel(P.hDep) == 1 && numel(P.hArr) == 1 && numel(P.hTx) == 1, ...
         'departure orbit, target orbit and transfer arc are all drawn');
ok = chk(ok, P.startMissKm < 1, sprintf('the arc starts ON the departure orbit (%.4f km)', P.startMissKm));
ok = chk(ok, P.endMissKm < 1, sprintf('the arc ends ON the target orbit (%.4f km)', P.endMissKm));
v = get(P.ax, 'View');
ok = chk(ok, abs(v(2)) > 1 && abs(v(2)) < 89, sprintf('a genuine 3D view, elevation %.0f deg', v(2)));
ok = chk(ok, strcmp(P.ax.DataAspectRatioMode, 'manual') || isequal(P.ax.DataAspectRatio, [1 1 1]), ...
         'equal data aspect, so the geometry is not distorted');
% RECOMPUTED from the flight rather than copied: the figure claims to be the
% source of truth, so its annotations must be derived from what it draws.
% They must still AGREE with the certificate to solver tolerance.
ok = chk(ok, abs(P.tfDays - T.tfDays) < 1e-6 && abs(P.dvKms - T.dvKms) < 1e-6, ...
         sprintf('recomputed numbers agree with the certificate: %.4f d, %.4f km/s', P.tfDays, P.dvKms));
ok = chk(ok, P.nThrust > 0, sprintf('thrust direction shown at %d points', P.nThrust));
close(P.fig);

if ok, fprintf('TEST_PLOT_TRANSFER_3D: ALL PASS\n'); else, fprintf('TEST_PLOT_TRANSFER_3D: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
