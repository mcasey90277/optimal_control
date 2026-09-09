function ok = test_plot_phase_sheet()
%% Purpose:
%
%   Tests plot_phase_sheet -- the (departure x arrival) phase-torus figure.
%   A plot test cannot judge whether a picture is good, so it checks the
%   things that make one WRONG: that the grid is drawn with departure on one
%   axis and arrival on the other at the right sizes, that uncertified cells
%   are drawn as gaps rather than as zeros or as interpolated colour, that
%   the reported extremes match the data, and that it writes the file it says
%   it wrote.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
nD = 12;  nA = 12;

Q = struct();
Q.sD = mod((0:nD-1)/nD, 1);
Q.sA = mod(0.0754 + (0:nA-1)/nA, 1);
Q.OK = false(nD, nA);  Q.TF = nan(nD, nA);
Q.TF(1, 1) = 4.0;  Q.OK(1, 1) = true;          % fastest
Q.TF(1, 3) = 5.5;  Q.OK(1, 3) = true;
Q.TF(4, 1) = 6.0;  Q.OK(4, 1) = true;          % slowest
Q.rungs = 0.070;
Q.meta = struct('tStar', 382981.289129055, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1, 'NpTulip', 7);

png = fullfile(tempdir, 'test_phase_sheet.png');
if isfile(png), delete(png); end
I = plot_phase_sheet(Q, png);

ok = chk(ok, isfile(png), sprintf('figure written to %s', png));
ok = chk(ok, I.nCert == 3 && I.nCells == nD*nA, ...
         sprintf('counts: %d certified of %d cells', I.nCert, I.nCells));
ok = chk(ok, abs(I.tfMinDays - 4.0*Q.meta.tStar/86400) < 1e-9 && ...
             abs(I.tfMaxDays - 6.0*Q.meta.tStar/86400) < 1e-9, ...
         sprintf('extremes from the data: %.4f .. %.4f d', I.tfMinDays, I.tfMaxDays));
ok = chk(ok, isequal(I.gridSize, [nD nA]), sprintf('grid drawn %s (departure x arrival)', mat2str(I.gridSize)));
ok = chk(ok, all(isnan(I.C(~Q.OK))) && ~any(isnan(I.C(Q.OK))), ...
         'uncertified cells are gaps (NaN), not zeros');
close all
if ok, fprintf('TEST_PLOT_PHASE_SHEET: ALL PASS\n'); else, fprintf('TEST_PLOT_PHASE_SHEET: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
