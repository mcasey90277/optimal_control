function smoke_gto_cell()
%% Purpose:
%
%   SMOKE TEST for the GTO -> tulip costate-catalog campaign: a single
%   (departure phase x arrival phase) cell at a single high-thrust rung
%   (15 N), run through the SAME production engine the full campaign will
%   use (DRO_tulip/indirect/thrust_ladder_library, family-agnostic via
%   costate_common/get_family_orbit). Proves the cold-start path -- direct
%   solve -> preflight screen -> certify -> ms_tfmin refine -> pumpkyn
%   tfMin acceptance -- before any sheet sweep is attempted.
%
%   Departure: GTO pseudo-family locus at orientDeg=0, sampled at PERIGEE
%   (sD0=0.0 -- time-since-perigee fraction). The spec's first choice,
%   apogee (sD0=0.5), was tried first and FAILED the flown-arrival gate by
%   five orders of magnitude (11.95e6 km, vs 100 km) -- the direct
%   transcription converged (defect ok, preflight ok) but to a solution
%   that never actually reaches the tulip; likely a bad local min for a
%   15 N cold start from that phase, not a scale/unit problem (the GTO
%   locus itself checks out physically: apogee radius and velocity both
%   land within 1% of two-body vis-viva). Perigee, tried per the brief's
%   one-alternative rule, converges cleanly. See task-2-report.md for the
%   apogee log.
%   Arrival: tulip Np=7, pm=-1 (the reference tulip used across every other
%   costate campaign in this program).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes results/smoke_gto_cell.mat + results/smoke.log)
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

setup_paths();

here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~isfolder(resDir), mkdir(resDir); end
outMat  = fullfile(resDir, 'smoke_gto_cell.mat');
logFile = fullfile(resDir, 'smoke.log');
if isfile(logFile), delete(logFile); end

opts = struct( ...
    'rungs',     15, ...
    'ispS',      1710, ...
    'm0kg',      150, ...
    'nD',        1, ...
    'nA',        1, ...
    'sD0',       0.0, ...             % perigee (apogee failed the gate)
    'sA0',       0, ...
    'depFamily', 'gto', ...
    'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'tulip', ...
    'arrParams', struct('Np', 7, 'pm', -1), ...
    'tauDRO',    0, ...                % sheet-key rule: departure-period
    ...                                % key carries orientDeg, not a
    ...                                % Kepler period, for the GTO family
    'logFile',   logFile);

fprintf('\n=== GTO->TULIP SMOKE: 1x1 cell, 15 N ===\n');
thrust_ladder_library(outMat, opts);

Q = load(outMat);
fprintf('\nSmoke result: nnz(OK) = %d (expect 1)\n', nnz(Q.OK));
if nnz(Q.OK) == 1
    [iD, iA, kr] = ind2sub(size(Q.OK), find(Q.OK, 1));
    fprintf('  tf=%.5f ND  fly=%.2f km  acc(dz)=%.1e  RES=%.1e\n', ...
        Q.TF(iD,iA,kr), Q.FLYKM(iD,iA,kr), Q.ACCDZ(iD,iA,kr), Q.RES(iD,iA,kr));
end
end
