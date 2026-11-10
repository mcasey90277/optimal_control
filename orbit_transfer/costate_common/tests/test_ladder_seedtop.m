function ok = test_ladder_seedtop()
% TEST_LADDER_SEEDTOP  opts.seedTop / opts.redo of thrust_ladder_library.
%
% Four checks on the stage 1 GTO -> DRO demo cell (orient 0, DRO tau 1,
% departure phase 0.10, arrival phase 0, Earth floor 200 km), the one cell
% whose two top-rung basins are MEASURED (GTO_DRO/README.md, stage 1):
%   (a) a bad seedTop mesh refuses before any solve
%   (b) redo without cells refuses
%   (c) a normal call solves rungs [15 12]; the 15 N cold start lands near
%       0.3040 ND and the 12 N rung near 0.3070 ND
%   (d) redo + seedTop built from the 12 N entry re-solves the cell and the
%       15 N rung lands on the faster basin (t_f < 0.29 ND), OK, with the
%       attempt counter reset to 1. The call raises directGateKm to 1000
%       (the direct control of that basin flies 193.6 km off, so the
%       engine's first flown gate would stop it before refinement) and
%       leaves gateKm at the campaign 100 km; the CALLER enforces gate and
%       acceptance on the refined entries, and meta.gateKm stays 100
%   Also: NaN seedTop refuses; redo on a missing sheet refuses (redoSheet)
% The option-off identity is covered by the existing bitwise regression
% DRO_tulip/direct/certify/tests/test_minenergy_objective (no seedTop,
% no redo -> same code path as before this change).
%
% INPUTS:  none
% OUTPUTS: ok - logical
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    od = cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(od);
end
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

outDir = fullfile(tempdir, 'test_ladder_seedtop');
if isfolder(outDir), rmdir(outDir, 's'); end
mkdir(outDir);
outMat = fullfile(outDir, 'cell.mat');

N    = 400;
base = struct('rungs', [15 12], 'ispS', 1710, 'm0kg', 150, 'N', N, ...
    'floorKm', 500, 'floorEarthKm', 200, 'maxIter', 6000, 'maxCpuSec', 600, ...
    'gateKm', 100, 'accTol', 1e-6, 'nD', 1, 'nA', 1, 'sD0', 0.10, 'sA0', 0, ...
    'depFamily', 'gto', 'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'dro', 'arrParams', struct('tau', 1), 'tauDRO', 0, ...
    'thrLock', true, 'tf0', 0.30, 'resume', false, ...
    'logFile', fullfile(outDir, 'cell.log'));

% (a) bad seedTop mesh refuses before any solve
o = base;  o.cells = [1 1];  o.redo = true;
o.seedTop = struct('X', zeros(7, N), 'U', zeros(4, N+1), 'tf', 0.3);
try
    thrust_ladder_library(outMat, o);
    error('test:noThrow', 'bad seedTop mesh did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:seedTop'), ...
        'wrong error for a bad seedTop: %s', err.identifier);
end
o.seedTop = struct('X', nan(7, N+1), 'U', zeros(4, N+1), 'tf', 0.3);
try
    thrust_ladder_library(outMat, o);
    error('test:noThrow', 'NaN seedTop did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:seedTop'), ...
        'wrong error for a NaN seedTop: %s', err.identifier);
end
assert(~isfile(outMat), 'a refused call must not write the sheet');

% (b) redo without cells refuses
o = base;  o.redo = true;
try
    thrust_ladder_library(outMat, o);
    error('test:noThrow', 'redo without cells did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:redoCells'), ...
        'wrong error for redo without cells: %s', err.identifier);
end

% (b2) redo + cells on a missing sheet refuses (nothing to redo into)
o = base;  o.cells = [1 1];  o.redo = true;  o.resume = true;
try
    thrust_ladder_library(outMat, o);
    error('test:noThrow', 'redo on a missing sheet did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:redoSheet'), ...
        'wrong error for redo on a missing sheet: %s', err.identifier);
end
assert(~isfile(outMat), 'a refused call must not write the sheet');

% (c) the normal ladder: 15 N cold start, 12 N warm from it
P = thrust_ladder_library(outMat, base);
assert(all(P.OK(1, 1, :)), 'both rungs must be OK on the demo cell');
tf15 = P.TF(1, 1, 1);  tf12 = P.TF(1, 1, 2);
fprintf('  cold ladder: tf(15 N) = %.5f, tf(12 N) = %.5f ND\n', tf15, tf12);
assert(abs(tf15 - 0.3040) < 0.005, '15 N cold start expected near 0.3040 ND, got %.5f', tf15);
assert(abs(tf12 - 0.3070) < 0.005, '12 N rung expected near 0.3070 ND, got %.5f', tf12);

% (d) redo + seedTop from the 12 N entry: the reverse warm start
ndp  = nd_propulsion(12, 1710, 150, P.meta.lStar, P.meta.tStar);
[tD, rvD] = get_family_orbit('gto', struct('orientDeg', 0));
rv0  = phase_state(tD, rvD);  rv0 = rv0(0.10);
z12  = P.Z8(:, 1, 1, 2);
[tt, yy] = pumpkyn.cr3bp.tfMinProp(z12(8), [rv0(:); 1; z12(1:7)], ndp.Tnd, ndp.cnd, P.meta.muStar);
tg   = linspace(0, z12(8), N+1);
Y    = interp1(tt, yy, tg, 'spline');                         % [N+1 x 14]
lamV = Y(:, 11:13).';
U    = [-lamV ./ max(vecnorm(lamV, 2, 1), eps); ones(1, N+1)];
seedTop = struct('X', Y(:, 1:7).', 'U', U, 'tf', z12(8)*(12/15)^0.6);
o = base;  o.resume = true;  o.cells = [1 1];  o.redo = true;  o.seedTop = seedTop;
% Measured: the direct control of the faster basin flies 193.6 km off at
% N = 400, so the engine's first flown gate (100 km) stops it before
% ms_tfmin runs. Refinement closes it (flown 7e-5 km, acceptance 5e-13).
% Only that first gate is raised; gateKm stays 100 and the caller enforces
% it on the refined entries (below).
o.directGateKm = 1000;
P2 = thrust_ladder_library(outMat, o);
tf15b = P2.TF(1, 1, 1);
fprintf('  redo from the 12 N seed: tf(15 N) = %.5f ND, OK = %d, ATT = %d\n', ...
    tf15b, P2.OK(1, 1, 1), P2.ATT(1, 1));
assert(P2.OK(1, 1, 1), 'the re-solved 15 N rung must be OK');
assert(tf15b < 0.29, 'the reverse warm start must land on the faster basin (< 0.29 ND), got %.5f', tf15b);
assert(P2.ATT(1, 1) == 1, 'redo must reset the attempt counter (got %d)', P2.ATT(1, 1));
assert(P2.FLYKM(1, 1, 1) < 100, 'refined 15 N flown miss %.3g km exceeds the campaign gate', P2.FLYKM(1, 1, 1));
assert(P2.ACCDZ(1, 1, 1) < 1e-6, 'refined 15 N acceptance %.3g too large', P2.ACCDZ(1, 1, 1));
assert(P2.ACCDZ(1, 1, 2) < 1e-6, 'refined 12 N acceptance %.3g too large', P2.ACCDZ(1, 1, 2));
assert(P2.OK(1, 1, 2) && P2.FLYKM(1, 1, 2) < 100, 'the 12 N rung must re-walk OK within the gate');
Q = load(outMat, 'TF', 'OK', 'meta');
assert(abs(Q.TF(1, 1, 1) - tf15b) < 1e-12 && Q.OK(1, 1, 1), 'the sheet on disk must carry the re-solve');

assert(Q.meta.gateKm == 100, 'meta.gateKm must stay at the campaign 100 km (got %g)', Q.meta.gateKm);

rmdir(outDir, 's');
ok = true;
fprintf('test_ladder_seedtop: PASS\n');
end
