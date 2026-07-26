% TEST_FOC_CHECK_10N  foc_check integration test against the certified 10 N
% MEE_M2 row. Recovery block copied verbatim from
% earth_elliptic_to_geo/direct/tests/test_dual_extraction.m Test 2 (returnModel
% not yet wired into refresh_duals_mee at this task, so casadi_lt_mee is called
% directly with returnModel=true, warmTight=true to get out.model.{opti,creg}).
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/tests/test_dual_extraction.m Test 2
%       (the recovery block this test copies).

thisRoot = fileparts(fileparts(mfilename('fullpath')));  % orbit_transfer/verify_common
orbitRoot = fileparts(thisRoot);                          % orbit_transfer
earthRoot = fullfile(orbitRoot, 'earth_elliptic_to_geo', 'direct');
row = fullfile(earthRoot, 'results', 'MEE_M2_10N.mat');
if ~isfile(row)
    fprintf('SKIPPED -- cache absent\n'); return;
end

cd(thisRoot); setup_verify_common;
cd(earthRoot); setup_paths;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
import casadi.*

saved = sosc_load_row(row);
p     = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
outRow = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, ...
    struct('par', p, 'mode', 'fixedtf', 'eps', 0, 'tfTarget', saved.tfTarget, ...
           'x0', saved.X(:,1), 'xf', saved.xf, 'maxIter', saved.maxIter, ...
           'warmTight', true, 'printLevel', 0, 'returnModel', true));
assert(outRow.success && outRow.maxDefect < 1e-8, ...
    ['test_foc_check_10N: warm re-solve of the 10 N row did not land ' ...
     'machine-tight (%s, defect %.3e) -- cannot judge stationarity'], ...
    outRow.ipoptStatus, outRow.maxDefect);

rep = foc_check(outRow, saved.sigma, foc_manifest('earth_mee'), struct());
assert(rep.kktStatInf < 1e-6);        assert(rep.dirTanMax < 1e-6);
assert(rep.signPct == 100);           assert(rep.lamMassEndMapped < 1e-3);
% lamMassEndMapped, not lamMassEndRel: b3c6eaf made the Hager-mapped terminal
% covector THE GATED VALUE and stopped producing the interval-dual variants,
% but did not update this line -- so this test threw 'Unrecognized field name'
% instead of gating anything, from 2026-07-25 until it was caught on 07-26.
assert(rep.singularArcNodes == 0);    assert(rep.pass);
fprintf('test_foc_check_10N: ALL PASS (Sdot min rel %.2e, nSw %d)\n', rep.sdotMinRel, rep.nSwitches);
