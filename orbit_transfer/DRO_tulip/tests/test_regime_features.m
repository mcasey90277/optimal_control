function ok = test_regime_features()
%% Purpose:
%
%   Tests regime_features -- the per-arm feature extractor of the REGIME MAP
%   (goal 2026-09-07: when does each continuation family work, and why).
%   Fixtures are REAL arms already on disk from the high-gamma race at
%   (2,5) gamma = 2.0, the case where only huberc reached the floor:
%
%     huberc  p -> 0.001, delta -> 0.0039, m_f 0.948585   (arrived)
%     eps     walled at p = 0.657                          (walled)
%     huber   walled at p = 0.315                          (walled)
%
%   The extractor must classify those three outcomes, put every family on
%   the SAME sharpness axis (effective ramp width in Q: eps 2p, huber 0,
%   huberc delta), and report the switch-structure and geometry features the
%   regime hypotheses H1-H3 are stated in terms of.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
Lh = load(fullfile(resDir, 'minenergy_highgamma.mat'));
rec = Lh.R(find(arrayfun(@(r) isequal([r.iD r.iA], [2 5]) && abs(r.gam - 2) < 1e-9, Lh.R), 1));

% (1) the ARRIVED arm -----------------------------------------------------
A = armOf(fullfile(resDir, 'minfuel_hg_c25_g200_huberc.mat'));
F = regime_features(A, rec, 'huberc');
ok = chk(ok, strcmp(F.outcome, 'floor'), sprintf('huberc outcome = %s (want floor)', F.outcome));
ok = chk(ok, abs(F.pFloor - 0.001) < 1e-9, sprintf('huberc p_floor = %.4g', F.pFloor));
ok = chk(ok, abs(F.rampWidth - F.deltaFloor) < 1e-12 && F.rampWidth > 0, ...
         sprintf('huberc ramp width = delta = %.4g', F.rampWidth));
ok = chk(ok, abs(F.mf - 0.948585) < 1e-5, sprintf('huberc m_f = %.6f', F.mf));
ok = chk(ok, F.nCross >= 1, sprintf('huberc switches at the floor: nCross = %d', F.nCross));
ok = chk(ok, F.revs > 0.5 && F.revs < 12, sprintf('huberc revs about the Moon = %.2f', F.revs));
ok = chk(ok, isfinite(F.minAbsDQdt) && F.minAbsDQdt > 0, ...
         sprintf('huberc min|dQ/dt| at crossings = %.3f', F.minAbsDQdt));

% (2) the WALLED arms, same cell and gamma --------------------------------
Fe = regime_features(armOf(fullfile(resDir, 'minfuel_hg_c25_g200_eps.mat')), rec, 'eps');
ok = chk(ok, strcmp(Fe.outcome, 'wall'), sprintf('eps outcome = %s (want wall)', Fe.outcome));
ok = chk(ok, abs(Fe.rampWidth - 2*Fe.pFloor) < 1e-12, ...
         sprintf('eps ramp width = 2p = %.4g', Fe.rampWidth));
ok = chk(ok, isnan(Fe.deltaFloor), 'eps has no delta');
Fh = regime_features(armOf(fullfile(resDir, 'minfuel_hg_c25_g200_huber.mat')), rec, 'huber');
ok = chk(ok, strcmp(Fh.outcome, 'wall'), sprintf('huber outcome = %s (want wall)', Fh.outcome));
ok = chk(ok, Fh.rampWidth == 0, 'huber ramp width = 0 (the jump)');

% sharpness ordering at the wall/floor: huber 0 < huberc delta < eps 2p
ok = chk(ok, Fh.rampWidth < F.rampWidth && F.rampWidth < Fe.rampWidth, ...
         sprintf('common sharpness axis: huber %.4g < huberc %.4g < eps %.4g', ...
                 Fh.rampWidth, F.rampWidth, Fe.rampWidth));

% (3) a BISECTED floor: the walk's last accepted rung is an inserted point
%     just above the schedule's 1e-3, not 1e-3 itself. Real case: eps on
%     (1,2) gamma = 1.223 stopped at p = 0.001311 and IS the shipped
%     catalog entry -- outcome must be 'floor', not 'wall'.
B = A;  B.p(end) = 0.001311;  B.delta(end) = 0.0026;
Fb = regime_features(B, rec, 'eps');
ok = chk(ok, strcmp(Fb.outcome, 'floor'), ...
         sprintf('bisected floor p = %.4g -> outcome %s (want floor)', Fb.pFloor, Fb.outcome));
C = A;  C.p(end) = 0.0021;
Fc = regime_features(C, rec, 'eps');
ok = chk(ok, strcmp(Fc.outcome, 'wall'), ...
         sprintf('p = %.4g is NOT the floor -> outcome %s (want wall)', Fc.pFloor, Fc.outcome));

% (4) an arm that never accepted a rung -----------------------------------
E = A;  E.p = [];  E.mf = [];  E.delta = [];  E.Y = {};  E.coastFrac = [];  E.condJ = [];
Fn = regime_features(E, rec, 'huber');
ok = chk(ok, strcmp(Fn.outcome, 'none') && isnan(Fn.pFloor), ...
         sprintf('empty arm -> outcome %s, no crash', Fn.outcome));

if ok, fprintf('TEST_REGIME_FEATURES: ALL PASS\n');
else,  fprintf('TEST_REGIME_FEATURES: FAILURE (see lines above)\n');
end
end

function A = armOf(f)
% ARMOF  The single arm stored in a race output file.  INPUTS: f. OUTPUTS: A.
L = load(f);  fn = fieldnames(L.out.arms);  A = L.out.arms.(fn{1});
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
