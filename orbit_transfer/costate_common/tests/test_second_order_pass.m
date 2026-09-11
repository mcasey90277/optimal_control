function ok = test_second_order_pass()
%% Purpose:
%
%   Tests second_order_pass on ONE catalog entry with a fresh sidecar: the
%   three instruments are recorded, the lift margin is measured with the
%   TIGHT setting pair (1e-12 / 1e-9 -- the loose 1e-7 second setting
%   inflated the error estimate and read 4-9x on seven long-arc entries
%   that measure 40x+ with the tight pair), and the spectrum's candidates
%   are recorded by class rather than as a count.
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
catMat = fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', 'costate_catalog_dro_tulip_70mN.mat');
if ~isfile(catMat), fprintf('  SKIP  no 70 mN catalog on disk\n'); return, end
side = fullfile(tempdir, sprintf('sop_test_%d.mat', round(1e6*rand)));
S = second_order_pass(catMat, struct('sideMat', side, 'maxEntries', 1, 'logFile', ''));
ok = chk(ok, S.nDone == 1 && ~S.done, sprintf('one entry measured, census %d done / %d to do', S.nDone, S.nTodo));
r = S.rows(find([S.rows.done], 1));
ok = chk(ok, isfield(r, 'relTolPair') && isequal(r.relTolPair, [1e-12 1e-9]), 'lift margin measured with the TIGHT pair');
ok = chk(ok, r.liftMargin >= 10 && r.liftCertified, sprintf('lift margin %.0fx certified', r.liftMargin));
ok = chk(ok, isfield(r, 'nInteriorCand') && isfield(r, 'nNearMiss') && isfield(r, 'nZero'), ...
         'spectrum candidates recorded by class');
ok = chk(ok, isfield(r, 'candidates') && isstruct(r.candidates), 'and the located candidates themselves are kept');
if isfile(side), delete(side); end
if ok, fprintf('TEST_SECOND_ORDER_PASS: ALL PASS\n'); else, fprintf('TEST_SECOND_ORDER_PASS: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
