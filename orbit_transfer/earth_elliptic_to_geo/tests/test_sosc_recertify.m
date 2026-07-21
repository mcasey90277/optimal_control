% test_sosc_recertify  recertify_table3 integration: one real rung (10 N),
% no heavy loop. Confirms the sidecar is written, the campaign .mat cache is
% left byte-unchanged (never clobber production caches, DESIGN_sosc.md sec 8),
% and the returned row is well-formed against the full 5-way verdict set
% {PASS, WEAK_MIN, FAIL, INCONCLUSIVE, ERROR} (sec 11.5).
%
% This runs one real IPOPT warm re-solve (~1-2 min).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

campaign = fullfile(module_root(),'results','MEE_M2_10N.mat');
info0  = dir(campaign); before = info0.bytes;

sidecarPath = fullfile(module_root(),'results','sosc','sosc_MEE_M2_10N.mat');
if isfile(sidecarPath), delete(sidecarPath); end   % force a fresh write this run

T = recertify_table3(10);

assert(isscalar(T), 'one row for 10 N');
assert(T(1).thrustN==10, 'thrustN echoed');
assert(strcmp(T(1).tag,'MEE_M2_10N'), 'tag matches headline row');
assert(ismember(T(1).verdict, {'PASS','WEAK_MIN','FAIL','INCONCLUSIVE','ERROR'}), ...
    'verdict in the full sec-11.5 taxonomy');
assert(isfile(sidecarPath), 'sidecar written');

sc = load(sidecarPath);
assert(isfield(sc,'sosc') && strcmp(sc.sosc.verdict, T(1).verdict), ...
    'sidecar verdict matches returned row');

info1 = dir(campaign);
assert(info1.bytes==before, 'campaign .mat must be untouched (byte-identical)');

fprintf('test_sosc_recertify PASSED (10 N verdict=%s, nFlat=%d, method=%s, robust=%d)\n', ...
    T(1).verdict, T(1).nFlat, T(1).method, T(1).robust);
