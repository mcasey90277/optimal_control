function ok = test_audit_fail_closed()
% TEST_AUDIT_FAIL_CLOSED  The catalog audit must fail CLOSED (FINDINGS 78):
% a re-polish that times out, a hypothesis call that fails, a recomputed
% gate that is not positive, or a polished root that moved, is a BAD row --
% never a clean one. One heavy call at a time is knocked out through
% opts.cappedWrap (which is handed the audit's real fence), so the
% test runs in seconds on one entry of the library of record.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
catMat = fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
assert(isfile(catMat), 'the library of record is missing: %s', catMat);
base = struct('idx', 1, 'pool', []);

% ---- a polish that times out is BAD, not skipped --------------------------
A = audit_phase_catalog(catMat, setfield(base, 'cappedWrap', @(realCap) @(varargin) failing(realCap, varargin, 'ms_tfmin'))); %#ok<SFLD>
ok = chk(ok, A.nBad == 1 && contains(A.rows(1).problem, 'polish'), sprintf('polish timeout -> BAD (%s)', A.rows(1).problem));

% ---- a hypothesis call that fails is BAD ----------------------------------
A = audit_phase_catalog(catMat, setfield(base, 'cappedWrap', @(realCap) @(varargin) failing(realCap, varargin, 'mintime_hypothesis_gates'))); %#ok<SFLD>
ok = chk(ok, A.nBad == 1 && contains(A.rows(1).problem, 'hypothesis'), sprintf('gates failure -> BAD (%s)', A.rows(1).problem));

% ---- a recomputed gate that is not positive is BAD ------------------------
A = audit_phase_catalog(catMat, setfield(base, 'cappedWrap', @(realCap) @(varargin) tampered(realCap, varargin, 'minLamV', -1))); %#ok<SFLD>
ok = chk(ok, A.nBad == 1 && contains(A.rows(1).problem, 'lam_v'), sprintf('min|lam_v| <= 0 -> BAD (%s)', A.rows(1).problem));

% ---- a polished root that moved away from the stored one is BAD -----------
A = audit_phase_catalog(catMat, setfield(base, 'cappedWrap', @(realCap) @(varargin) movedRoot(realCap, varargin))); %#ok<SFLD>
ok = chk(ok, A.nBad == 1 && contains(A.rows(1).problem, 'moved'), sprintf('polished root moved -> BAD (%s)', A.rows(1).problem));

if ok, fprintf('test_audit_fail_closed: ALL PASS\n'); else, fprintf('test_audit_fail_closed: FAIL\n'); end
end

% ---- the injected capped: real for everything but the named call ----------
function varargout = failing(realCap, args, name)
% FAILING  The named capped call fails (ok = false, empty outputs); every
% other call is the real fence.  INPUTS: realCap; args (its arguments); name.
% OUTPUTS: as capped.
fn = args{3};  nout = args{4};
if strcmp(func2str(fn), name) || strcmp(func2str(fn), ['@' name])
    varargout = [{false}, repmat({[]}, 1, nout)];
else
    [varargout{1:nout+1}] = realCap(args{:});
end
end

function varargout = tampered(realCap, args, field, value)
% TAMPERED  The hypothesis-gates call returns its real result with one field
% overwritten.  INPUTS: realCap; args; field; value.  OUTPUTS: as capped.
fn = args{3};  nout = args{4};
[varargout{1:nout+1}] = realCap(args{:});
if contains(func2str(fn), 'mintime_hypothesis_gates') && varargout{1}
    varargout{2}.(field) = value;
end
end

function varargout = movedRoot(realCap, args)
% MOVEDROOT  The polish returns a root 1e-3 away from the stored one.
% INPUTS: realCap; args.  OUTPUTS: as capped.
fn = args{3};  nout = args{4};
[varargout{1:nout+1}] = realCap(args{:});
if contains(func2str(fn), 'ms_tfmin') && varargout{1}
    varargout{2} = varargout{2} + 1e-3;                       % z, shifted
end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
