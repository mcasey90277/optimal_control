function [ok, dz, z8out, info] = tfmin_accepts(cat, entry, opts)
%% Purpose:
%
%   The "pumpkyn tfMin accepts this entry unchanged" check: hand the
%   entry's z8 to pumpkyn.cr3bp.tfMin as its initial guess, at the entry's
%   own endpoints, and require the returned vector to stay within tolDz of
%   it. The call runs under a hard wall-clock cap (run_capped) when a pool
%   exists, directly otherwise.
%
%  ASSUMPTIONS / NOTES:
%
% • dz is the certifier's gate quantity (certify_root's C.dz, the tfMin
%   witness): the 2-norm of (z8out - z8), written sqrt(sum(.^2)).
% • Ruling R9: a status-1 entry ("flies to the target, necessary conditions
%   not established") is not expected to be accepted, so tfMin is NOT run
%   on it (it would cost the full cap): ok = false, dz = NaN, z8out = [],
%   info.expected = true, info.ran = false. opts.force = true runs it
%   anyway (expected stays true). Never thrown.
% • A cap timeout, a thrown tfMin and an unusable return are all refusals
%   (ok = false) with dz = NaN, named in info; nothing is thrown for them.
% • Entry forms as fly_catalog_entry: a transfers_at row, struct('iD',,
%   'iA',) (the primary of that cell of sheet 1), or struct('z8',,'sD',,
%   'sA',). Phases are FRACTIONS of a period.
% • run_capped reports a worker error and a timeout alike as ok = false; the
%   two are told apart by the elapsed time against the cap.
% • The default B is catalog_endpoints(cat): sheet 1 at rungs_N(1) (this
%   library has one sheet and one rung); pass opts.B for anything else.
%
%% Inputs:
%
%  cat                      struct                  library catalog
%  entry                    struct                  entry in any of the
%                                                   three forms above
%  opts                     struct (optional)       .tolDz [1e-6]; .pool
%                                                   [current_pool()];
%                                                   .capSec [300]; .B []
%                                                   (a catalog_endpoints
%                                                   result); .tfminFcn
%                                                   [@pumpkyn.cr3bp.tfMin]
%                                                   test seam; .force
%                                                   [false] run tfMin even
%                                                   on a status-1 entry
%
%% Outputs:
%
%  ok                       logical                 isfinite(dz) && dz <= tolDz
%  dz                       double                  scalar, 2-norm(z8out-z8); NaN
%                                                   when tfMin gave nothing
%  z8out                    double                  [8x1] tfMin's answer ([]
%                                                   when none)
%  info                     struct                  .expected .ran .timedOut
%                                                   .threw .fenced .message
%                                                   .status .wallSec .tolDz
%                                                   .capSec
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  M. Casey  ASSUMPTIONS: the default B is sheet 1 at rungs_N(1)
%                                                             (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(opts), opts = struct(); end
tolDz   = 1e-6;  if isfield(opts, 'tolDz')   && ~isempty(opts.tolDz),   tolDz   = opts.tolDz;   end
capSec  = 300;   if isfield(opts, 'capSec')  && ~isempty(opts.capSec),  capSec  = opts.capSec;  end
tfminFh = @pumpkyn.cr3bp.tfMin;
if isfield(opts, 'tfminFcn') && ~isempty(opts.tfminFcn), tfminFh = opts.tfminFcn; end
if isfield(opts, 'pool'), pool = opts.pool; else, pool = current_pool(); end

%% Resolve the entry form:
if isfield(entry, 'iD') && isfield(entry, 'iA') && ~isfield(entry, 'z8')
    sh = cat.sheets(1);
    T  = transfers_at(cat, sh.sD_frac(entry.iD), sh.sA_frac(entry.iA));
    ip = find(strcmp({T.kind}, 'primary'), 1);
    if isempty(ip), error('tfmin_accepts:noPrimary', 'cell (%d,%d) has no primary', entry.iD, entry.iA); end
    entry = T(ip);
end
if ~isfield(entry, 'z8') || numel(entry.z8) ~= 8
    error('tfmin_accepts:entry', 'the entry needs a z8 of 8 elements');
end
z8 = entry.z8(:);

status = NaN;
if isfield(entry, 'status') && ~isempty(entry.status) && isfinite(entry.status), status = entry.status; end

if isfield(opts, 'B') && ~isempty(opts.B), B = opts.B; else, B = catalog_endpoints(cat); end
rv0 = B.stateD(entry.sD);  rv0 = rv0(1:6);
rvf = B.stateA(entry.sA);  rvf = rvf(1:6);

info = struct('expected', status == 1, 'ran', false, 'timedOut', false, 'threw', '', 'fenced', ~isempty(pool), ...
              'message', '', 'status', status, 'wallSec', NaN, 'tolDz', tolDz, 'capSec', capSec);
dz = NaN;  z8out = [];  ok = false;
force = false;  if isfield(opts, 'force') && ~isempty(opts.force), force = logical(opts.force); end
if info.expected && ~force
    info.message = ['the certification witness failed on this entry at build time; tfMin is not expected to accept it ' ...
                    'and was not run (opts.force = true runs it anyway)'];
    return
end
info.ran = true;

%% Run tfMin, fenced when a pool exists:
t0 = tic;
try
    if isempty(pool)
        za = tfminFh(rv0', rvf', z8, B.Tnd, B.cnd, B.mu);
        done = true;
    else
        [done, za] = run_capped(pool, tfminFh, 1, capSec, rv0', rvf', z8, B.Tnd, B.cnd, B.mu);
    end
catch ME
    info.wallSec = toc(t0);
    info.threw = ME.message;
    info.message = sprintf('tfMin threw: %s', ME.message);
    return
end
info.wallSec = toc(t0);

if ~done
    if info.wallSec >= capSec
        info.timedOut = true;
        info.message = sprintf('tfMin exceeded its %g s cap', capSec);
    else
        info.threw = 'tfMin failed on the worker';
        info.message = 'tfMin threw an error on the worker before its cap';
    end
    return
end
if ~(isnumeric(za) && numel(za) == 8 && all(isfinite(za(:))))
    info.message = 'tfMin returned an unusable vector';
    return
end

z8out = za(:);
d  = z8out - z8;
dz = sqrt(sum(d.^2));
ok = isfinite(dz) && dz <= tolDz;

%% Message:
if ok
    info.message = sprintf('tfMin accepts the entry unchanged (|dz| = %.2e <= %g).', dz, tolDz);
elseif info.expected
    info.message = sprintf('status-1 entry: the certification witness failed on it, so tfMin is not expected to accept it (|dz| = %.2e).', dz);
else
    info.message = sprintf('tfMin moved the entry (|dz| = %.2e > %g).', dz, tolDz);
end
end
