function pool = capped_pool(nWorkers)
%% Purpose:
%
%   The parallel pool that run_capped's hard wall-clock fence needs, with a
%   PRIVATE JobStorageLocation so that concurrently running `matlab -batch`
%   sessions do not collide over the default one.
%
%   Workers inherit the client's path at pool creation, so this must be
%   called AFTER the campaign's startup has put pumpkyn and the library
%   folders on the path -- otherwise every fenced call fails on the worker
%   for a missing function rather than on its own merits.
%
%% Inputs:
%
%  nWorkers                 double (optional)       [4]
%
%% Outputs:
%
%  pool                     parallel.Pool           existing pool if one is
%                                                   already open, else a new
%                                                   one; [] if the Parallel
%                                                   Computing Toolbox is
%                                                   unavailable
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1 || isempty(nWorkers), nWorkers = 4; end
pool = [];
if isempty(ver('parallel')), warning('capped_pool:noPCT', 'Parallel Computing Toolbox unavailable: calls will run UNFENCED'); return, end
pool = current_pool();
if ~isempty(pool) && isvalid(pool), return, end
% A MISSING FENCE MUST DEGRADE, NOT CRASH. The licence can be held by
% another MATLAB session on this machine (a shared desktop session is
% enough), and a campaign that dies at pool creation loses the whole run
% rather than the fence. Consumers already branch on an empty pool.
try
    c = parcluster('Processes');
    tmp = fullfile(tempdir, sprintf('capped_pool_%d_%s', feature('getpid'), datestr(now, 'HHMMSSFFF')));
    if ~isfolder(tmp), mkdir(tmp); end
    c.JobStorageLocation = tmp;
    pool = c.parpool(min(nWorkers, c.NumWorkers));
catch ME
    pool = [];
    warning('capped_pool:noPool', 'could not open a pool (%s): calls will run UNFENCED', ME.message);
end
end
