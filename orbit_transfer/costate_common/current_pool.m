function pool = current_pool()
%% Purpose:
%
%   The open parallel pool, or [] when there is none -- and, unlike a bare
%   gcp('nocreate'), it does not THROW when the Parallel Computing Toolbox
%   is missing or its licence is held by another MATLAB session on this
%   machine. A shared desktop session holding the PCT seat is enough to make
%   every `matlab -batch` job on the same machine poolless, so the fence
%   must degrade to an unfenced call rather than crash the campaign.
%
%   Callers that already hold a pool must pass it. Callers that want a pool
%   CREATED should use capped_pool; this function never creates one.
%
%  ASSUMPTIONS / NOTES:
%
% • [] means UNFENCED: run_capped hands its first argument straight to
%   parfeval, which refuses an empty pool, so every consumer must branch on
%   emptiness before calling it.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  pool                     parallel.Pool or []     the open pool, or [] if
%                                                   there is none available
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

pool = [];
if isempty(ver('parallel')), return, end
try
    pool = gcp('nocreate');
catch
    pool = [];
end
end
