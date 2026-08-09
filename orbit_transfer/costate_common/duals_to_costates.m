function [lam, tStations, diag_] = duals_to_costates(spec)
%% Purpose:
%
%   DELEGATE. The covector mapping -- defect-constraint KKT multipliers to
%   continuous-costate samples, with the scheme-specific station
%   association, sign vote, and lambda_t check -- was promoted to the
%   cross-folder optimal-control library on 2026-08-09 and lives at
%
%       oclib/+oc/duals_to_costates.m        (call as oc.duals_to_costates)
%
%   because a SECOND top-level consumer exists: booster_landing's G5
%   primer gate independently re-invented (and independently re-fixed) the
%   Hermite-Simpson midpoint-vs-node association bug this function owns.
%   One home per subtle rule. This delegate keeps every costate_common
%   caller working unchanged.
%
%% Inputs:
%
%  spec                     struct                  See oc.duals_to_costates
%
%% Outputs:
%
%  lam, tStations, diag_                            See oc.duals_to_costates
%
%% Revision History:
%  M. Casey                                                   (c) 08/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('oc.duals_to_costates'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'oclib'));
end
if nargin == 0
    oc.duals_to_costates;                      % forward the self-demo
    return
end
[lam, tStations, diag_] = oc.duals_to_costates(spec);
end
