function casadi_path()
%% Purpose:
%
%   Put CasADi on the path if casadi.Opti is not already there, from
%   ~/casadi-3.7.0 under the user's home. No absolute path is written down.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('casadi.Opti'))
    addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
end
if isempty(which('casadi.Opti'))
    error('casadi_path:missing', 'CasADi not found (expected ~/casadi-3.7.0)');
end
end
