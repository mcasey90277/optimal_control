function out = ms_conjugate_test(info, spec)
%% Purpose:
%
%   DELEGATE since 2026-09-17: the Jacobi (conjugate-point) test was promoted
%   to the cross-folder library, oclib/+oc/ms_conjugate_test, when the
%   cart-pole study scripts (optimal_control_examples) became its second
%   TOP-LEVEL consumer. This file keeps every costate_common caller working;
%   the instrument, its conventions and its tests live with the
%   implementation.
%
%% Inputs:
%
%  info                     struct                  see oc.ms_conjugate_test
%
%  spec                     struct (optional)       see oc.ms_conjugate_test
%
%% Outputs:
%
%  out                      struct                  see oc.ms_conjugate_test
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('oc.ms_conjugate_test'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'oclib'));
end
if nargin < 2, spec = struct(); end
out = oc.ms_conjugate_test(info, spec);
end
