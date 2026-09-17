function [p, info] = ms_bvp(prob, seed, opts)
%% Purpose:
%
%   DELEGATE since 2026-09-16. The generic multiple-shooting two-point BVP
%   engine was promoted to the cross-folder optimal-control library and
%   lives at
%
%       oclib/+oc/ms_bvp.m                   (call as oc.ms_bvp)
%
%   because a SECOND top-level consumer exists: the cart-pole PMP-BVP demo
%   (collocation_examples/ex3_cart_pole_pmp) solves its Pontryagin boundary
%   value problem through the same engine, with no orbit, no CR3BP quantity
%   and no campaign anywhere in it -- which is the property the engine's
%   problem-agnostic contract claims. One home per contract. This delegate
%   keeps every costate_common caller working unchanged; the contract, the
%   maths and the tests live with the implementation.
%
%% Inputs:
%
%  prob, seed, opts                                 See oc.ms_bvp
%
%% Outputs:
%
%  p, info                                          See oc.ms_bvp
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('oc.ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'oclib'));
end
if nargin == 0
    oc.ms_bvp;                                 % forward the self-demo
    return
end
if nargin < 3, opts = struct(); end
[p, info] = oc.ms_bvp(prob, seed, opts);
end
