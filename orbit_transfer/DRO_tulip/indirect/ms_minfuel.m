function [z, info] = ms_minfuel(rv0, rvf, tf, seed, Tmax, c, muStar, smooth, opts)
% MS_MINFUEL  Multiple-shooting solve of the CR3BP fixed-tf SMOOTHED
% energy->fuel PMP problem -- the homotopy sibling of ms_minenergy.
%
% Cost J = Int L(s) dt with L selected by the smoothing spec (see
% cr3bp_minfuel_pmp: 'eps' Bertrand-Epenoy or 'huber' PLQ), final time
% FIXED at tf, final mass free. Terminal conditions identical to
% ms_minenergy: r(tf) = rf, v(tf) = vf, lambda_m(tf) = 0. H is a first
% integral and is reported (info.H/.Hdrift).
%
% Thin PROBLEM DEFINITION bound to costate_common/ms_bvp (opts.fixedTf):
% only the closures differ from ms_minenergy -- the field and propagator
% carry the smoothing spec. Optionally runs ss_bvp_accept on the same
% closures (opts.accept).
%
% INPUTS:
%   rv0, rvf - endpoint rotating-frame pos/vel, ND [6x1]
%   tf       - FIXED time of flight, ND [scalar]
%   seed     - struct .tGrid [1 x K+1], .Y [14 x K+1] (col-1 state part
%              overwritten with [rv0; 1])
%   Tmax, c, muStar - as cr3bp_minfuel_pmp
%   smooth   - struct('family','eps'|'huber', 'p', value)
%   opts     - (optional) ms_bvp options + .accept [false], .tolDz [1e-6],
%              .tolRss [1e-6] (as ms_minenergy)
%
% OUTPUTS:
%   z    - lambda0 [7x1]
%   info - as ms_minenergy (.Y junction starts, .H, .Hdrift, .s throttle
%          at junctions, .coastFrac = fraction of junctions with s < 1e-3,
%          .accept when requested)
%
% REFERENCES:
%   [1] costate_common/cr3bp_minfuel_pmp (the field + family notes).
%   [2] ms_minenergy.m (the eps = 1 anchor this generalizes).

if nargin < 9, opts = struct(); end
if isempty(which('ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts( ...
        mfilename('fullpath')))), 'costate_common'));
end
rv0 = rv0(:);  rvf = rvf(:);

prob = struct('ny', 14, 'freeIdx0', 8:14, ...
    'prop',     @(dt, y0, needSTM) cr3bp_minfuel_prop(dt, y0, needSTM, Tmax, c, muStar, smooth), ...
    'rhs',      @(y) cr3bp_minfuel_pmp(y, Tmax, c, muStar, smooth), ...
    'terminal', @(y, needJ) terminalMinFuel(y, rvf));

doAccept = isfield(opts, 'accept') && opts.accept;
tolDz = 1e-6;   tolRss = 1e-6;
if isfield(opts, 'tolDz'),  tolDz  = opts.tolDz;  end
if isfield(opts, 'tolRss'), tolRss = opts.tolRss; end
mo = opts;
for f = {'accept', 'tolDz', 'tolRss'}
    if isfield(mo, f{1}), mo = rmfield(mo, f{1}); end
end
mo.fixedTf = true;

conjTest = isfield(opts, 'conjTest') && opts.conjTest;
if conjTest, mo.keepSTMs = true; end
if isfield(mo, 'conjTest'), mo = rmfield(mo, 'conjTest'); end

seed.tf = tf;
seed.Y(1:7,1) = [rv0; 1];
[p, info] = ms_bvp(prob, seed, mo);
z = p(1:7);

if conjTest
    % Fixed-tf Jacobi spec as in ms_minenergy (validated on the analytic
    % LQ case): no flow column, no quotient, rows [1:6 14], cols 8:14.
    info.conj = ms_conjugate_test(info, struct( ...
        'stateRows', [1:6 14], 'costateCols', 8:14, ...
        'quotientDir', [], 'freeTime', false));
end

K = size(info.Y, 2);
info.H = zeros(1, K);  info.s = zeros(1, K);
for k = 1:K
    [~, ~, ax] = cr3bp_minfuel_pmp(info.Y(:,k), Tmax, c, muStar, smooth);
    info.H(k) = ax.H;  info.s(k) = ax.s;
end
info.Hdrift = max(abs(info.H - info.H(1)));
info.coastFrac = nnz(info.s < 1e-3) / K;

if doAccept
    y1 = info.Y(:,1);
    ao = struct('fixedTf', true, 'tolDz', tolDz, 'tolR', tolRss);
    [~, info.accept] = ss_bvp_accept(prob, y1, tf, ao);
end
end

% ------------------------------------------------------------------------
function [g, dgdy] = terminalMinFuel(y, rvf)
% TERMINALMINFUEL  Fixed-tf terminal conditions: r, v matched;
% lambda_m(tf) = 0 (final mass free).
% INPUTS: y [14x1]; rvf [6x1].  OUTPUTS: g [7x1]; dgdy [7x14].
g = [y(1:6) - rvf; y(14)];
dgdy = [eye(6), zeros(6, 8);
        zeros(1, 13), 1];
end
