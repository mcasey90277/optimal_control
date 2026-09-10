function R = build_ribs(sheetMat, opts)
%% Purpose:
%
%   Walk a DEPARTURE-phase rib off every certified point of an arrival
%   sheet, and save them. This is the second axis of the (departure x
%   arrival) torus: the arrival direction is walked by continuation
%   (arclength_arrival, where the folds are), the departure direction by a
%   bisecting walker (rib_from_crossing, where the steps are cheap), and
%   every point on either axis goes through the same gate stack.
%
%   Ribs are walked one per certified arrival phase, independently, so a
%   rib that stalls costs only its own row.
%
%% Inputs:
%
%  sheetMat                 char                    an arrival sheet .mat
%                                                   (build_arrival_sheet)
%  opts                     struct (optional)
%   .direction [-1] the departure sense that works out of the anchor,
%   .nD [12] .nPts [nD-1] .wallSec [900] per point, .only [] grid columns
%   to walk (default: every certified one), .out [results/arrival_ribs.mat]
%
%% Outputs:
%
%  R                        struct array            one per rib: .j .sA
%                                                   .pts .stop .nSolve
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
L = load(sheetMat);  S = L.S;
nD = d('nD', 12);  nPts = d('nPts', nD - 1);  dirn = d('direction', -1);
out = d('out', fullfile(here, 'results', 'arrival_ribs.mat'));
pool = capped_pool();
% RECONSTRUCT FROM THE SHEET, not from defaults. The builder used to rebuild
% the default setup, so a sheet certified at a non-default operating point
% could acquire ribs generated at a different one -- valid trajectories, but
% not valid additions to that sheet. (Astra chain review 2026-09-10.)
so = d('setupOpts', struct());
if isfield(S, 'problem')
    P = S.problem;
    so.thrustN = P.thrustN;  so.ispS = P.ispS;  so.m0kg = P.m0kg;
    so.tauDRO = P.tauDRO;    so.NpTulip = P.NpTulip;  so.sD = P.sD;
end
[B, anc] = arclength_arrival('setup', so);
if isfield(S, 'problem')
    assert(abs(B.problem.sD - S.problem.sD) < 1e-12 && ...
           abs(B.problem.thrustN - S.problem.thrustN) < 1e-12, ...
           'rib setup does not reproduce the sheet''s problem identity');
end

cols = d('only', find(isfinite(S.TF)));
R = struct('j', {}, 'sA', {}, 'pts', {}, 'stop', {}, 'nSolve', {});
for j = cols(:)'
    c = S.cand{j};
    k = find([c.ok] & abs([c.tfDays] - S.TF(j)) < 1e-9, 1);
    if isempty(k), continue, end
    fprintf('rib at grid %d (sA %.4f, t_f %.4f d), %d points, direction %+d\n', ...
        j, S.sA(j), S.TF(j), nPts, dirn);
    t0 = tic;
    Rj = rib_from_crossing(c(k), B, anc, struct('nD', nD, 'direction', dirn, ...
        'nPts', nPts, 'wallSec', d('wallSec', 900), 'copts', struct('pool', pool)));
    R(end+1) = struct('j', j, 'sA', S.sA(j), 'pts', Rj.pts, 'stop', Rj.stop, ...
                      'nSolve', Rj.nSolve); %#ok<AGROW>
    fprintf('  -> %d certified points, %d solves, %.0f s, %s\n', ...
        numel(Rj.pts), Rj.nSolve, toc(t0), Rj.stop);
    if isfield(S, 'problem'), problem = S.problem; else, problem = struct(); end %#ok<NASGU>
    save(out, 'R', 'problem');            % after every rib, not at the end
end
fprintf('build_ribs: %d ribs, %d certified points -> %s\n', ...
    numel(R), sum(arrayfun(@(r) numel(r.pts), R)), out);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
