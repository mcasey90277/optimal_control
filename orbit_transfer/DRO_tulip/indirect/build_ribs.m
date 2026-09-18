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
%   .nD [12] .nPts [nD-1] .wallSec [900] per point, .targets [] explicit
%   unwrapped departure offsets from the spine (rib_targets; beats nD/nPts),
%   .only [] grid columns
%   to walk (default: every certified one), .out [results/arrival_ribs.mat]
%   (published atomically: written to an exclusive temp name, then one rename),
%   .progress [] handle called after every solve (rib_from_crossing),
%   .checkpoint '' resumable checkpoint path, or a handle of the column j
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

% TEST SEAM: handles to this file's local functions (tests/test_fill_holes_physics)
if ischar(sheetMat) && strcmp(sheetMat, 'localfunctions'), R = localHandles(localfunctions);  return, end
if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'), fullfile(fileparts(here), '..', 'campaign_common'));
L = load(sheetMat);  S = L.S;
nD = d('nD', 12);  nPts = d('nPts', nD - 1);  dirn = d('direction', -1);
targets = d('targets', []);
if ~isempty(targets), nPts = numel(targets); end
out = d('out', fullfile(here, 'results', 'arrival_ribs.mat'));
pool = capped_pool();
% RECONSTRUCT FROM THE SHEET, not from defaults. The builder used to rebuild
% the default setup, so a sheet certified at a non-default operating point
% could acquire ribs generated at a different one -- valid trajectories, but
% not valid additions to that sheet. (Astra chain review 2026-09-10.)
% CLOSURES ONLY. A rib starts from a root the SHEET already certified, so it
% needs the endpoint closures and the propulsion constants and nothing else.
% It used to take the full setup, which re-polishes the shipped 70 mN anchor
% and so fails for any other engine, orbit pair or departure phase
% (FINDINGS 78). The walker reads one thing from `anc`: the spine's phase.
so = setupFromSheet(S, d('setupOpts', struct()));
[B, ~] = arclength_arrival('setup', so);
anc = struct('sD', B.problem.sD);
if isfield(S, 'problem')
    % check the WHOLE identity, not two of its fields: the assert existed to
    % catch a rib built at another operating point, and thrust plus phase
    % would have passed a tulip of the wrong petal count or branch.
    idf = {'sD', 'thrustN', 'ispS', 'm0kg', 'tauDRO', 'NpTulip', 'pmTulip'};
    for kf = 1:numel(idf)
        f = idf{kf};
        if ~isfield(S.problem, f), continue, end
        assert(abs(B.problem.(f) - S.problem.(f)) <= 1e-12*max(1, abs(S.problem.(f))), ...
               'build_ribs:identity', ...
               ['rib setup does not reproduce the sheet''s problem identity: ' ...
                '%s is %g here and %g in the sheet'], f, B.problem.(f), S.problem.(f));
    end
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
    ck = d('checkpoint', '');
    if isa(ck, 'function_handle'), ck = ck(j); end
    Rj = rib_from_crossing(c(k), B, anc, struct('nD', nD, 'direction', dirn, ...
        'nPts', nPts, 'targets', targets, 'wallSec', d('wallSec', 900), 'copts', struct('pool', pool), ...
        'progress', d('progress', []), 'checkpoint', ck, ...
        'problem', pickField(S, 'problem', struct())));
    R(end+1) = struct('j', j, 'sA', S.sA(j), 'pts', Rj.pts, 'stop', Rj.stop, ...
                      'nSolve', Rj.nSolve); %#ok<AGROW>
    fprintf('  -> %d certified points, %d solves, %.0f s, %s\n', ...
        numel(Rj.pts), Rj.nSolve, toc(t0), Rj.stop);
    if isfield(S, 'problem'), problem = S.problem; else, problem = struct(); end
    % after every rib, not at the end -- and PUBLISHED ATOMICALLY: the
    % campaign queue reads "this file exists" as "this unit is done", so a
    % save interrupted half-way must not leave a file behind
    tmp = sprintf('%s.%s.part', out, char(java.util.UUID.randomUUID()));   % exclusive name
    save(tmp, 'R', 'problem');
    publish_atomic(tmp, out);                       % one rename, or an error
end
fprintf('build_ribs: %d ribs, %d certified points -> %s\n', ...
    numel(R), sum(arrayfun(@(r) numel(r.pts), R)), out);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function v = pickField(s, f, d_)
% PICKFIELD  Field with default (present even if empty).  INPUTS: s; f; d_.
% OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = d_; end
end

function so = setupFromSheet(S, so)
% SETUPFROMSHEET  The setup request that reproduces the SHEET's problem:
% engine, orbits, branch and the spine's departure phase, closures only.
% INPUTS: S (sheet; .problem optional); so (caller's setup options).
% OUTPUTS: so struct for arclength_arrival('setup', so).
if isfield(S, 'problem')
    P = S.problem;
    so.thrustN = P.thrustN;  so.ispS = P.ispS;  so.m0kg = P.m0kg;
    so.tauDRO = P.tauDRO;    so.NpTulip = P.NpTulip;  so.sD = P.sD;
    % the BRANCH travels with the petal count. It was a hardcoded -1 in the
    % setup, so carrying Np without pm used to be harmless; it is not now.
    if isfield(P, 'pmTulip'), so.pmTulip = P.pmTulip; end
end
so.physicsOnly = true;
end

function H = localHandles(fh)
% LOCALHANDLES  This file's local functions as a struct of handles keyed by
% name -- the TEST SEAM: a test calls the real helper, not a copy of it.
% INPUTS: fh (cell of handles, from localfunctions).  OUTPUTS: H struct.
H = struct();
for k = 1:numel(fh), H.(func2str(fh{k})) = fh{k}; end
end
