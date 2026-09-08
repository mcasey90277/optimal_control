function V = regime_verdicts(F, cfg)
%% Purpose:
%
%   Classify each (cell, gamma) case of the REGIME MAP from its family arms:
%   which families SOLVED it, which failed, and whether the failures sit at
%   a different switch count than the solution (the signature the regime
%   hypotheses turn on).
%
%   Two questions, two answers, and only the second is about the family:
%
%     procedural  did the p-walk reach the schedule's last rung?
%     PHYSICAL    did the walk find the optimum? -- measured as agreement of
%                 the final mass with the best arm of the same case.
%
%   They differ. Measured 2026-09-07: on (2,5) gamma 1.1 eps stops at
%   p = 0.0016 while huber and huberc reach 0.001, yet all three agree on
%   m_f to 2.4e-6 -- eps's remaining rungs were unnecessary, and calling it
%   "walled" would invent a family difference. On (1,2) gamma 1.223 huber
%   and huberc stop with m_f low by 9.4e-3 (1.4 kg of a 150 kg spacecraft)
%   -- that is a real failure. SOLVED therefore means dm_f <= tolSolve
%   against the best arm, whatever p the walk stopped at.
%
%% Inputs:
%
%  F                        struct array            regime_features output
%                                                   (.family .cellIdx .gamma
%                                                   .mf .pFloor .nCross ...)
%
%  cfg                      struct (optional)       .tolSolve [1e-3] mass
%                                                   agreement that counts as
%                                                   solving the case
%                                                   .tolTight [1e-5] the
%                                                   "indistinguishable" band
%
%% Outputs:
%
%  V                        struct array            one per case: .cellIdx
%                                                   .gamma .families .dMf
%                                                   .solved .failed .verdict
%                                                   .nSwitchSolved
%                                                   .nSwitchFailed
%                                                   .structureChange .mfBest
%                                                   .allTight
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, cfg = struct(); end
tolSolve = 1e-3;  if isfield(cfg, 'tolSolve'), tolSolve = cfg.tolSolve; end
tolTight = 1e-5;  if isfield(cfg, 'tolTight'), tolTight = cfg.tolTight; end

key = arrayfun(@(f) sprintf('%d%d_%.4f', f.cellIdx(1), f.cellIdx(2), f.gamma), ...
               F, 'UniformOutput', false);
[u, ~, ic] = unique(key);
V = struct([]);
for k = 1:numel(u)
    g = F(ic == k);
    mfBest = max([g.mf]);
    dMf = mfBest - [g.mf];
    isSolved = dMf <= tolSolve;
    v = struct('cellIdx', g(1).cellIdx, 'gamma', g(1).gamma, ...
        'families', {{g.family}}, 'dMf', dMf, 'mfBest', mfBest, ...
        'solved', {{g(isSolved).family}}, 'failed', {{g(~isSolved).family}}, ...
        'nSwitchSolved', NaN, 'nSwitchFailed', [], 'structureChange', false, ...
        'allTight', all(dMf <= tolTight), 'nArms', numel(g), 'verdict', '');
    % the switch structure of the BEST arm, and of the ones that failed
    [~, ib] = min(dMf);
    v.nSwitchSolved = g(ib).nCross;
    v.nSwitchFailed = [g(~isSolved).nCross];
    if ~isempty(v.nSwitchFailed) && ~isnan(v.nSwitchSolved)
        v.structureChange = any(v.nSwitchFailed ~= v.nSwitchSolved);
    end
    if numel(g) < 3
        v.verdict = sprintf('INCOMPLETE (%d/3)', numel(g));
    elseif numel(v.solved) == 3
        v.verdict = 'all three';
    elseif numel(v.solved) == 1
        v.verdict = sprintf('%s ONLY', v.solved{1});
    elseif numel(v.solved) == 2
        v.verdict = sprintf('%s + %s', v.solved{1}, v.solved{2});
    else
        v.verdict = 'NONE solved';
    end
    if isempty(V), V = v; else, V(end+1) = v; end %#ok<AGROW>
end
end
