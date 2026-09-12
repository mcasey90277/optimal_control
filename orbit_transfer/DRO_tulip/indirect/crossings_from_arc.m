function X = crossings_from_arc(A, levels, B, opts)
%% Purpose:
%
%   RE-SCAN a saved continuation arc for crossings of levels it was not
%   asked about when it was walked.
%
%   Why this exists. `arclength_arrival` records crossings of the grid
%   levels it is GIVEN, during the walk. Refining the arrival grid would
%   therefore appear to need the arcs walked again -- hours each. It does
%   not: the arc stores its full solution vector at every step (`A.p`, a
%   cell of [n x 1] over 1839-4001 steps), so any level inside the arc's
%   q-range can be bracketed after the fact and a starting guess produced
%   by interpolating between two adjacent steps. Adjacent steps differ by
%   dq ~ 2e-4, far inside the Newton basin, so `certify_crossing` polishes
%   the guess onto the level exactly as it would have during the walk.
%
%   A CROSSING IS A CONVERGED ROOT AT ITS LEVEL, not a bracket. The walk
%   interpolates and then corrects with a fixed-phase Newton, and so does
%   this -- through the same library corrector (newton_fixed_q), with the
%   same locality guard that rejects a correction which has wandered out of
%   the bracket it came from onto another branch. Without that step the
%   records carry converged = false and `sheet_from_arcs` refuses every one
%   of them, which is exactly what a first version of this did: 120 of 120
%   crossings reported "not converged (|R| = NaN)" -- its own flags read
%   back, not a solver failure.
%
%   The output is in the same shape `sheet_from_arcs` already consumes, so
%   a refined sheet is the existing pipeline with a different level list.
%
%  ASSUMPTIONS / NOTES:
%
% • FOLDS ARE KEPT. The arc is not monotone in q (it turns at folds), so a
%   level can be crossed several times and each crossing is a distinct
%   candidate -- that is the whole reason a grid point can hold more than
%   one local minimum. Every bracket is returned.
% • Levels are matched on the UNWRAPPED q the arc walks, and also at
%   q +/- 1, +/- 2 ... so a level is found however many periods the arc has
%   travelled. `sheet_from_arcs` folds them back with mod 1.
% • `.converged` is false and `.normR` NaN by construction: these are
%   SEEDS at the level, not roots. The certifier makes them roots.
%
%% Inputs:
%
%  A                        struct                  a saved arclength_arrival
%                                                   arc (.q [1 x N], .p
%                                                   {1 x N})
%  levels                   [1 x L]                 arrival phases wanted
%                                                   (fractions; matched on
%                                                   the unwrapped q)
%  B                        struct                  arclength_arrival setup:
%                                                   needs .res (resFactory),
%                                                   .Dx. Pass [] to skip the
%                                                   correction and return
%                                                   brackets only (for
%                                                   inspection, not for a
%                                                   sheet)
%  opts                     struct (optional)
%   .nWrap [3] how many periods either side to match a level at
%   .dedupe [1e-9] two brackets whose interpolated q differ by less than
%   this are the same crossing, .newtonTol [1e-9] .newtonMax [12] the
%   correction, .maxCorrFrac [2] reject a correction longer than this many
%   bracket spans (it has left the branch)
%
%% Outputs:
%
%  X                        struct array            one per crossing, in the
%                                                   order the arc meets them:
%                                                   .level (as given, folded)
%                                                   .q (unwrapped, = the
%                                                   matched level) .p
%                                                   [n x 1] interpolated
%                                                   .converged (of the
%                                                   correction) .normR
%                                                   .afterIndex (the step
%                                                   before the bracket)
%                                                   .frac (position in the
%                                                   bracket, 0..1)
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, B = []; end
if nargin < 4, opts = struct(); end
nWrap = fieldd(opts, 'nWrap', 3);  dedupe = fieldd(opts, 'dedupe', 1e-9);
nTol = fieldd(opts, 'newtonTol', 1e-9);  nMax = fieldd(opts, 'newtonMax', 12);
maxCorrFrac = fieldd(opts, 'maxCorrFrac', 2);
correct = ~isempty(B);
if correct
    assert(isfield(B, 'res') && isfield(B, 'Dx'), 'crossings_from_arc:setup', ...
           'B must carry .res (resFactory) and .Dx; pass [] to return brackets only');
end
assert(isfield(A, 'q') && isfield(A, 'p') && iscell(A.p) && numel(A.p) == numel(A.q), ...
       'crossings_from_arc:arc', 'A needs .q [1 x N] and .p {1 x N} from the same walk');
q = A.q(:).';  N = numel(q);
levels = levels(:).';

X = struct('level', {}, 'q', {}, 'p', {}, 'converged', {}, 'normR', {}, ...
           'afterIndex', {}, 'frac', {});
hit = zeros(0, 2);                        % [afterIndex, qCross] for de-duplication

for L = levels
    % the same phase reachable at any number of whole periods travelled
    for w = -nWrap:nWrap
        Lw = L + w;
        if Lw < min(q) - 1e-12 || Lw > max(q) + 1e-12, continue, end
        s = q - Lw;
        for i = 1:N-1
            if s(i) == 0 || (s(i) < 0) ~= (s(i+1) < 0)     % bracket or exact hit
                if s(i) == s(i+1), f = 0; else, f = s(i)/(s(i) - s(i+1)); end
                f = min(max(f, 0), 1);
                qc = q(i) + f*(q(i+1) - q(i));
                if ~isempty(hit) && any(hit(:,1) == i & abs(hit(:,2) - qc) < dedupe), continue, end
                hit(end+1, :) = [i, qc]; %#ok<AGROW>
                % linear interpolation between adjacent continuation steps
                % (~2e-4 apart in q), then the SAME fixed-phase correction
                % the walk applies, with the same locality guard
                pPred = (1 - f)*A.p{i} + f*A.p{i+1};
                pC = pPred;  cv = false;  nr = NaN;
                if correct
                    [pC, cv, nr] = newton_fixed_q(B.res, Lw, pPred, B.Dx, nTol, nMax);
                    span = max(norm((A.p{i+1} - A.p{i})./B.Dx), realmin);
                    if cv && norm((pC - pPred)./B.Dx) > maxCorrFrac*span
                        cv = false;                 % left its bracket: another branch
                    end
                end
                X(end+1) = struct('level', mod(L, 1), 'q', qc, 'p', pC, ...
                                  'converged', cv, 'normR', nr, ...
                                  'afterIndex', i, 'frac', f); %#ok<AGROW>
            end
        end
    end
end

% in the order the arc meets them
if ~isempty(X)
    [~, ord] = sort([X.afterIndex] + [X.frac]);
    X = X(ord);
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
