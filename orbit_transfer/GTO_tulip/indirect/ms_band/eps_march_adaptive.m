function best = eps_march_adaptive(Zseed, probBase, epsStart, epsFloor, tolR, stateMat)
% EPS_MARCH_ADAPTIVE  Sharp-start smoothing continuation with adaptive steps.
%
% Continuation controller for SHARP dual seeds (native near-bang costates;
% GPT-5.6 review + coordinator directive 2026-07-10): start at a small
% epsStart, and instead of a fixed ladder, choose each next eps
% geometrically (target ratio 1/3) with SWITCH-DISPLACEMENT feedback:
% after every converged step, the S = 0 crossings of the new solution
% (sampled on a fixed 4000-point sigma grid) are compared to the previous
% accepted solution's (the seed's, for the first step). If a solve FAILS
% (after eps_march-style relays) or the crossings move by more than
% DISP_TOL grid nodes / change count, the step is DISCARDED (guard
% discipline) and eps is geometrically bisected toward the last accepted
% eps; a first-step failure moves SHARPER (epsTry/3 — for a sharp seed,
% smaller eps is closer to the seed's native structure; up to 3 such
% reductions). Kill-robustness: state is saved to stateMat after every
% solve; a rerun resumes from it (an external ~70-min watchdog killed
% three prior long runs).
%
% success requires an accepted eps <= 1e-3 (march targets epsFloor).
%
% INPUTS:
%   Zseed    - seed at epsStart [(16M-8)x1]
%   probBase - problem struct with sJ set (epsSmooth overridden per step)
%   epsStart - first eps to solve at [scalar, e.g. 1e-2]
%   epsFloor - march target [scalar, e.g. 1e-4]
%   tolR     - per-solve success threshold [scalar, default 1e-9]
%   stateMat - state-file path for resume [char]
%
% OUTPUTS:
%   best - struct: Z, eps (smallest accepted), resNorm, success
%          (eps <= 1e-3), history (struct array per attempted step:
%          eps, resNorm, converged, relays, accepted, maxDisp)
%
% REFERENCES:
%   [1] Bertrand, Epenoy, OCAM 23(4), 2002 (smoothing continuation).
%   [2] .superpowers/sdd/gpt56_review_S1.md (adaptive fine continuation).

if nargin < 5 || isempty(tolR), tolR = 1e-9; end

MAX_RELAY   = 4;
MAX_LM_ITER = 200;
MAX_SOLVES  = 20;
DISP_TOL    = 100;          % grid nodes (of 4000 over sigf; 2.5%)
NGRID       = 4000;
RATIO       = 1/3;          % nominal geometric step
MIN_PROG    = 0.97;         % abandon if bisection stalls (epsTry/epsAcc)

best = struct('Z', [], 'eps', Inf, 'resNorm', Inf, 'success', false, ...
              'history', struct('eps', {}, 'resNorm', {}, 'converged', {}, ...
                                'relays', {}, 'accepted', {}, 'maxDisp', {}));

sGrid   = linspace(0, probBase.sigf, NGRID);
crossRef = s_crossings(Zseed, probBase, epsStart, sGrid);
fprintf('eps_march_adaptive: seed has %d S-crossings on the reference grid\n', ...
        numel(crossRef));

Zwarm    = Zseed(:);
epsAcc   = Inf;             % last ACCEPTED eps (Inf = none yet)
epsTry   = epsStart;
nSolves  = 0;
nSharper = 0;
if nargin >= 6 && ~isempty(stateMat) && isfile(stateMat)
    S0 = load(stateMat);
    best = S0.best;  Zwarm = S0.Zwarm;  epsAcc = S0.epsAcc;
    epsTry = S0.epsTry;  nSolves = S0.nSolves;  nSharper = S0.nSharper;
    crossRef = S0.crossRef;
    fprintf('eps_march_adaptive RESUME: %d solves done, epsAcc=%.3g, epsTry=%.3g\n', ...
            nSolves, epsAcc, epsTry);
end

while nSolves < MAX_SOLVES
    prob = probBase;  prob.epsSmooth = epsTry;
    Rstart = norm(prob.resFun(Zwarm, prob));
    out    = ms_solve(Zwarm, prob, tolR, MAX_LM_ITER);
    nSolves = nSolves + 1;

    nRelay = 0;
    while ~out.success && out.flag == 0 && out.resNorm <= 0.9*Rstart ...
            && nRelay < MAX_RELAY
        nRelay = nRelay + 1;
        fprintf('adaptive: relay %d/%d at eps=%.3g (||R||=%.3e vs %.3e)\n', ...
                nRelay, MAX_RELAY, epsTry, out.resNorm, Rstart);
        Rstart = out.resNorm;
        out    = ms_solve(out.Z, prob, tolR, MAX_LM_ITER);
        nSolves = nSolves + 1;
    end

    accepted = false;  maxDisp = NaN;
    if out.success
        crossNew = s_crossings(out.Z, probBase, epsTry, sGrid);
        maxDisp  = cross_displacement(crossRef, crossNew);
        accepted = maxDisp <= DISP_TOL;
        fprintf('adaptive: eps=%.3g CONVERGED ||R||=%.3e crossings %d disp %g -> %s\n', ...
                epsTry, out.resNorm, numel(crossNew), maxDisp, ...
                pick(accepted, 'ACCEPT', 'REJECT (structure moved)'));
    end
    best.history(end+1) = struct('eps', epsTry, 'resNorm', out.resNorm, ...
        'converged', out.success, 'relays', nRelay, 'accepted', accepted, ...
        'maxDisp', maxDisp); %#ok<AGROW>

    if accepted
        Zwarm = out.Z;  epsAcc = epsTry;
        crossRef = s_crossings(out.Z, probBase, epsTry, sGrid);
        best.Z = out.Z;  best.eps = epsTry;  best.resNorm = out.resNorm;
        if epsTry <= epsFloor, break; end
        epsTry = max(epsFloor, epsTry*RATIO);
    elseif isinf(epsAcc)
        % first step never accepted: move SHARPER (toward the seed's
        % native near-bang structure), up to 3 times
        nSharper = nSharper + 1;
        if nSharper > 3
            fprintf('adaptive: abandoning — first step unconverged after %d sharper tries\n', ...
                    nSharper - 1);
            break;
        end
        epsTry = epsTry*RATIO;
        fprintf('adaptive: first step failed — trying SHARPER eps=%.3g\n', epsTry);
    else
        epsNew = sqrt(epsAcc*epsTry);              % geometric bisection
        if epsNew/epsAcc > MIN_PROG
            fprintf('adaptive: abandoning — bisection stalled (eps %.3g ~ accepted %.3g)\n', ...
                    epsNew, epsAcc);
            break;
        end
        epsTry = epsNew;
        fprintf('adaptive: bisecting to eps=%.3g\n', epsTry);
    end
    if nargin >= 6 && ~isempty(stateMat)
        save(stateMat, 'best', 'Zwarm', 'epsAcc', 'epsTry', 'nSolves', ...
             'nSharper', 'crossRef');
    end
end
best.success = best.eps <= 1e-3;
end

% -------------------------------------------------------------------------
function idx = s_crossings(Z, probBase, epsK, sGrid)
% S_CROSSINGS  Grid indices of the S = 0 crossings of an MS iterate.
%
% INPUTS:
%   Z        - MS unknowns [(16M-8)x1]
%   probBase - problem struct with sJ set
%   epsK     - smoothing value for the trajectory integration [scalar]
%   sGrid    - fixed sigma sampling grid [1xG]
%
% OUTPUTS:
%   idx - grid indices where sign(S) changes [1xK]
prob = probBase;  prob.epsSmooth = epsK;
traj = sms_traj(Z, prob);
[sU, iu] = unique(traj.sig);
Sg  = interp1(sU, traj.S(iu), sGrid, 'linear', 'extrap');
idx = find(diff(sign(Sg)) ~= 0);
end

function d = cross_displacement(cA, cB)
% CROSS_DISPLACEMENT  Max node shift between two crossing sets (Inf if the
% counts differ — the switching structure changed).
%
% INPUTS:
%   cA - reference crossing grid indices [1xK1]
%   cB - new crossing grid indices [1xK2]
%
% OUTPUTS:
%   d - max |shift| in grid nodes, or Inf when K1 ~= K2 [scalar]
if numel(cA) ~= numel(cB)
    d = Inf;
else
    d = max(abs(cA - cB));
    if isempty(d), d = 0; end
end
end

function s = pick(cond, a, b)
% PICK  Ternary string select.
%
% INPUTS:
%   cond - condition [logical scalar]
%   a    - string returned when cond is true
%   b    - string returned when cond is false
%
% OUTPUTS:
%   s - a or b
if cond, s = a; else, s = b; end
end
