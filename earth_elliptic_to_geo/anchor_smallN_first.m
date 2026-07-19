function anchor = anchor_smallN_first(T, par, warmAnchor, aopts)
% ANCHOR_SMALLN_FIRST  1 N min-time anchor strategy (Task 3): certify a
% min-time anchor at a LOW nodes-per-rev density first (nprLo), then
% mesh-refine it up to production density (nprHi) with a single tight
% solve. Harvested VERBATIM (not reinvented) from two real campaign hand
% scripts:
%
%   Stage 1 (nprLo) - results/task7c_step1_manual.m's manual relaxed-stall
%     continuation: chain casadi_lt_mee(...,'mode','mintime',...) calls,
%     choosing warmTight/maxIter each round from the INCOMING point's own
%     measured status/defect (smallN_warmtight.m / smallN_maxiter.m, the
%     two pure helpers this loop calls -- see those files for the exact
%     harvested predicates), with NO decadeMin stall-floor enforcement:
%     any forward progress is retained (out = outNew, always advance the
%     warm start), the best-ever point is tracked for diagnostics, and
%     every round is cached to its own file so a crash loses at most one
%     round. Stops when isGood(out) certifies or a round/wall budget is hit.
%   Stage 2 (nprLo -> nprHi) - results/task7c_step1b_refine.m's mesh
%     refine: interp_warmstart the certified nprLo anchor onto the nprHi
%     grid, then ONE warmTight=true, maxIter=1500 casadi_lt_mee solve (no
%     homotopy re-sweep -- the nprLo point is already essentially exact,
%     so the denser mesh only needs to re-satisfy tighter defect/terminal
%     tolerances from an excellent starting shape).
%
% Stage 1's seed is either interp_warmstart of the PREVIOUS rung's
% converged min-time anchor (warmAnchor, a self-similar-shape start, far
% better conditioned than a cold spiral per the campaign) onto the nprLo
% grid, or -- if no warmAnchor is available -- a cold mee_seed tangential
% seed at the SAME throttle/rev-count defaults run_mintime_mee.m's own
% Stage B cold fallback uses (thr=0.4, nRev=3), so the fallback shape is
% never a fresh invention either.
%
% INPUTS:
%   T          - max thrust level [N]                                  [scalar]
%   par        - kepler_lt_params(T, m0kg, ispS) struct (dimensional
%                context: TU_s, Tmax, ...); caller-supplied so this
%                function never re-derives m0kg/ispS on its own          [struct]
%   warmAnchor - previous rung's CONVERGED min-time anchor, fields
%                .solverOut.X [7x(N+1)] .solverOut.U [4x(N+1)] .dL_mt
%                [scalar] .N [scalar] (exactly run_mintime_mee's/this
%                function's own output shape), used to seed Stage 1 via
%                interp_warmstart; pass [] to fall back to a cold mee_seed
%                tangential seed instead                            [struct|[]]
%   aopts      - optional struct:
%     .nprLo       - Stage 1 (low-density) nodes-per-rev; default 15      [scalar]
%     .nprHi       - Stage 2 (production-density) nodes-per-rev; default 25 [scalar]
%     .tag         - REPRO_-family cache tag stem for Stage 1's per-round
%                    cache files (results/repro/<tag>_smallN_round%02d.mat);
%                    default built from T (see local ttag_local)           [char]
%     .maxRoundsManual - Stage 1 continuation round budget; default 60     [scalar]
%     .maxWallS        - Stage 1 continuation wall-clock budget [s];
%                    default 2.5*3600                                      [scalar]
%
% OUTPUTS:
%   anchor - struct MATCHING run_mintime_mee's output shape, so
%            reproduce_row.m's fuel stage consumes it uniformly:
%              .tfmin        - converged min transfer time [ND]           [scalar]
%              .tfmin_h      - tfmin in hours                              [scalar]
%              .dL_mt        - converged total true-longitude span [rad]   [scalar]
%              .revs         - dL_mt/(2*pi)                                [scalar]
%              .N            - Stage 2 (nprHi) node count                  [scalar]
%              .solverOut    - full casadi_lt_mee output at nprHi          [struct]
%              .certified    - true (Stage 2 assert already enforced this) [logical]
%              .anchorSource - 'smallN_first'                              [char]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/results/task7c_step1_manual.m (Stage 1
%       source, harvested verbatim).
%   [2] earth_elliptic_to_geo/results/task7c_step1b_refine.m (Stage 2
%       source, harvested verbatim).
%   [3] earth_elliptic_to_geo/run_mintime_mee.m (Stage-B cold-seed
%       defaults thr=0.4/nRev=3 mirrored by this function's cold fallback;
%       the output-struct shape this function matches).
%   [4] earth_elliptic_to_geo/interp_warmstart.m / casadi_lt_mee.m /
%       mee_seed.m (composed, not edited, by this function).
%   [5] .superpowers/sdd/task-3-brief.md (this task's spec).

if nargin < 3, warmAnchor = []; end
if nargin < 4, aopts = struct(); end
d = @(f,v) optdef(aopts, f, v);

nprLo           = d('nprLo', 15);
nprHi           = d('nprHi', 25);
tag             = d('tag', sprintf('REPRO_smallN_%s', ttag_local(T)));
maxRoundsManual = d('maxRoundsManual', 60);
maxWallS        = d('maxWallS', 2.5*3600);

here     = fileparts(mfilename('fullpath'));
reproDir = fullfile(here, 'results', 'repro');
if ~exist(reproDir, 'dir'), mkdir(reproDir); end

nRevSeed = 3;     % run_mintime_mee.m's Stage-B cold-seed default nRev
seedThr  = 0.4;   % run_mintime_mee.m's Stage-B cold-seed default throttle

% ---------------------------------------------------------------------------
% STAGE 1 SEED (nprLo grid)
% ---------------------------------------------------------------------------
if ~isempty(warmAnchor)
    revsGuess  = warmAnchor.dL_mt / (2*pi);
    N1         = round(nprLo * revsGuess);
    sigmaPrev  = linspace(0, 1, warmAnchor.N + 1).';
    sigma1     = linspace(0, 1, N1 + 1).';
    W  = interp_warmstart(warmAnchor.solverOut.X, warmAnchor.solverOut.U, ...
        warmAnchor.dL_mt, sigmaPrev, sigma1);
    X1 = W.X;  U1 = W.U;  dL1 = W.dL;
    fprintf(['ANCHOR_SMALLN_FIRST T=%g N Stage 1 seed: warm-started from prior ' ...
             'rung''s anchor (revs_guess=%.4f, N=%d nodes, %d nodes/rev)\n'], ...
            T, revsGuess, N1, nprLo);
else
    N1 = round(nprLo * nRevSeed);
    [sigma1, X1, U1, dL1, seedInfo] = mee_seed(par, struct('thr', seedThr, ...
        'betaMode', 'tangential', 'nRev', nRevSeed, 'N', N1));
    fprintf(['ANCHOR_SMALLN_FIRST T=%g N Stage 1 seed: cold tangential thr=%.2g, ' ...
             'nRev=%g target (N=%d, %d nodes/rev, achieved seed nRev=%.4f)\n'], ...
            T, seedThr, nRevSeed, N1, nprLo, seedInfo.nRev);
end
x0state = X1(:,1);

% ---------------------------------------------------------------------------
% STAGE 1: manual relaxed-stall continuation (harvested from
% results/task7c_step1_manual.m)
% ---------------------------------------------------------------------------
isGood = @(o) o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;

% Round 0: bridge the raw ODE seed (no ipoptStatus/maxDefect of its own)
% into a casadi_lt_mee output so the harvested continuation loop below has
% an incoming point to evaluate smallN_warmtight/smallN_maxiter against.
% task7c_step1_manual.m itself had no equivalent bridging step because it
% resumed from an ALREADY-SOLVED round produced by a separate driver; here
% the loose/adaptive-mu convention (warmTight=false) and smallN_maxiter's
% own no-prior-defect default (Inf -> 75) are reused rather than inventing
% an unharvested maxIter constant.
round0File = fullfile(reproDir, sprintf('%s_smallN_round00.mat', tag));
if isfile(round0File)
    S0 = load(round0File);  out = S0.out;
    fprintf('  [cached] %s round 0\n', tag);
else
    out = casadi_lt_mee(sigma1, X1, U1, dL1, struct('par', par, 'mode', 'mintime', ...
        'x0', x0state, 'maxIter', smallN_maxiter(Inf), 'warmTight', false, 'printLevel', 3));
    save(round0File, 'out');
end
fprintf('  round 0: status=%s defect=%.3e termErr=%.3e\n', ...
        out.ipoptStatus, out.maxDefect, out.termErr);

rnd  = 0;
best = out;
t0   = tic;
while ~isGood(out) && rnd < maxRoundsManual && toc(t0) < maxWallS
    rnd = rnd + 1;
    rf  = fullfile(reproDir, sprintf('%s_smallN_round%02d.mat', tag, rnd));
    if isfile(rf)
        S = load(rf);  outNew = S.outNew;
        fprintf('  [cached] round %d\n', rnd);
    else
        warmTight = smallN_warmtight(out.ipoptStatus, out.maxDefect);
        maxIterR  = smallN_maxiter(out.maxDefect);
        rt0 = tic;
        outNew = casadi_lt_mee(linspace(0,1,size(out.X,2)).', out.X, out.U, out.dL, ...
            struct('par', par, 'mode', 'mintime', 'x0', x0state, 'maxIter', maxIterR, ...
            'warmTight', warmTight, 'printLevel', 3));
        rtElapsed = toc(rt0);
        save(rf, 'outNew', 'rtElapsed', 'warmTight', 'maxIterR');
    end
    fprintf(['  round %d (warmTight=%d, maxIter=%d): defect %.4e -> %.4e, termErr=%.4e, ' ...
             'status=%s (%.1f s elapsed total)\n'], rnd, warmTight, maxIterR, ...
            out.maxDefect, outNew.maxDefect, outNew.termErr, outNew.ipoptStatus, toc(t0));
    if outNew.maxDefect < best.maxDefect
        best = outNew;
    end
    out = outNew;   % always advance the warm-start chain (no floor gate -- see header)
end

assert(isGood(out), 'anchor_smallN_first:stage1NotCertified', ...
    ['T=%g N: Stage 1 (nprLo=%d) manual continuation did not certify after %d ' ...
     'round(s), %.1f s wall (best-ever defect=%.3e termErr=%.3e status=%s)'], ...
    T, nprLo, rnd, toc(t0), best.maxDefect, best.termErr, best.ipoptStatus);

anchorLo = out;
fprintf(['ANCHOR_SMALLN_FIRST T=%g N Stage 1 CERTIFIED: nprLo=%d, tf=%.6f ND, ' ...
         'revs=%.4f, defect=%.3e, %d round(s)\n'], ...
        T, nprLo, anchorLo.tf, anchorLo.dL/(2*pi), anchorLo.maxDefect, rnd);

% ---------------------------------------------------------------------------
% STAGE 2: mesh-refine nprLo -> nprHi (harvested from
% results/task7c_step1b_refine.m)
% ---------------------------------------------------------------------------
revsHi     = anchorLo.dL / (2*pi);
Nhi        = round(nprHi * revsHi);
sigmaSrcHi = linspace(0, 1, size(anchorLo.X,2)).';
sigmaDstHi = linspace(0, 1, Nhi + 1).';
Whi = interp_warmstart(anchorLo.X, anchorLo.U, anchorLo.dL, sigmaSrcHi, sigmaDstHi);
x0hi = Whi.X(:,1);

fprintf('ANCHOR_SMALLN_FIRST T=%g N Stage 2 refine: nprLo=%d (N=%d) -> nprHi=%d (N=%d), revs=%.4f\n', ...
        T, nprLo, size(anchorLo.X,2)-1, nprHi, Nhi, revsHi);

o = casadi_lt_mee(sigmaDstHi, Whi.X, Whi.U, Whi.dL, struct('par', par, 'mode', 'mintime', ...
    'x0', x0hi, 'maxIter', 1500, 'warmTight', true, 'printLevel', 3));
fprintf('  Stage 2 refine solve: status=%s defect=%.4e termErr=%.4e tf=%.6f revs=%.4f\n', ...
        o.ipoptStatus, o.maxDefect, o.termErr, o.tf, o.dL/(2*pi));

certified = o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;
assert(certified, 'anchor_smallN_first:stage2NotCertified', ...
    ['T=%g N: Stage 2 refine (nprLo=%d -> nprHi=%d, N=%d) did NOT certify: ' ...
     'status=%s defect=%.3e termErr=%.3e'], T, nprLo, nprHi, Nhi, o.ipoptStatus, o.maxDefect, o.termErr);

anchor = struct('tfmin', o.tf, 'tfmin_h', o.tf*par.TU_s/3600, 'dL_mt', o.dL, ...
    'revs', o.dL/(2*pi), 'N', Nhi, 'solverOut', o, 'certified', true, ...
    'anchorSource', 'smallN_first');

fprintf(['ANCHOR_SMALLN_FIRST T=%g N DONE: tfmin=%.6f ND (%.2f h), revs=%.4f, ' ...
         'N=%d (%d nodes/rev), defect=%.3e\n'], ...
        T, anchor.tfmin, anchor.tfmin_h, anchor.revs, anchor.N, nprHi, o.maxDefect);

end

% ---------------------------------------------------------------------------
function s = ttag_local(T)
% TTAG_LOCAL  Numeric-only tag fragment for a thrust level (mirrors
% reproduce_row.m's own ttag / mee_fuel_tag.m's convention): integer -> plain
% digits, non-integer -> decimal point replaced with 'p'. Used only when the
% caller does not supply aopts.tag explicitly.
%
% INPUTS:  T - max thrust [N]                                       [scalar]
% OUTPUTS: s - numeric tag fragment, e.g. 1->'1', 2.5->'2p5'          [char]
if abs(T - round(T)) < 1e-9
    s = sprintf('%d', round(T));
else
    s = strrep(sprintf('%g', T), '.', 'p');
end
end
