function res = run_mintime(thrustN, hx0, N)
% RUN_MINTIME  Free-L (manifold) min-time anchor at one thrust level.
%
% t_f,min sets every c_tf scale downstream (paper's TfMin is free-longitude).
%
% METHODOLOGICAL FINDING (Task 8, 2026-07-16): a direct one-shot 'mintime'
% solve of casadi_lt_2body against the free-longitude insertion MANIFOLD,
% warm-started straight from the cold tangential seed_2body seed, does NOT
% converge (tested coplanar case: 3000 iters, Maximum_Iterations_Exceeded,
% maxDefect ~2.6e-3, termErr ~0.14 -- gates missed by 5+ orders of magnitude).
% The manifold terminal set is too weakly constraining relative to the
% free-time mintime objective for the cold seed to find the right basin.
%
% What DOES converge (used here as the real algorithm, not a fallback):
%   Stage 1 - 'mintime' vs a FIXED terminal rendezvous at the seed's own
%             arrival longitude (geo_terminal('fixed', p, sinfo.Larr)) --
%             zero terminal gap by construction, so this is an easy warm-up
%             even if IPOPT does not fully converge it (max-iter here is
%             fine; only the trajectory SHAPE is used downstream).
%   Stage 2 - 'mintime' vs the free-longitude MANIFOLD, warm-started from
%             stage 1's (possibly non-converged) trajectory, warmTight=true.
%   Continuation - if stage 2 itself does not reach gate tolerance, repeat
%             the manifold solve warm-started from the previous iterate
%             (warmTight=true) up to CONTINGENCY_MAX_ROUNDS times. Empirically
%             this needed 0 extra rounds for the 3D case (stage 2 converged
%             directly) and 1 extra round for the coplanar case (stage 2
%             stalled at maxDefect 2.9e-4, one more warm continuation reached
%             2.6e-12). A round that improves maxDefect by less than one
%             decade without reaching the gate tolerance is treated as a
%             stall and raises an error (never grind indefinitely).
%
% Caching is CONVERGED-ONLY: a result is written to the canonical cache path
% only if success && maxDefect<1e-8 && termErr<1e-8, so a non-converged
% iterate can never sit in a path that downstream tasks silently load.
%
% INPUTS:  thrustN - max thrust [N];  hx0 - initial hx (0 coplanar | 0.0612);
%          N - mesh segments (default 600)
% OUTPUTS: res - .out .tfmin .tfmin_h .dL_mt .revs (also saved/cached in
%          results/); .continuationRounds records how many extra warm
%          continuation rounds beyond stage 1+2 were needed.
%
% REFERENCES: [1] DESIGN.md sec 4 step 1.  [2] PLAN.md Task 8.
%   [3] task-8-brief.md contingency comment (origin of the stage-1/stage-2
%       idea; promoted here from documented fallback to the default path
%       after the one-shot manifold solve was shown not to converge).
if nargin < 3, N = 600; end
CONTINGENCY_MAX_ROUNDS = 3;

here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
tag = sprintf('mintime_T%d_i%d', round(10*thrustN), round(hx0 > 0)*7);
fn  = fullfile(resDir, [tag '.mat']);
if isfile(fn), S = load(fn); res = S.res; fprintf('cached %s\n', fn); return; end

p  = kepler_lt_params(thrustN, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, hx0, 0, pi, p.mu);
rv0 = [r0; v0];
[sg, X0, U0, tauf0, sinfo] = seed_2body(p, rv0, struct('sbar',1,'tDur',[],'N',N));
term_fixed    = geo_terminal('fixed', p, sinfo.Larr);
term_manifold = geo_terminal('manifold', p, []);
baseOpts = struct('par',p, 'mode','mintime', 'rv0',rv0, 'maxIter',3000, 'printLevel',3);
warmOpts = baseOpts;  warmOpts.warmTight = true;
isGood = @(o) o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;

fprintf('MINTIME stage 1 (fixed rendezvous @ seed arrival L=%.4f)...\n', sinfo.Larr);
out1 = casadi_lt_2body(sg, X0, U0, tauf0, term_fixed, baseOpts);
fprintf('  stage1: status=%s defect=%.3e termErr=%.3e (max-iter here is fine -- warm-up only)\n', ...
        out1.ipoptStatus, out1.maxDefect, out1.termErr);

fprintf('MINTIME stage 2 (manifold, warm-started from stage 1)...\n');
out = casadi_lt_2body(sg, out1.X, out1.U, tauf0, term_manifold, warmOpts);
fprintf('  stage2: status=%s defect=%.3e termErr=%.3e\n', out.ipoptStatus, out.maxDefect, out.termErr);

round_ = 0;
while ~isGood(out) && round_ < CONTINGENCY_MAX_ROUNDS
    round_ = round_ + 1;
    prevDefect = out.maxDefect;
    outNew = casadi_lt_2body(sg, out.X, out.U, tauf0, term_manifold, warmOpts);
    fprintf('  continuation round %d: defect %.3e -> %.3e, termErr=%.3e, status=%s\n', ...
            round_, prevDefect, outNew.maxDefect, outNew.termErr, outNew.ipoptStatus);
    if ~isGood(outNew)
        decadeImprove = log10(max(prevDefect, realmin)) - log10(max(outNew.maxDefect, realmin));
        if decadeImprove < 1
            error('run_mintime:stall', ['continuation stalled at round %d for thrustN=%g hx0=%g: ' ...
                'defect %.3e -> %.3e (%.2f decades, need >=1), termErr=%.3e, status=%s'], ...
                round_, thrustN, hx0, prevDefect, outNew.maxDefect, decadeImprove, outNew.termErr, outNew.ipoptStatus);
        end
    end
    out = outNew;
end
if ~isGood(out)
    error('run_mintime:noConverge', ...
        ['thrustN=%g hx0=%g failed to converge after stage1+stage2+%d continuation round(s): ' ...
         'defect=%.3e termErr=%.3e status=%s'], ...
        thrustN, hx0, round_, out.maxDefect, out.termErr, out.ipoptStatus);
end

Lun = unwrap(atan2(out.X(2,:), out.X(1,:)));
res = struct('out', out, 'tfmin', out.tf, 'tfmin_h', out.tf*p.TU_s/3600, ...
             'dL_mt', Lun(end)-Lun(1), 'revs', (Lun(end)-Lun(1))/(2*pi), ...
             'thrustN', thrustN, 'hx0', hx0, 'N', N, 'seedInfo', sinfo, ...
             'continuationRounds', round_);
save(fn, 'res');
fprintf('MINTIME T=%g N: tf=%.4f ND = %.1f h, dL=%.1f rad (%.2f revs), defect %.2e, %s\n', ...
        thrustN, res.tfmin, res.tfmin_h, res.dL_mt, res.revs, out.maxDefect, out.ipoptStatus);
end
