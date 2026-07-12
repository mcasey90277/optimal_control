function summary = verify_direct_pmp(matFile, opts)
% VERIFY_DIRECT_PMP  Independent PMP-consistency certificate for a direct solve.
%
% VERIFIER, not solver (no LM anywhere): checks that a certified direct
% -method bang-bang solution satisfies the continuous first-order PMP
% conditions along the whole trajectory, at the resolution its own mesh
% permits. Costates come from the KKT defect duals via the adjudicated
% midpoint-principled map (SMS_SEED_DUALS mode 'd'); each arc of the
% 16-dim Sundman system SMS_EOM is propagated from the direct solution's
% OWN (state; costate) at the arc joint, with the PMP control law from
% the PROPAGATED costates, and the terminal defect against the direct
% solution at the next joint is reported per block. Structure checks
% compare the dual-implied switching function S to the direct throttle,
% and off-throttle segments are classified as true coasts (reach the
% s < 0.05 bound) vs intermediate DIPS (threshold crossings that never
% reach the bound) so counting artifacts are visible.
%
% CERTIFICATE GATING: the "consistent with a continuous PMP extremal"
% paragraph is printed ONLY when (i) state-block defects pass the 1e-2
% line on all non-adjudicated arcs, (ii) every direct switch is matched
% within nodeTol or explicitly adjudicated, (iii) primer mean <= 0.2 deg,
% (iv) |lamM(sigf)| <= 1e-3. Otherwise: "checks incomplete". Adjudicated
% rows (opts.adjArcs / opts.adjSwitches) must be justified by a dig and
% are listed with the certificate. Expected floor (dual-map table,
% test_sms_dualmap): ~5e-3 on worst (perigee) arcs.
%
% NOTE on normalization: costate-block defects are normalized by the
% GLOBAL max of the block over all joints — arcs where a block is locally
% small look better under this normalization than a local one would show.
%
% INPUTS:
%   matFile - direct-solution .mat: out.X [8x(N+1)], out.U [4x(N+1)],
%             out.lamDef [8xN]; sigma/tauf0/factor top-level
%   opts    - (optional) struct: M arcs [default 40], mode dual map
%             [default 'd'], epsEval smoothing for propagation
%             [default 1e-4], nodeTol switch-match tolerance in node
%             indices [default 1], makeFig [default true], adjArcs arc
%             indices adjudicated by an external dig [default []],
%             adjSwitches switch indices adjudicated [default []]
%
% OUTPUTS:
%   summary - struct: perArc [1xM struct array], defTab [Mx8],
%             worstStateDef, nSwitches, nCrossings, nMatched,
%             matchDistNode, minSepTau, offSegMinS/offSegCoast (dwell
%             classification), primerMeanDeg/primerP95Deg/primerMaxDeg,
%             lamMend, certOK, adjArcs, adjSwitches, tauCr, tauSw,
%             Snode, tauN; saved to verify_pmp_<name>.mat (+ .png figure)
%
% REFERENCES:
%   [1] .superpowers/sdd/gpt56_review_S1.md (dual-map adjudication).
%   [2] MS_BAND_CAMPAIGN.md 2026-07-10 entries (verifier reframe +
%       switch-count adjudication).

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'M'),           opts.M = 40;            end
if ~isfield(opts, 'mode'),        opts.mode = 'd';        end
if ~isfield(opts, 'epsEval'),     opts.epsEval = 1e-4;    end
if ~isfield(opts, 'nodeTol'),     opts.nodeTol = 1;       end
if ~isfield(opts, 'makeFig'),     opts.makeFig = true;    end
if ~isfield(opts, 'adjArcs'),     opts.adjArcs = [];      end
if ~isfield(opts, 'adjSwitches'), opts.adjSwitches = [];  end
[~, baseName] = fileparts(matFile);

[~, prob, info] = sms_seed_duals(matFile, opts.M, opts.epsEval, opts.mode);
M    = opts.M;
tauN = info.tauN;  Y16 = info.Y16;  X = info.X;  U = info.U;
nN   = size(X, 2);
% direct-solution (state; costate) at ALL M+1 joints (source-grid samples)
yJ = interp1(tauN.', Y16.', prob.sJ.', 'linear').';    % [16 x (M+1)]

% global block magnitudes for costate normalization (over all joints; see
% the normalization NOTE in the header)
blk  = {1:3, 4:6, 7, 8, 9:11, 12:14, 15, 16};
bMag = cellfun(@(b) max(max(abs(yJ(b, :)), [], 'all'), 1e-12), blk);

% ---- (2)+(3) per-arc propagation defects + along-arc |Ht + lamT| ----------
perArc = struct('arc', {}, 'def', {}, 'HtMax', {}, 'HtRms', {});
defTab = zeros(M, 8);
for k = 1:M
    [~, Yk] = ode113(@(ss, y) sms_eom(ss, y, prob.Tmax, prob.c, ...
              prob.muStar, opts.epsEval, prob.pSund), ...
              [prob.sJ(k) prob.sJ(k+1)], yJ(:, k), prob.odeOpts);
    d = abs(Yk(end, :).' - yJ(:, k+1));
    def = zeros(1, 8);
    for b = 1:4, def(b) = max(d(blk{b})); end                 % absolute ND
    for b = 5:8, def(b) = max(d(blk{b}))/bMag(b); end         % relative
    HtArc = zeros(1, size(Yk, 1));
    for q = 1:size(Yk, 1)
        [~, Htq] = sms_eom(0, Yk(q, :).', prob.Tmax, prob.c, ...
                           prob.muStar, opts.epsEval, prob.pSund);
        HtArc(q) = Htq + Yk(q, 16);
    end
    perArc(k) = struct('arc', k, 'def', def, 'HtMax', max(abs(HtArc)), ...
                       'HtRms', rms(HtArc));
    defTab(k, :) = def;
end

% ---- (4) switching-structure check on the FULL node grid ------------------
lamVn = Y16(12:14, :);  lamMn = Y16(15, :);
Snode = 1 - sqrt(sum(lamVn.^2, 1))*prob.c./X(7, :) - lamMn;
crossI = find(diff(sign(Snode)) ~= 0);                 % node index before crossing
tauCr  = zeros(1, numel(crossI));
for q = 1:numel(crossI)
    kq = crossI(q);
    tauCr(q) = tauN(kq) + (0 - Snode(kq))*(tauN(kq+1) - tauN(kq)) ...
               /(Snode(kq+1) - Snode(kq));
end
s     = U(4, :);
swI   = find(diff(double(s > 0.5)) ~= 0);              % direct switch intervals
tauSw = (tauN(swI) + tauN(swI+1))/2;
matched = false(1, numel(swI));  matchDistNode = inf(1, numel(swI));
for q = 1:numel(swI)
    if isempty(crossI), break; end
    [dNode, ~] = min(abs(crossI - swI(q)));
    matchDistNode(q) = dNode;
    matched(q) = dNode <= opts.nodeTol;
end
% resolution dig: separations between adjacent S crossings
if numel(tauCr) >= 2
    minSepTau = min(diff(sort(tauCr)));
else
    minSepTau = Inf;
end
% coast-bound dwell: classify each off-throttle segment (s <= 0.5 between
% threshold crossings) as a TRUE COAST (reaches s < 0.05) or a DIP
segStart = [1, swI + 1];
segEnd   = [swI, nN];
segOn    = arrayfun(@(a, b) mean(s(a:b)) > 0.5, segStart, segEnd);
offSeg   = find(~segOn);
offSegMinS = arrayfun(@(q) min(s(segStart(q):segEnd(q))), offSeg);
offSegCoast = offSegMinS < 0.05;

% ---- (5) primer alignment on burn arcs ------------------------------------
burnN = find(s > 0.5);
% exclude nodes within 3 of a switch
for w = -3:3, burnN = setdiff(burnN, swI + w); end
primer = -lamVn(:, burnN)./sqrt(sum(lamVn(:, burnN).^2, 1));
aDir   = U(1:3, burnN)./sqrt(sum(U(1:3, burnN).^2, 1));
cang   = min(max(sum(primer.*aDir, 1), -1), 1);
angDeg = acosd(cang);
primerMean = mean(angDeg);
primerP95  = prctile(angDeg, 95);
primerMax  = max(angDeg);

% ---- (6) transversality ----------------------------------------------------
lamMend = abs(Y16(15, end));

% ---- print -----------------------------------------------------------------
fprintf('\n=== verify_direct_pmp: %s  (mode %s, M = %d, epsEval = %.0e) ===\n', ...
        baseName, opts.mode, M, opts.epsEval);
fprintf('beta = %.5f  node1Err = %.3e\n', info.beta, info.node1Err);
fprintf('%-4s %-9s %-9s %-9s %-9s | %-9s %-9s %-9s %-9s | %-9s %s\n', ...
        'arc', 'r', 'v', 'm', 't', 'lamR%', 'lamV%', 'lamM%', 'lamT%', ...
        '|Ht+lamT|', 'flag');
for k = 1:M
    d = defTab(k, :);
    flag = 'PASS';
    if max(d(1:4)) > 1e-2
        flag = 'ATTN';
        if ismember(k, opts.adjArcs), flag = 'ADJ'; end
    end
    fprintf('%-4d %-9.2e %-9.2e %-9.2e %-9.2e | %-9.2e %-9.2e %-9.2e %-9.2e | %-9.2e %s\n', ...
            k, d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), ...
            perArc(k).HtMax, flag);
end
worstState = max(max(defTab(:, 1:4)));
nonAdjArcs = setdiff(1:M, opts.adjArcs);
stateOK    = all(max(defTab(nonAdjArcs, 1:4), [], 2) <= 1e-2);
swAdj      = ismember(1:numel(swI), opts.adjSwitches);
swOK       = all(matched | swAdj);
primerOK   = primerMean <= 0.2;
transOK    = lamMend <= 1e-3;
certOK     = stateOK && swOK && primerOK && transOK;

fprintf('\nSummary (%s):\n', baseName);
fprintf('  worst state-block arc defect      : %.3e all arcs; %.3e non-adjudicated (line 1e-2, floor ~5e-3)  %s\n', ...
        worstState, max(max(defTab(nonAdjArcs, 1:4))), pass_attn(stateOK));
fprintf('  direct switches (s>0.5 crossings) : %d;  S-crossings on node grid: %d\n', ...
        numel(swI), numel(crossI));
fprintf('  switches matched (<= %d node)      : %d/%d   (max node dist %g)  %s\n', ...
        opts.nodeTol, nnz(matched), numel(swI), max(matchDistNode), ...
        pass_attn(swOK));
for q = find(~matched)
    win = max(1, swI(q)-20):min(nN, swI(q)+20);
    fprintf(['    UNMATCHED switch #%d at node %d (tau %.3f, t %.4f): dual-S ' ...
             'stays one-signed, S@switch %.3e, min|S| within +-20 nodes ' ...
             '%.3e%s\n'], q, swI(q), tauN(swI(q)), X(8, swI(q)), ...
            Snode(swI(q)), min(abs(Snode(win))), ...
            pick(swAdj(q), '  [ADJUDICATED]', ''));
end
fprintf('  off-throttle segments             : %d, of which %d reach the coast bound (s<0.05); %d dip(s), min s = [%s]\n', ...
        numel(offSeg), nnz(offSegCoast), nnz(~offSegCoast), ...
        sprintf('%.2f ', offSegMinS(~offSegCoast)));
fprintf('  min separation of S-crossings     : %.4f tau units (grid-resolution dig)\n', minSepTau);
fprintf('  primer alignment on burn nodes    : mean %.4f deg, p95 %.4f deg, max %.4f deg   %s\n', ...
        primerMean, primerP95, primerMax, pass_attn(primerOK));
fprintf('  |lamM(sigma_f)| (transversality)  : %.3e   %s\n', ...
        lamMend, pass_attn(transOK));
if certOK
    fprintf(['CERTIFICATE: consistent with a continuous PMP extremal at ' ...
             'the transcription''s\nO(h^2) resolution (dual-map floor ' ...
             '~5e-3 on worst perigee arcs); NOT a\nshooting-converged ' ...
             'extremal.\n']);
    if ~isempty(opts.adjArcs) || ~isempty(opts.adjSwitches)
        fprintf('  with adjudications: arcs [%s], switches [%s] (justifications in the run record)\n', ...
                sprintf('%d ', opts.adjArcs), sprintf('%d ', opts.adjSwitches));
    end
else
    fprintf('CERTIFICATE NOT ISSUED: checks incomplete — see ATTN rows.\n');
end

summary = struct('matFile', matFile, 'opts', opts, 'beta', info.beta, ...
    'perArc', perArc, 'defTab', defTab, 'worstStateDef', worstState, ...
    'nSwitches', numel(swI), 'nCrossings', numel(crossI), ...
    'nMatched', nnz(matched), 'matchDistNode', matchDistNode, ...
    'minSepTau', minSepTau, 'offSegMinS', offSegMinS, ...
    'offSegCoast', offSegCoast, 'primerMeanDeg', primerMean, ...
    'primerP95Deg', primerP95, 'primerMaxDeg', primerMax, ...
    'lamMend', lamMend, 'certOK', certOK, 'adjArcs', opts.adjArcs, ...
    'adjSwitches', opts.adjSwitches, 'tauCr', tauCr, 'tauSw', tauSw, ...
    'Snode', Snode, 'tauN', tauN);
save(sprintf('verify_pmp_%s.mat', baseName), 'summary');

% ---- figure -----------------------------------------------------------------
if opts.makeFig
    fig = figure('Visible', 'off', 'Position', [50 50 1100 720]);
    subplot(2, 1, 1);
    semilogy(1:M, max(defTab(:, 1:4), [], 2), 'o-', 'LineWidth', 1.2); hold on;
    semilogy(1:M, max(defTab(:, 5:8), [], 2), 's-', 'LineWidth', 1.2);
    yline(1e-2, 'k--'); yline(5e-3, 'k:');
    grid on; xlabel('arc'); ylabel('max defect');
    legend('state blocks (abs ND)', 'costate blocks (rel)', ...
           'heuristic 1e-2', 'dual-map floor 5e-3', 'Location', 'best');
    title(sprintf('%s: per-arc PMP propagation defects (mode %s, M = %d)', ...
          baseName, opts.mode, M), 'Interpreter', 'none');
    subplot(2, 1, 2);
    tPhys = X(8, :);
    yyaxis left;  plot(tPhys, Snode, 'LineWidth', 1.0); ylabel('S (from duals)');
    hold on; plot(interp1(tauN, tPhys, tauCr), zeros(size(tauCr)), 'go');
    yline(0, 'k-');
    yyaxis right; stairs(tPhys, s, 'LineWidth', 0.8); ylabel('direct throttle s');
    ylim([-0.05 1.05]);
    xlabel('t (ND)'); grid on;
    title('switching function vs direct throttle (o = S zero crossings)');
    exportgraphics(fig, sprintf('verify_pmp_%s.png', baseName), 'Resolution', 140);
    close(fig);
end
end

% -------------------------------------------------------------------------
function s = pass_attn(cond)
% PASS_ATTN  PASS/ATTN string for verdict rows.
%
% INPUTS:
%   cond - condition [logical scalar]
%
% OUTPUTS:
%   s - 'PASS' or 'ATTN'
if cond, s = 'PASS'; else, s = 'ATTN'; end
end

function out = pick(cond, a, b)
% PICK  Ternary string select.
%
% INPUTS:
%   cond - condition [logical scalar]
%   a    - string returned when cond is true
%   b    - string returned when cond is false
%
% OUTPUTS:
%   out - a or b
if cond, out = a; else, out = b; end
end
