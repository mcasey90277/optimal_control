function T = mesh_report(S, tag, resDir)
% MESH_REPORT  Print and save the mesh-convergence study table for one row.
%
% One line per tracked quantity: its value at each level, the observed order
% where one is admissible, the Richardson limit where permitted, and a verdict.
%
% THE VERDICTS, AND WHAT EACH IS ALLOWED TO CLAIM:
%   CONVERGED       physical quantity, differences monotone and shrinking, and
%                   observed p in [1,3]. The only verdict that asserts the
%                   discrete solution is approaching a continuous limit.
%   NOT-CONVERGED   physical quantity that is non-monotone, or whose p falls
%                   outside [1,3].
%   STABLE          integer quantity unchanged across every level.
%   NOT-STABLE      integer quantity that moved. For a switch count this is a
%                   FINDING about the solution's topology, not a solver fault.
%   O(h)-CONSISTENT a tracked diagnostic shrinking at about the refinement
%                   ratio -- the signature of a representation artifact rather
%                   than a property of the solution.
%   DIAGNOSTIC      tracked, reported, and given NO convergence verdict.
%                   Solver/discretization readouts measure the NLP and its
%                   solution process, not the continuous OCP; ranking them
%                   beside physical quantities would imply otherwise.
%   INSUFFICIENT    fewer than three levels, so no order exists. Two levels
%                   support a stability statement only.
%
% THE HYPOTHESIS COLUMN. For the quantities where the study's two external
% reviewers made competing predictions -- final mass, delta-V, switch movement
% -- the measured p is named against both: H1 (p ~ 1, the switch-crossing
% interval dominates) or H2 (p ~ 2, sub-grid switch placement and
% alternating-sign cancellation). Discriminating between them is the study's
% purpose, so the report states which the measurement supports rather than
% quietly confirming whichever was written down first.
%
% INPUTS:
%   S      - mesh_collect output [struct]
%   tag    - row label for the header and filename, e.g. 'MEE_M2_10N' [char]
%   resDir - directory for mesh_<tag>.mat; '' or omitted to skip saving [char]
%
% OUTPUTS:
%   T - table of the printed rows [table], also saved in the .mat
%
% SIDE EFFECTS: prints the table; writes <resDir>/mesh_<tag>.mat when resDir
% is given. The .mat carries S and T only -- never a live CasADi model, which
% cannot be serialized (MATLAB warns "Serializing SWIG objects not supported").
%
% REFERENCES:
%   [1] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Task 4.

if nargin < 3, resDir = ''; end
names = fieldnames(S.q);
nL = S.nLevels;

fprintf('\n===== MESH-CONVERGENCE STUDY: %s =====\n', tag);
fprintf('levels    : %s  (factors)\n', num2str(S.factors, '%-12g'));
fprintf('intervals : %s\n', num2str(S.N, '%-12d'));
fprintf('h_phys med: %s   <- the relevant h; node-count multipliers are not h\n', ...
        num2str(S.hPhys, '%-12.4g'));
fprintf('h_phys max: %s   (ratio max/med = %.1f on the coarsest level)\n', ...
        num2str(S.hPhysMax, '%-12.4g'), S.hPhysMax(1)/S.hPhys(1));

if nL < 3
    fprintf(['\nNOTE: %d levels. No order can be reported -- a p from two levels ' ...
             'is a bug. Stability statements only.\n'], nL);
end

hdr = sprintf('%-22s %-10s', 'quantity', 'class');
for k = 1:nL, hdr = [hdr sprintf(' %-14s', sprintf('x%g', S.factors(k)))]; end %#ok<AGROW>
hdr = [hdr sprintf(' %-8s %-14s %-16s %-6s', 'p', 'Richardson', 'verdict', 'hyp')];
fprintf('\n%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));

rowName = {};  rowVerdict = {};  rowP = [];  rowRich = [];  rowHyp = {};
for q = 1:numel(names)
    nm = names{q};  Q = S.q.(nm);
    [verdict, hyp] = local_verdict(Q, nL);
    line = sprintf('%-22s %-10s', Q.label, Q.class);
    for k = 1:nL
        v = Q.vals(k);
        if isnan(v), line = [line sprintf(' %-14s', '-')]; %#ok<AGROW>
        else,        line = [line sprintf(' %-14.7g', v)]; end %#ok<AGROW>
    end
    pStr = local_num(Q.order.p);
    if ~Q.richardsonOK, pStr = 'policy'; end
    rStr = local_num(Q.order.rich);
    if ~Q.richardsonOK, rStr = 'forbidden'; end
    line = [line sprintf(' %-8s %-14s %-16s %-6s', pStr, rStr, verdict, hyp)]; %#ok<AGROW>
    fprintf('%s\n', line);
    rowName{end+1} = Q.label;     rowVerdict{end+1} = verdict; %#ok<AGROW>
    rowP(end+1) = Q.order.p;      rowRich(end+1) = Q.order.rich; %#ok<AGROW>
    rowHyp{end+1} = hyp; %#ok<AGROW>
end

% --- switch structure gets its own block: it is not a scalar series --------
st = S.switchStability;
fprintf('\n--- switch structure ---\n');
fprintf('counts per level : %s   -> %s\n', num2str(st.counts, '%-6d'), ...
        local_tf(st.stable, 'STABLE', 'NOT-STABLE (a topology finding)'));
fprintf('matched vs finest: %s\n', num2str(st.matched, '%-6d'));
fprintf('unmatched (self) : %s\n', num2str(st.unmatchedSelf, '%-6d'));
fprintf('unmatched (fine) : %s   <- switches the finest mesh has and this one does not\n', ...
        num2str(st.unmatchedFine, '%-6d'));
fprintf('max |dt| matched : %s\n', num2str(st.maxAbsDt, '%-10.4g'));
fprintf('window (median)  : %s   <- PER-SWITCH, 2x that switch''s own local step\n', ...
        num2str(st.tolMedian, '%-10.4g'));
for k = 1:nL-1
    if ~isempty(st.extraFine{k})
        fprintf('  x%g is missing %d fine switch(es) at t = %s\n', ...
                S.factors(k), numel(st.extraFine{k}), num2str(st.extraFine{k}, '%.5f '));
    end
end
meth = cellfun(@(m) m.method, S.switchMeta, 'UniformOutput', false);
fprintf('extraction       : %s\n', strjoin(meth, ', '));
fracEdge = cellfun(@(m) sum(m.subGridFrac < 1e-9 | m.subGridFrac > 1-1e-9), S.switchMeta);
fprintf('at a bracket edge: %s   <- nonzero means extraction collapsed onto nodes\n', ...
        num2str(fracEdge, '%-6d'));

T = table(rowName(:), rowVerdict(:), rowP(:), rowRich(:), rowHyp(:), ...
    'VariableNames', {'quantity','verdict','p','richardson','hypothesis'});

if ~isempty(resDir)
    if ~isfolder(resDir), mkdir(resDir); end
    f = fullfile(resDir, sprintf('mesh_%s.mat', tag));
    save(f, 'S', 'T', '-v7.3');
    fprintf('\nsaved -> %s\n', f);
end
end

% ---------------------------------------------------------------------------
function [v, hyp] = local_verdict(Q, nL)
% LOCAL_VERDICT  Verdict and hypothesis label for one quantity.
% INPUTS:  Q - one S.q entry; nL - number of levels
% OUTPUTS: v - verdict [char]; hyp - 'H1'|'H2'|'-' [char]
hyp = '-';
if all(isnan(Q.vals))
    v = 'ABSENT';  return
end
switch Q.class
    case 'integer'
        v = local_tf(all(Q.vals == Q.vals(1)), 'STABLE', 'NOT-STABLE');
    case 'diagnostic'
        % An O(h) shrink is the signature of a representation artifact. Judged
        % from the value ratio against the refinement ratio, never from a
        % fitted order -- Richardson is forbidden on these.
        if nL >= 3 && all(Q.vals > 0) && all(diff(Q.vals) < 0)
            rr = Q.vals(1:end-1) ./ Q.vals(2:end);
            if all(abs(rr - Q.order.r) < 0.5*Q.order.r)
                v = 'O(h)-CONSISTENT';  return
            end
        end
        v = 'DIAGNOSTIC';
    otherwise                                   % physical
        if nL < 3 || ~Q.richardsonOK
            v = local_tf(nL >= 3, 'REPORTED', 'INSUFFICIENT');
            if ~Q.richardsonOK, v = 'REPORTED'; end
            return
        end
        if Q.order.monotone && Q.order.p >= 1 && Q.order.p <= 3
            v = 'CONVERGED';
        else
            v = 'NOT-CONVERGED';
        end
        if ~isnan(Q.order.p)
            if abs(Q.order.p - 1) < abs(Q.order.p - 2), hyp = 'H1'; else, hyp = 'H2'; end
        end
end
end

% ---------------------------------------------------------------------------
function s = local_num(x)
% LOCAL_NUM  Compact numeric label, '-' for NaN. INPUTS: x  OUTPUTS: s [char]
if isnan(x), s = '-'; else, s = sprintf('%.4g', x); end
end

% ---------------------------------------------------------------------------
function s = local_tf(c, a, b)
% LOCAL_TF  Pick a label by condition. INPUTS: c,a,b  OUTPUTS: s [char]
if c, s = a; else, s = b; end
end
