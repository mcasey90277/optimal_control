function foc_report(rep, tag, resDir)
% FOC_REPORT  Standard first-order optimality report (print + optional sidecar save).
%
% Prints a fixed-format FOC/KKT advisory report from foc_check output struct.
% This is the standard reporting block for all campaign drivers (deliverable b).
% All verdicts are advisory (report-only burn-in); does not alter certified status.
%
% INPUTS:
%   rep   - report struct from foc_check [struct], fields:
%           .kktStatInf       max-norm Lagrangian stationarity residual [scalar]
%           .sLag             resolved sign convention (+1 or -1) [scalar]
%           .dirTanMax        max tangential dL/dbeta residual [scalar]
%           .dirTanMed        median tangential dL/dbeta residual [scalar]
%           .signPct          switching-function sign-law agreement [%] [scalar]
%           .Sd               switching function samples [1 x (N+1)] [numeric]
%           .lam              nodal costates [nx x (N+1)] [numeric]
%           .lamTimeCoV       time-costate coefficient of variation [scalar]
%           .lamTimeEnd       time-costate value at final node [scalar]
%           .lamMassEndMapped free-mass transversality residual (relative),
%                             Hager-mapped terminal covector [scalar]
%           .derivedFromKKT   names of reported lines that are PROJECTIONS of
%                             kktStatInf, not independent tests [cellstr]
%           .singularArcNodes count of >=3-node near-zero switching runs [scalar]
%           .sdotMinRel       minimum relative |Sdot| at switches [scalar]
%           .nSwitches        number of burn/coast sign changes [scalar]
%           .horizonNote      horizon-condition caveat text [char]
%           .checksRun        which checks executed [cellstr]
%           .pass             advisory overall verdict [logical]
%           (optional) .ipopt struct with .verdict [char] for 2nd-order check
%   tag    - file basename for sidecar naming [char or string]
%   resDir - results directory for optional sidecar save [char or string];
%            if empty or omitted, only prints (no sidecar)
%
% OUTPUTS:
%   (none) Prints report block to stdout; optionally saves sidecar .mat.
%
% REFERENCES:
%   [1] orbit_transfer/verify_common/foc_check.m — FOC gate producing rep struct.
%   [2] .superpowers/sdd/2026-07-25-foc-gate-layer/task-4-brief.md — report format.
%   [3] Hager, W.W., "Runge-Kutta Methods in Optimal Control and the
%       Transformed Adjoint System," Numer. Math. 87, 247-282, 2000 (the
%       mapped terminal covector printed as the gated transversality line).

% Handle defaults
if nargin < 3, resDir = ''; end
if nargin < 2
    error('foc_report: tag and rep required');
end

% Assert tag is char/string
assert(ischar(tag) || isstring(tag), ...
    'foc_report: tag must be char or string');
tag = char(tag);  % ensure char

% Thresholds (must match foc_check defaults)
tolStat = 1e-6;   % KKT stationarity / tangential tolerance
tolSign = 99;     % sign-law pass threshold [%]
tolTrans = 1e-3;  % free-mass transversality relative tolerance
sdotMin = 1e-3;   % minimum regular-switching relative Sdot

% --- Header ----------------------------------------------------------------
fprintf('\n========== FIRST-ORDER OPTIMALITY REPORT: %s ==========\n', tag);

% --- MASTER RESIDUAL, and the two lines that are projections of it ----------
% Grouping matters and is not cosmetic. The tangential direction residual and
% the mapped terminal covector are, by construction, a projection and a single
% ENTRY of the same vector whose max-norm is kktStatInf (see foc_check.m
% sections (3) and (6)). Neither can fail unless kktStat already has, so
% printing them flush with the structural checks would overstate how much
% independent evidence this report carries. They are indented under it.
s_kkt = make_status_str(rep.kktStatInf, tolStat, '<=');
fprintf(' KKT stationarity  ||grad_x L||_inf : %10.3e   (sign s=%+d)   %s\n', ...
    rep.kktStatInf, rep.sLag, s_kkt);

s_dir = make_status_str(rep.dirTanMax, tolStat, '<=');
fprintf('   |- min condition, direction tan  : %10.3e                 %s\n', ...
    rep.dirTanMax, s_dir);

s_trans = '--';
if ~isnan(rep.lamMassEndMapped)
    s_trans = make_status_str(rep.lamMassEndMapped, tolTrans, '<=');
end
fprintf('   |- transversality |lam_m(tf)|    : %10.3e                 %s\n', ...
    rep.lamMassEndMapped, s_trans);

% --- STRUCTURAL CHECKS (genuinely independent of the master residual) -------
if isnan(rep.signPct)
    s_sign = '--';
else
    s_sign = make_status_str(rep.signPct, tolSign, '>=');
end
fprintf(' Bang-bang sign law (S<0 <=> burn)  : %9.2f %%                %s\n', ...
    rep.signPct, s_sign);

% --- INFORMATIONAL (reported, never gated) ---------------------------------
if isnan(rep.lamTimeCoV)
    fprintf(' Time-costate CoV (H-const dual)    :         --   end %+.3e\n', ...
        rep.lamTimeEnd);
else
    fprintf(' Time-costate CoV (H-const dual)    : %10.3e   end %+.3e\n', ...
        rep.lamTimeCoV, rep.lamTimeEnd);
end

% --- Singular arcs ---------------------------------------------------------
if isnan(rep.singularArcNodes)
    s_sarc = '--';
else
    s_sarc = make_count_status_str(rep.singularArcNodes, 0);
end
fprintf(' Singular-arc nodes                 : %10d                 %s\n', ...
    rep.singularArcNodes, s_sarc);

% --- Regular switching: Sdot -----------------------------------------------
if isnan(rep.sdotMinRel)
    s_sdot = '--';
else
    s_sdot = make_status_str(rep.sdotMinRel, sdotMin, '>');
end
fprintf(' Regular switching min|Sdot| rel    : %10.3e   (%d switches) %s\n', ...
    rep.sdotMinRel, rep.nSwitches, s_sdot);

% --- notes ------------------------------------------------------------------
if isfield(rep,'bangBangChecksRun') && ~rep.bangBangChecksRun && ~isempty(rep.Sd)
    fprintf([' NOTE: throttle cost is %s, not affine -- sign law, singular-arc and\n' ...
             '       regular-switching checks SKIPPED (bang-bang law does not apply).\n'], ...
        rep.throttleCostKind);
end

% --- Horizon note ----------------------------------------------------------
fprintf(' Horizon: %s\n', rep.horizonNote);

% --- Optional IPOPT 2nd-order check ----------------------------------------
if isfield(rep, 'ipopt') && ~isempty(rep.ipopt)
    if isfield(rep.ipopt, 'verdict') && ~isempty(rep.ipopt.verdict)
        fprintf(' 2nd-order (IPOPT inertia, delta_w) : %s\n', rep.ipopt.verdict);
    end
end

% --- ADVISORY verdict line -------------------------------------------------
advisory_status = 'PASS';
if ~rep.pass
    advisory_status = 'FAIL';
end
fprintf(' ADVISORY verdict: %s   (report-only burn-in: does not alter certified status)\n', ...
    advisory_status);

% --- Save sidecar if resDir nonempty ----------------------------------------
if ~isempty(resDir)
    sidecar_file = fullfile(resDir, sprintf('foc_%s.mat', tag));
    save(sidecar_file, 'rep');
end

end

%--------------------------------------------------------------------------
function s = make_status_str(val, threshold, cmp_op)
% Determine status string from value and comparison.
% cmp_op: '<=', '>=', '>', '<'

if isnan(val)
    s = '--';
    return
end

switch cmp_op
    case '<='
        if val <= threshold
            s = 'PASS';
        else
            s = 'FAIL';
        end
    case '>='
        if val >= threshold
            s = 'PASS';
        else
            s = 'FAIL';
        end
    case '>'
        if val > threshold
            s = 'PASS';
        else
            s = 'FAIL';
        end
    case '<'
        if val < threshold
            s = 'PASS';
        else
            s = 'FAIL';
        end
    otherwise
        error('make_status_str: unknown comparison operator');
end

end

%--------------------------------------------------------------------------
function s = make_count_status_str(count, threshold)
% Status from integer count (e.g., singularArcNodes).
% Returns PASS if count == threshold, FAIL otherwise.

if count == threshold
    s = 'PASS';
else
    s = 'FAIL';
end

end
