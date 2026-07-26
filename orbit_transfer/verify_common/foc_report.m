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
%           .lamMassEndMapped free-mass transversality residual (relative), Hager-
%                             mapped terminal covector -- THE GATED VALUE [scalar]
%           .lamMassEndRel    same, endpoint-extrapolated interval dual (I5),
%                             printed as an attribution companion [scalar]
%           .lamMassEndRelOneSided same, raw one-sided interval dual (legacy),
%                             printed as an attribution companion [scalar]
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

% --- KKT Stationarity ------------------------------------------------------
s_kkt = make_status_str(rep.kktStatInf, tolStat, '<=');
fprintf(' KKT stationarity  ||grad_x L||_inf : %10.3e   (sign s=%+d)   %s\n', ...
    rep.kktStatInf, rep.sLag, s_kkt);

% --- Direction: Tangential min condition -----------------------------------
s_dir = make_status_str(rep.dirTanMax, tolStat, '<=');
fprintf(' Min condition (direction) tan max  : %10.3e                 %s\n', ...
    rep.dirTanMax, s_dir);

% --- Throttle sign-law (sign-pct) ------------------------------------------
if isnan(rep.signPct)
    s_sign = '--';
else
    s_sign = make_status_str(rep.signPct, tolSign, '>=');
end
fprintf(' Min condition (throttle sign law)  : %9.2f %%                %s\n', ...
    rep.signPct, s_sign);

% --- Transversality: free-mass at terminal time (mapped terminal covector) -
% GATED VALUE is the Hager-mapped terminal covector (finding I5's principled
% fix): its mass component is exact by the discrete KKT system when the
% final mass is genuinely free, unlike the raw/extrapolated interval duals
% below, which carry an O(h) endpoint-representation offset by construction.
if isnan(rep.lamMassEndMapped)
    s_trans = '--';
else
    s_trans = make_status_str(rep.lamMassEndMapped, tolTrans, '<=');
end
fprintf(' Transversality |lam_m(tf)| mapped   : %10.3e                 %s\n', ...
    rep.lamMassEndMapped, s_trans);

% --- Time-costate behavior -------------------------------------------------
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

% --- attribution lines (external-review findings I1/I5, 2026-07-25) ---------
% Print the two superseded statistics beside the gated one so a marginal miss
% can be ATTRIBUTED rather than merely moved: a large gap between mapped and
% raw/extrapolated means the old number was an endpoint-representation
% artifact, a small gap means the finding is real (all three converge).
if isfield(rep,'lamMassEndRel') && ~isnan(rep.lamMassEndRel)
    fprintf('   ^ transversality, extrapolated (I5): %10.3e\n', rep.lamMassEndRel);
end
if isfield(rep,'lamMassEndRelOneSided') && ~isnan(rep.lamMassEndRelOneSided)
    fprintf('   ^ transversality, legacy one-sided : %10.3e\n', ...
        rep.lamMassEndRelOneSided);
end
if isfield(rep,'sdotMinRelLegacy') && ~isnan(rep.sdotMinRelLegacy) && ~isnan(rep.sdotMinRel)
    fprintf('   ^ regular switching, legacy raw    : %10.3e   [I1 mesh fix]\n', ...
        rep.sdotMinRelLegacy);
end
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
