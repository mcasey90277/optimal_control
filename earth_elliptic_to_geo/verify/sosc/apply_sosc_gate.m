function report = apply_sosc_gate(report, sosc)
% APPLY_SOSC_GATE  Tiered SOSC gate: only a PROVEN saddle (verdict FAIL)
% demotes a feasibility-certified row. PASS/WEAK_MIN/INCONCLUSIVE/ERROR all
% keep report.certified true -- WEAK_MIN is a POSITIVE certificate (reduced
% Hessian PSD with no descent direction, the expected outcome for bang-bang
% min-fuel), and INCONCLUSIVE/ERROR are non-demoting (the 2nd-order
% certificate simply could not be established, which is not evidence AGAINST
% the point, only an open question). See process/DESIGN_sosc.md sec 11.6.
%
% INPUTS:
%   report - run_transfer_mee (or reproducer) report struct; must have
%            .certified [logical]                                    [struct]
%   sosc    - verify_sosc_mee output struct: .verdict in {PASS,WEAK_MIN,
%             FAIL,INCONCLUSIVE,ERROR}, .status, .reason (+ other fields
%             attached as-is)                                        [struct]
%
% OUTPUTS:
%   report - with .sosc = sosc attached, and .certified demoted to false
%            IFF sosc.verdict == 'FAIL' (a proven, non-degenerate saddle)
%
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 11.6 (revised tiered gate).
%   [2] .superpowers/sdd/task-9-brief.md (original apply_sosc_gate spec).
report.sosc = sosc;
if strcmp(sosc.verdict, 'FAIL')
    report.certified = false;
    warning('apply_sosc_gate:fail', ...
        'SOSC FAIL (proven saddle) -> demoted to feasible-only: %s', sosc.reason);
end
end
