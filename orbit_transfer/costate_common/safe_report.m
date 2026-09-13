function ok = safe_report(fcn, label)
%% Purpose:
%
%   Run a REPORTING block so that it can never fail the work it reports on.
%
%   Why this exists (2026-09-13). Twice in one campaign a completed,
%   correctly saved stage was reported as a FAILED job because its summary
%   print threw: once on char(string(NaN)) for a phase with no certified
%   transfer, once on [c.ok] where a column with no candidates holds a plain
%   double rather than an empty struct. In both cases the science was on
%   disk and the job said it had failed.
%
%   Printing is not the work. A report that throws is a bug in the report,
%   and it must be visible as exactly that: this catches, says so loudly,
%   and returns false without disturbing the caller.
%
%   Use fmtNum/fmtCount beside it for the two shapes that actually broke.
%
%% Inputs:
%
%  fcn                      fhandle                 the reporting block
%  label                    char                    what it was reporting
%
%% Outputs:
%
%  ok                       logical                 the block ran cleanly
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
try
    fcn();
catch ME
    ok = false;
    fprintf(2, ['\n*** REPORT FAILED (%s): %s\n' ...
                '    This is a bug in the REPORT, not in the work: the stage''s\n' ...
                '    output is on disk and unaffected. %s\n\n'], ...
            label, ME.message, ME.identifier);
end
end
