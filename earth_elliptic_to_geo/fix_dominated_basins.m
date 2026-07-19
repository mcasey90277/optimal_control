function fix_dominated_basins()
% FIX_DOMINATED_BASINS  Recover two dominated basins exposed by the M3 front
% densify pass (Task 14): c_tf=1.35 and c_tf=2.00 each land BELOW their
% c_tf=1.20/1.75 left neighbor, violating monotonicity (more time can never
% cost fuel).
%
% c_tf=2.00 (1377.10 kg, from the ORIGINAL 5-point chain, seeded from
% c_tf=1.50) is 2.76 kg below c_tf=1.75 (1379.86 kg, a densify-pass point).
% It predates the densify targets and was never given a fresh/alternate-
% neighbor attempt. Two attempts, seeded from BOTH now-better neighbors:
%   alt1 <- sweep_T100_c175.mat
%   alt2 <- sweep_T100_c225.mat
%
% c_tf=1.35 (1360.21 kg, seeded DOWNWARD from c150, the harder direction)
% is 0.16 kg below c_tf=1.20 (1360.37 kg) -- small, but a definite
% violation, and it sits ~16 kg below c150 on a steep stretch (smells like
% a bad basin, not a flat curve). One attempt seeded UPWARD (the direction
% never tried):
%   alt1 <- sweep_T100_c120_fresh.mat
% If alt1 still doesn't beat 1360.37, also try a FRESH full-schedule (no
% seedMat, N=600) cold solve: sweep_T100_c135_fresh.
%
% Resume-safe (skip-if-exists), per-attempt try/catch, appends to
% results/densify_run.log (same log densify_front.m uses).
%
% INPUTS:  none
% OUTPUTS: none (results/sweep_T100_c200_alt{1,2}.mat,
%   results/sweep_T100_c135_alt1.mat, optionally
%   results/sweep_T100_c135_fresh.mat, + results/densify_run.log)
%
% REFERENCES: [1] densify_front.m (pattern this follows). [2] run_transfer.m
%   (certified-only save gate, no-resample seedMat warm start).
resDir = fullfile(module_root(), 'results');
logFn  = fullfile(resDir, 'densify_run.log');

attempt(resDir, logFn, 2.00, 'sweep_T100_c200_alt1', 'sweep_T100_c175.mat');
attempt(resDir, logFn, 2.00, 'sweep_T100_c200_alt2', 'sweep_T100_c225.mat');
attempt(resDir, logFn, 1.35, 'sweep_T100_c135_alt1', 'sweep_T100_c120_fresh.mat');

% conditional fresh attempt for c135: only if alt1 didn't beat c120's 1360.37
alt1Fn    = fullfile(resDir, 'sweep_T100_c135_alt1.mat');
needFresh = true;
if isfile(alt1Fn)
    S = load(alt1Fn);
    if S.res.report.m_f_kg > 1360.37
        needFresh = false;
        logmsg(logFn, sprintf(['sweep_T100_c135_alt1 (%.2f kg) beats 1360.37 -- ' ...
            'skipping fresh attempt'], S.res.report.m_f_kg));
    end
end
if needFresh
    tag = 'sweep_T100_c135_fresh';
    fn  = fullfile(resDir, [tag '.mat']);
    if isfile(fn)
        logmsg(logFn, sprintf('skip %s (already exists)', tag));
    else
        try
            res = run_transfer(struct('thrustN',10, 'ctf',1.35, 'hx0',0.0612, ...
                        'term','manifold', 'tag',tag, 'N',600));
            if res.report.certified
                logmsg(logFn, sprintf('DONE %s: certified=1 mf=%.2f kg sw=%d (fresh, no seed)', ...
                    tag, res.report.m_f_kg, res.report.switches));
            else
                logmsg(logFn, sprintf('UNCERTIFIED %s: mf=%.2f kg defect=%.2e (fresh, no seed)', ...
                    tag, res.report.m_f_kg, res.report.defect));
            end
        catch ME
            logmsg(logFn, sprintf('ERROR %s: %s', tag, ME.message));
        end
    end
end

logmsg(logFn, 'fix_dominated_basins: pass complete');
end

% ---------------------------------------------------------------------------
function attempt(resDir, logFn, cf, tag, seed)
% ATTEMPT  One resume-safe, try/catch-guarded seeded solve.
fn = fullfile(resDir, [tag '.mat']);
if isfile(fn)
    logmsg(logFn, sprintf('skip %s (already exists)', tag));
    return;
end
try
    res = run_transfer(struct('thrustN',10, 'ctf',cf, 'hx0',0.0612, ...
                'term','manifold', 'tag',tag, 'seedMat',fullfile(resDir, seed)));
    if res.report.certified
        logmsg(logFn, sprintf('DONE %s: certified=1 mf=%.2f kg sw=%d (seed=%s)', ...
            tag, res.report.m_f_kg, res.report.switches, seed));
    else
        logmsg(logFn, sprintf('UNCERTIFIED %s: mf=%.2f kg defect=%.2e (seed=%s)', ...
            tag, res.report.m_f_kg, res.report.defect, seed));
    end
catch ME
    logmsg(logFn, sprintf('ERROR %s: %s', tag, ME.message));
end
end

% ---------------------------------------------------------------------------
function logmsg(logFn, msg)
% LOGMSG  Timestamp + append a line to the densify run log, echo to stdout.
fid = fopen(logFn, 'a');
fprintf(fid, '%s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg);
fclose(fid);
fprintf('%s\n', msg);
end
