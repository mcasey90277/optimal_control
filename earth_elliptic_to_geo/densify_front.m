function densify_front()
% DENSIFY_FRONT  Add 4 densifying points to the c_tf front (T_max=10 N).
%
% Adds c_tf in {1.35, 1.75, 2.25, 2.75} to the existing {1.2, 1.5, 2.0, 2.5,
% 3.0} front (run_ctf_sweep.m leg 1), each seeded from its nearest EXISTING
% banked neighbor (no-resample rule, run_transfer.m):
%   1.35 <- sweep_T100_c150.mat  (rescale DOWN, the harder direction; on
%           certification failure, retry seeded from sweep_T100_c120_fresh.mat
%           going UP instead)
%   1.75 <- sweep_T100_c150.mat
%   2.25 <- sweep_T100_c200.mat
%   2.75 <- sweep_T100_c250_fresh.mat
%
% Resume-safe: skips any point whose results/<tag>.mat already exists
% (run_transfer only saves CERTIFIED results, campaign rule -- so an
% existing file always means "done and good"). Per-point try/catch so one
% point's solver error does not abort the other three. Progress and outcome
% for every point (skip / done / uncertified-retry / failed / error) is
% appended to results/densify_run.log, one line per event, so a relaunch
% after the documented CasADi MEX init crash (dies at the IPOPT banner, ~1
% in 10 solves, see run_ctf_sweep.m header) picks up exactly where it left
% off just by re-running this function.
%
% INPUTS:  none
% OUTPUTS: none (results/sweep_T100_c{135,175,225,275}.mat + densify_run.log)
%
% REFERENCES: [1] run_ctf_sweep.m (front sweep, seedMat no-resample warm
%   start, certified-only save gate). [2] run_transfer.m.
resDir = fullfile(module_root(), 'results');
logFn  = fullfile(resDir, 'densify_run.log');

pts = struct('cf', {1.35, 1.75, 2.25, 2.75}, ...
             'seed', {'sweep_T100_c150.mat', 'sweep_T100_c150.mat', ...
                      'sweep_T100_c200.mat', 'sweep_T100_c250_fresh.mat'}, ...
             'fallback', {'sweep_T100_c120_fresh.mat', '', '', ''});

for k = 1:numel(pts)
    pt  = pts(k);
    tag = sprintf('sweep_T100_c%03d', round(100*pt.cf));
    fn  = fullfile(resDir, [tag '.mat']);
    if isfile(fn)
        logmsg(logFn, sprintf('skip %s (already exists)', tag));
        continue;
    end
    try
        res = run_transfer(struct('thrustN',10, 'ctf',pt.cf, 'hx0',0.0612, ...
                    'term','manifold', 'tag',tag, ...
                    'seedMat',fullfile(resDir, pt.seed)));
        if res.report.certified
            logmsg(logFn, sprintf('DONE %s: certified=1 mf=%.2f kg sw=%d (seed=%s)', ...
                tag, res.report.m_f_kg, res.report.switches, pt.seed));
        else
            logmsg(logFn, sprintf('UNCERTIFIED %s: mf=%.2f kg defect=%.2e (seed=%s)', ...
                tag, res.report.m_f_kg, res.report.defect, pt.seed));
            if ~isempty(pt.fallback)
                logmsg(logFn, sprintf('retry %s seeded from fallback %s', tag, pt.fallback));
                res2 = run_transfer(struct('thrustN',10, 'ctf',pt.cf, 'hx0',0.0612, ...
                            'term','manifold', 'tag',tag, ...
                            'seedMat',fullfile(resDir, pt.fallback)));
                if res2.report.certified
                    logmsg(logFn, sprintf('DONE %s: certified=1 mf=%.2f kg sw=%d (fallback seed=%s)', ...
                        tag, res2.report.m_f_kg, res2.report.switches, pt.fallback));
                else
                    logmsg(logFn, sprintf('FAILED %s: both primary and fallback seeds uncertified', tag));
                end
            else
                logmsg(logFn, sprintf('FAILED %s: no fallback seed configured for this point', tag));
            end
        end
    catch ME
        logmsg(logFn, sprintf('ERROR %s: %s', tag, ME.message));
    end
end
logmsg(logFn, 'densify_front: pass complete');
end

% ---------------------------------------------------------------------------
function logmsg(logFn, msg)
% LOGMSG  Timestamp + append a line to the densify run log, echo to stdout.
fid = fopen(logFn, 'a');
fprintf(fid, '%s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg);
fclose(fid);
fprintf('%s\n', msg);
end
