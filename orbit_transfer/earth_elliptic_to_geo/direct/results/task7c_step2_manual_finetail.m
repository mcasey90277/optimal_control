% TASK7C_STEP2_MANUAL_FINETAIL  Manual fine eps-walk toward 0 for the 1 N
% fuel solve (npr=25, N=1740), resuming from the warm-tail's step11 (eps=
% 0.001, defect=5.1154e-14 -- essentially exact) rather than jumping
% straight to eps=0 as homotopy_mee.m's fixed schedule does.
%
% ROOT-CAUSE NOTE: homotopy_mee.m's warm-start chain (Xk/Uk/dLk) is ONLY
% advanced on a step that clears defect<1e-8 ("a loose iterate must not
% poison the chain" -- a deliberate, documented design choice, not a bug).
% This is CORRECT in spirit but means a run of several non-certifying steps
% collapses into one large eps jump from the last good point. Here, the
% eps=0.2->eps=0.001 jump (skipping 8 schedule entries) happened to
% succeed, but the immediate NEXT jump (eps=0.001 -> eps=0, warmTight=true,
% a SMALL jump from an essentially-exact point) landed in a completely
% different, uncertified basin (defect=6.31e-02) -- bit-identical to the
% totally cold direct-eps=0 attempt's own result, both at npr=25 and
% (separately) at npr=15, suggesting eps=0's true bang-bang landscape has a
% strongly attracting bad local optimum near this warm-start neighborhood,
% not merely a "too big a jump" problem. This script tries a genuinely
% GRADUAL final approach (many small eps decrements, ALWAYS advancing the
% warm start regardless of certification -- deliberately reversing the
% freeze-on-failure guard for just this manual tail, since the risk here is
% under-exploring near eps=0, not chain-poisoning from a wildly bad step).
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here); addpath(here);
resDir = fullfile(here, 'results');

S11 = load(fullfile(resDir, 'MEE_M2_1N_warmtail_step11.mat'));
o = S11.o;
fprintf('RESUME from warmtail step11: eps=%.4g defect=%.4e sw=%d N=%d\n', ...
        S11.e, o.maxDefect, o.switches, size(o.X,2)-1);

par = kepler_lt_params(1, 1500, 2000);
tfTarget = 1.5 * 223.808136;   % ctf * step1b anchor tfmin
x0state = o.X(:,1);
sigma = linspace(0,1,size(o.X,2)).';

schedFine = [7e-4 5e-4 3e-4 2e-4 1e-4 5e-5 2e-5 1e-5 5e-6 2e-6 1e-6 0];

isGood = @(oo) oo.success && oo.maxDefect < 1e-8 && oo.termErr < 1e-8;
rf0 = fullfile(resDir, 'task7c_manual_finetail_state.mat');

results = struct('e', {}, 'defect', {}, 'sw', {}, 'status', {});
for k = 1:numel(schedFine)
    e = schedFine(k);
    rf = fullfile(resDir, sprintf('task7c_finetail_round%02d.mat', k));
    if isfile(rf)
        S = load(rf); onew = S.onew;
        fprintf('  [cached] fine round %d\n', k);
    else
        t0 = tic;
        onew = casadi_lt_mee(sigma, o.X, o.U, o.dL, struct('par', par, 'mode', 'fixedtf', ...
            'eps', e, 'tfTarget', tfTarget, 'x0', x0state, 'maxIter', 800, ...
            'warmTight', true, 'printLevel', 3));
        rt = toc(t0);
        save(rf, 'onew', 'rt', 'e');
    end
    fprintf('FINE round %d (eps=%.4g): defect %.4e -> %.4e, sw=%d->%d, status=%s\n', ...
            k, e, o.maxDefect, onew.maxDefect, o.switches, onew.switches, onew.ipoptStatus);
    o = onew;   % ALWAYS advance (see ROOT-CAUSE NOTE above)
    results(end+1) = struct('e', e, 'defect', onew.maxDefect, 'sw', onew.switches, 'status', onew.ipoptStatus); %#ok<SAGROW>
    save(rf0, 'o', 'results');
    if isGood(o) && e == 0
        fprintf('CERTIFIED at eps=0!\n');
        break;
    end
end

fprintf('\nFINAL: eps=%.4g defect=%.4e termErr=%.4e status=%s certified(e==0 & good)=%d\n', ...
        schedFine(min(k,numel(schedFine))), o.maxDefect, o.termErr, o.ipoptStatus, ...
        isGood(o) && schedFine(min(k,numel(schedFine)))==0);
save(fullfile(resDir, 'task7c_step2_finetail_final.mat'), 'o', 'results');
