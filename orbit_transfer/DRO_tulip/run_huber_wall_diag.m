function D = run_huber_wall_diag(cfg)
% RUN_HUBER_WALL_DIAG  MfMax-review steps 1 + 2 on the two Huber walls of
% FINDINGS 23: (6,8)@1.2 (walled at p ~ 0.033) and (1,2)@1.2 with the
% lambda/2 seed (walled at p ~ 0.77).
%
% Two configurations per wall, run sequentially (one parpool):
%   TIGHT  -- the 09-02/09-05 protocol (every rung converged to ms tolR), now
%             logging cond(J) at every accepted rung AND every failed iterate.
%             Fold diagnosis: cond(J) blowing up on approach = a turning
%             point in p (arclength needed); flat = the field / grazing.
%   LOOSE  -- MfMax's homCI acceptance: a rung is accepted at normR < 1e-3
%             (Hdrift < 1e-2), the FLOOR is refined to full tolerance and a
%             loose floor that will not tighten falls back to the deepest
%             tight rung. Tests whether the walls are a gate artifact.
%
% INPUTS:
%   cfg - (optional) struct: .logFile [''], .wallSec [300]
%
% OUTPUTS:
%   D - struct array, one row per (wall, config): .cell .gamma .config
%       .pDeepest .mf .nFail .nBisect .wallMin .condJ (accepted rungs)
%       .failP .failCondJ .failNormR .nLoose
%
% REFERENCES:
%   [1] doc/mfmax_ideas_review.md sections 1.1, 1.2, 3 (steps 1-2).
%   [2] DRO_tulip/FINDINGS.md section 23 (the walls being diagnosed).

if nargin < 1, cfg = struct(); end
logFile = '';  wallSec = 300;
if isfield(cfg, 'logFile'), logFile = cfg.logFile; end
if isfield(cfg, 'wallSec'), wallSec = cfg.wallSec; end
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
outMat = fullfile(resDir, 'huber_wall_diag.mat');

walls = {[6 8], 1.2, 1.0; [1 2], 1.2, 0.5};
cfgs  = {'tight', struct(); 'loose', struct('rungTolR', 1e-3, 'rungHdriftTol', 1e-2)};

D = struct([]);
for kw = 1:size(walls, 1)
    for kc = 1:size(cfgs, 1)
        cell_ = walls{kw,1};  gam = walls{kw,2};  ss = walls{kw,3};
        tag = sprintf('c%d%d_g%03d_%s', cell_, round(100*gam), cfgs{kc,1});
        c2 = cfgs{kc,2};
        c2.cell = cell_;  c2.gamma = gam;  c2.families = {'huber'};
        c2.seedScale = ss;  c2.wallSec = wallSec;  c2.logFile = logFile;
        c2.outMat = fullfile(resDir, ['minfuel_huber_' tag '.mat']);
        lg('WALL DIAG: (%d,%d) gamma %.2f seed %.2g config %s', cell_, gam, ss, cfgs{kc,1});
        try
            out = run_minfuel_race(c2);  A = out.arms.huber;
        catch ME
            lg('WALL DIAG ERROR: %s', ME.message);  continue
        end
        d = struct('cell', cell_, 'gamma', gam, 'config', cfgs{kc,1}, ...
            'pDeepest', lastOr(A.p, NaN), 'mf', lastOr(A.mf, NaN), ...
            'nFail', A.nFail, 'nBisect', A.nBisect, 'wallMin', A.wallTotal/60, ...
            'p', A.p, 'condJ', A.condJ, 'nLoose', nnz(~A.tight), ...
            'failP', A.failP, 'failCondJ', A.failCondJ, 'failNormR', A.failNormR);
        if isempty(D), D = d; else, D(end+1) = d; end %#ok<AGROW>
        save(outMat, 'D');
        lg('WALL DIAG (%d,%d) %s: deepest p=%.4g mf=%.6f fails=%d loose-accepted=%d', ...
           cell_, cfgs{kc,1}, d.pDeepest, d.mf, d.nFail, d.nLoose);
        lg('   condJ along accepted rungs: %s', mat2str(A.condJ, 3));
        lg('   condJ at failed iterates  : %s', mat2str(A.failCondJ, 3));
    end
end
lg('WALL DIAG DONE -> %s', outMat);
end

% ------------------------------------------------------------------------
function v = lastOr(x, d)
% LASTOR  x(end) or d if empty.  INPUTS: x; d.  OUTPUTS: v.
if isempty(x), v = d; else, v = x(end); end
end

function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
