function H = run_highgamma_race(cfg)
% RUN_HIGHGAMMA_RACE  The "new line of attack" test proper: on every record
% the gamma-walk reached in the high-gamma band (where the direct solve
% fails, so eps never had a seed), run THREE continuation arms from the
% min-energy solution and compare -- huberc at fixed delta (0.03, then a
% delta-walk to 0.003), eps, and huber (lambda/2 seed). Each deepest
% solution is re-tested with the fixed-tf conjugate instrument.
%
% INPUTS:
%   cfg - (optional) struct: .logFile [''], .wallSec [300], .recMat
%         [direct/results/minenergy_highgamma.mat]
%
% OUTPUTS:
%   H - struct array, one row per (record, arm): .cellIdx .gamma .family
%       .pDeepest .delta .mf .mfEnergy .fails .bisects .wallMin .coast
%       .acceptOk .conjVerdict .nRungs; saved to direct/results/highgamma_race.mat
%
% REFERENCES:
%   [1] run_gamma_walk.m (the records), run_minfuel_race.m (the arms).
%   [2] FINDINGS 25-26 (the fixed-delta huberc recipe being tested here).

if nargin < 1, cfg = struct(); end
logFile = '';  wallSec = 300;  recMat = '';
if isfield(cfg, 'logFile'), logFile = cfg.logFile; end
if isfield(cfg, 'wallSec'), wallSec = cfg.wallSec; end
if isfield(cfg, 'recMat'),  recMat  = cfg.recMat;  end
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
if isempty(recMat), recMat = fullfile(resDir, 'minenergy_highgamma.mat'); end
Lr = load(recMat);  R = Lr.R;
outMat = fullfile(resDir, 'highgamma_race.mat');

arms = { 'huberc', struct('delta', 0.03, 'deltaSched', [0.02 0.012 0.008 0.005 0.003], 'seedScale', 0.5); ...
         'eps',    struct('seedScale', 1); ...
         'huber',  struct('seedScale', 0.5) };

H = struct([]);
for k = 1:numel(R)
    rec = R(k);
    for ka = 1:size(arms, 1)
        fam = arms{ka,1};  extra = arms{ka,2};
        tag = sprintf('c%d%d_g%03d_%s', rec.iD, rec.iA, round(100*rec.gam), fam);
        c2 = extra;
        c2.cell = [rec.iD rec.iA];  c2.gamma = rec.gam;  c2.families = {fam};
        c2.pilotMat = recMat;  c2.wallSec = wallSec;  c2.logFile = logFile;
        c2.outMat = fullfile(resDir, ['minfuel_hg_' tag '.mat']);
        lg('HIGH-GAMMA (%d,%d) gamma %.2f arm %s', rec.iD, rec.iA, rec.gam, fam);
        try
            if isfile(c2.outMat)
                Lo = load(c2.outMat);  out = Lo.out;
                if ~(isfield(out.arms, fam) && isfield(out.arms.(fam), 'wallTotal')), out = run_minfuel_race(c2); end
            else
                out = run_minfuel_race(c2);
            end
            A = out.arms.(fam);
        catch ME
            lg('HIGH-GAMMA ERROR (%d,%d) %.2f %s: %s', rec.iD, rec.iA, rec.gam, fam, ME.message);  continue
        end
        cv = 'n/a';
        if ~isempty(A.p)
            try
                sm = struct('family', fam, 'p', A.p(end));
                if ~isnan(A.delta(end)), sm.delta = A.delta(end); end
                Yj = A.Y{end};  Kj = size(Yj, 2);
                yT = cr3bp_minfuel_prop(rec.tf/Kj, Yj(:,end), false, rec.Tmax, rec.c, rec.muStar, sm);
                seed = struct('tf', rec.tf, 'tGrid', linspace(0, rec.tf, Kj+1), 'Y', [Yj, yT]);
                [~, it] = ms_minfuel(rec.rv0(:), rec.rvf(:), rec.tf, seed, rec.Tmax, rec.c, rec.muStar, sm, ...
                                     struct('conjTest', true, 'wallSec', wallSec));
                if it.converged, cv = it.conj.verdict; else, cv = 'UNCONVERGED'; end
            catch ME
                cv = ['ERR: ' ME.message];
            end
        end
        h = struct('cellIdx', [rec.iD rec.iA], 'gamma', rec.gam, 'family', fam, ...
            'pDeepest', lastOr(A.p, NaN), 'delta', lastOr(A.delta, NaN), ...
            'mf', lastOr(A.mf, NaN), 'mfEnergy', rec.direct.mf, ...
            'fails', A.nFail, 'bisects', A.nBisect, 'wallMin', A.wallTotal/60, ...
            'coast', lastOr(A.coastFrac, NaN), 'acceptOk', lastOr(getf(A, 'acceptOk'), false), ...
            'conjVerdict', cv, 'nRungs', numel(A.p));
        if isempty(H), H = h; else, H(end+1) = h; end %#ok<AGROW>
        save(outMat, 'H');
        lg('HIGH-GAMMA (%d,%d) %.2f %-6s: p=%.4g mf=%.6f (energy %.6f) fails=%d wall=%.1f min conj=%s', ...
           rec.iD, rec.iA, rec.gam, fam, h.pDeepest, h.mf, h.mfEnergy, h.fails, h.wallMin, cv);
    end
end
lg('HIGH-GAMMA RACE DONE -> %s', outMat);
end

% ------------------------------------------------------------------------
function v = lastOr(x, d)
% LASTOR  x(end) or d if empty.  INPUTS: x; d.  OUTPUTS: v.
if isempty(x), v = d; else, v = x(end); end
end

function v = getf(s, f)
% GETF  s.(f) or [] if absent.  INPUTS: s; f.  OUTPUTS: v.
if ~isfield(s, f), v = []; else, v = s.(f); end
end

function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
