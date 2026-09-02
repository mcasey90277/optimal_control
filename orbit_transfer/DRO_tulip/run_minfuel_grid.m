function G = run_minfuel_grid(cfg)
% RUN_MINFUEL_GRID  Step-5 Task 5: the eps continuation (race winner) over
% every PASSING min-energy backbone record -- the program's FIRST MIN-FUEL
% record set.
%
% For each (cell, gamma) record in minenergy_pilot.mat with verdict PASS,
% runs the eps arm of run_minfuel_race (same schedule, gates, hard caps;
% per-record out .mat), then aggregates: eps-deepest, m_f (fuel) vs m_f
% (energy), coast fraction, acceptance, and the FULL ms junction states of
% the deepest solution (identifiability rule). Saves
% direct/results/minfuel_grid.mat and prints the summary table.
%
% INPUTS:
%   cfg - (optional) struct: .logFile ['']; other fields forwarded to
%         run_minfuel_race (e.g. .sched, .wallSec).
%
% OUTPUTS:
%   G - struct array, one per record: .cellIdx .gamma .pDeepest .mfFuel
%       .mfEnergy .dMfPct .coastFrac .acceptOk .acceptDz .z .Y .tGrid
%       .nFail .nBisect .wallMin
%
% REFERENCES:
%   [1] run_minfuel_race.m (the walk; eps arm = the 2026-09-02 race winner).
%   [2] run_minenergy_pilot.m (the backbone records).

if nargin < 1, cfg = struct(); end
logFile = '';
if isfield(cfg, 'logFile'), logFile = cfg.logFile; cfg = rmfield(cfg, 'logFile'); end
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'direct', 'results');
L = load(fullfile(resDir, 'minenergy_pilot.mat'));
sel = find(arrayfun(@(r) (ischar(r.verdict) && strcmp(r.verdict, 'PASS')) || ...
                         (~ischar(r.verdict) && ~isempty(r.gates) && all(r.gates)), L.R));
lg('GRID: %d passing backbone records', numel(sel));

G = struct([]);
for ks = 1:numel(sel)
    rec = L.R(sel(ks));
    tag = sprintf('c%d%d_g%03d', rec.iD, rec.iA, round(100*rec.gam));
    om  = fullfile(resDir, ['minfuel_eps_' tag '.mat']);
    lg('GRID: (%d,%d) gamma %.2f -> %s', rec.iD, rec.iA, rec.gam, om);
    c2 = cfg;
    c2.cell = [rec.iD rec.iA];  c2.gamma = rec.gam;
    c2.families = {'eps'};  c2.outMat = om;  c2.logFile = logFile;
    try
        out = run_minfuel_race(c2);
        A = out.arms.eps;
        g = struct('cellIdx', c2.cell, 'gamma', rec.gam, ...
            'pDeepest', A.p(end), 'mfFuel', A.mf(end), ...
            'mfEnergy', rec.direct.mf, ...
            'dMfPct', 100*(A.mf(end) - rec.direct.mf), ...
            'coastFrac', A.coastFrac(end), 'acceptOk', A.acceptOk, ...
            'acceptDz', A.acceptDz, 'z', A.z, 'Y', A.Y{end}, ...
            'tGrid', [], 'nFail', A.nFail, 'nBisect', A.nBisect, ...
            'wallMin', A.wallTotal/60);
        if isempty(G), G = g; else, G(end+1) = g; end %#ok<AGROW>
    catch ME
        lg('GRID: (%d,%d) gamma %.2f ERROR: %s', rec.iD, rec.iA, rec.gam, ME.message);
    end
    save(fullfile(resDir, 'minfuel_grid.mat'), 'G');
end

lg('GRID SUMMARY:');
lg('cell   gam   eps_min    mf_fuel   mf_energy  dmf[%%m0]  coast  accept  fails');
for k = 1:numel(G)
    lg('(%d,%d) %.2f  %.4g   %.6f  %.6f   +%.3f    %.2f   %d      %d', ...
       G(k).cellIdx(1), G(k).cellIdx(2), G(k).gamma, G(k).pDeepest, ...
       G(k).mfFuel, G(k).mfEnergy, G(k).dMfPct, G(k).coastFrac, ...
       G(k).acceptOk, G(k).nFail);
end
lg('GRID DONE: %d/%d records', numel(G), numel(sel));
end

% ------------------------------------------------------------------------
function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
end
