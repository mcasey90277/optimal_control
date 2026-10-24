function tbl = drone_table1(opts)
%% Purpose:
%
%   Reproduce Table 1 of Zhang, Liu, Darbon & Karniadakis (arXiv:
%   2610.30388) with certified classical solves: for each M, the hard-
%   clearance optimum (the reference) and the paper's barrier form, each
%   through verify_drones_direct, written to results/table1.csv beside the
%   paper's OCFormer and SympOCNet costs, with the M = 4 and M = 8 fixtures
%   saved for Route 2. The paper's own numbers (its Table 1):
%
%       M     J_LB    OCFormer  SympOCNet
%       4     0.200   0.2003    0.2001
%       8     0.400   0.4036    0.4046
%       16    0.800   0.8307    0.8345
%       32    1.156   1.2531    1.4413
%       64    2.297   2.448     --
%       128   4.704   5.1424    13.3 (did not converge)
%
%   The paper's numbers are REPORTED, never gated: the rigorous antipodal-
%   pair bound (column pairBound) exceeds the paper's SympOCNet 0.2001 at
%   M = 4. M = 64 and 128 are STRETCH: run them from batch/run_table1.sh
%   with a log file, never from the test suite.
%
%   Hard rows escalate the mesh N = 100 -> 200 -> 400 while the mesh gate
%   (G4) fails. The barrier row warm-starts from the hard solution, so it
%   needs a hard row ('hard' must be in opts.modes) and solves ONCE at the
%   hard row's final N; its gates are recorded as they come out (gatesOk
%   false if its G4 fails). A skipped mirror gate records mirrorOk false.
%   When the mirror swirl is cheaper it is adopted and re-verified without
%   a mirror, and the original G5 block is carried over with .out = the
%   previous primary (the adopted solution's mirror).
%
%% Inputs:
%
%  opts                     struct (optional)       .Ms [4 8 16 32]
%                                                   .modes {'hard',
%                                                   'barrier'} .mesh
%                                                   [true for M <= 16]
%                                                   .saveFixtures [true]
%                                                   .csv [results/
%                                                   table1.csv] .logFile
%                                                   [''] .N [100]; no
%                                                   argument runs the demo
%                                                   (M = 4, both modes, csv
%                                                   to a tempname, no
%                                                   fixtures)
%
%% Outputs:
%
%  tbl                      table                   One row per (M, mode)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: M = 4, both modes, mesh gate on; never touches results/ or data/:
     tbl = drone_table1(struct('Ms', 4, 'modes', {{'hard', 'barrier'}}, 'mesh', true, ...
                               'saveFixtures', false, 'csv', [tempname '.csv']));
     return
end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
d = @(f, v) fieldd(opts, f, v);
Ms     = d('Ms', [4 8 16 32]);
modes  = d('modes', {'hard', 'barrier'});
doMesh = d('mesh', []);
saveFx = d('saveFixtures', true);
csvOut = d('csv', fullfile(here, 'results', 'table1.csv'));
logF   = d('logFile', '');
N      = d('N', 100);
lg = @(varargin) logmsg(logF, sprintf(varargin{:}));
if ~iscellstr(modes) || ~all(ismember(modes, {'hard', 'barrier'})) || ~any(strcmp(modes, 'hard'))
    error('drone_table1:modes', ...
          'opts.modes must include ''hard'' (barrier rows warm-start from the hard solution) and hold only ''hard''/''barrier''');
end
modes = [{'hard'}, setdiff(modes, {'hard'}, 'stable')];   % the hard row first: the barrier row needs it

paperMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
paperMap(4)   = [0.2003, 0.2001];   paperMap(8)   = [0.4036, 0.4046];
paperMap(16)  = [0.8307, 0.8345];   paperMap(32)  = [1.2531, 1.4413];
paperMap(64)  = [2.448,  NaN];      paperMap(128) = [5.1424, 13.3];

rows = {};
for M = Ms
    [x0, xT, geom] = drone_layouts('paper', M);
    mesh = doMesh;  if isempty(mesh), mesh = M <= 16; end
    store = struct();
    for mm = 1:numel(modes)
        mode = modes{mm};
        if strcmp(mode, 'hard')
            % hard rows escalate the mesh while G4 fails: N -> 2N -> 4N
            Ncur = N;
            while true
                lg('M = %d, mode hard: solving N = %d', M, Ncur);
                o = solve_drones_direct(x0, xT, geom, struct('mode', 'hard', 'N', Ncur));
                rep = verify_drones_direct(o, struct('mirror', true, 'mesh', mesh));
                if rep.G4.ran && ~rep.G4.ok && Ncur < 4*N
                    Ncur = 2*Ncur;
                    lg('M = %d, mode hard: mesh gate failed (%.2e), escalating to N = %d', M, rep.G4.relChange, Ncur);
                else
                    break
                end
            end
        else
            % the barrier row: ONE solve at the hard row's final N, warm-started
            % from the hard solution; its gates are recorded as they come out
            Nb = store.out.N;
            lg('M = %d, mode barrier: solving N = %d (the hard row''s final N), warm-started from the hard solution', M, Nb);
            o = solve_drones_direct(x0, xT, geom, struct('mode', 'barrier', 'N', Nb, 'X0', store.out.X, 'V0', store.out.V));
            vopts = struct('mirror', true, 'mesh', mesh);
            if store.rep.G5.ran && ~isempty(store.rep.G5.out)
                vopts.mirrorX0 = store.rep.G5.out.X;  vopts.mirrorV0 = store.rep.G5.out.V;
                lg('M = %d, mode barrier: barrier mirror warm-started from the hard mirror solution', M);
            end
            rep = verify_drones_direct(o, vopts);
        end
        mirrorOk = rep.G5.ok;  mirrorRel = rep.G5.relDiff;
        if rep.G5.ran && ~isempty(rep.G5.out) && rep.G5.out.J < o.J - 1e-12
            % adopt the cheaper mirror swirl; re-verify it without a mirror, then
            % carry the original G5 block over with .out = the previous primary
            % (the adopted solution's mirror), so the barrier warm start and
            % the stored report still see a mirror that ran
            G5 = rep.G5;  prev = o;
            o = rep.G5.out;
            rep = verify_drones_direct(o, struct('mirror', false, 'mesh', mesh));
            G5.out = prev;  G5.ran = true;  G5.Jmirror = prev.J;
            G5.status = prev.ipoptStatus;  G5.success = prev.success;
            rep.G5 = G5;
            lg('M = %d, mode %s: adopted the mirror swirl (lower cost)', M, mode);
        end
        gatesOk = rep.G1.ok && rep.G2.ok && rep.G3.ok && rep.G4.ok;
        pp = paperMap(M);
        pb = antipodal_pair_bound(x0, xT, o.T, geom.r);
        rows(end+1, :) = {M, mode, o.N, o.JLB, pb, o.J, 100*(o.J - o.JLB)/o.JLB, pp(1), pp(2), ...
                          rep.G1.minHfine, o.maxDefect, rep.G4.relChange, mirrorRel, mirrorOk, ...
                          o.iterations, o.wall, o.ipoptStatus, gatesOk};
        lg('M = %d, mode %s: N = %d, J = %.6f (J_LB %.4f, pairBound %.6f, OCFormer %.4f), %s, gates %s, mirror %s (%.2e), wall %.1f s', ...
           M, mode, o.N, o.J, o.JLB, pb, pp(1), o.ipoptStatus, mat2str(gatesOk), mat2str(mirrorOk), mirrorRel, o.wall);
        if strcmp(mode, 'hard'), store.out = o;  store.rep = rep;
        else,                    store.outB = o; store.repB = rep;
        end
    end
    if saveFx && any(M == [4 8]) && isfield(store, 'out') && isfield(store, 'outB')
        out = store.out;  outB = store.outB;  rep = store.rep;  repB = store.repB;
        pp = paperMap(M);
        paper = struct('OCFormer', pp(1), 'SympOCNet', pp(2));
        if ~isfolder(fullfile(here, 'data')), mkdir(fullfile(here, 'data')); end
        save(fullfile(here, 'data', sprintf('drones_direct_M%d.mat', M)), ...
             'out', 'outB', 'rep', 'repB', 'geom', 'x0', 'xT', 'paper');
        lg('M = %d: fixture saved', M);
    end
end
tbl = cell2table(rows, 'VariableNames', {'M', 'mode', 'N', 'J_LB', 'pairBound', 'J', 'gapPct', 'OCFormer', 'SympOCNet', ...
                 'minHfine', 'maxDefect', 'meshRel', 'mirrorRel', 'mirrorOk', 'iters', 'wallSec', 'status', 'gatesOk'});
if ~isfolder(fileparts(csvOut)), mkdir(fileparts(csvOut)); end
writetable(tbl, csvOut);
lg('wrote %s', csvOut);
disp(tbl);
end

% ------------------------------------------------------------------------
function logmsg(f, s)
%% Purpose:
%
%   Append one line to the log file, or print it when no file is set.
%
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a');  fprintf(fid, '%s\n', s);  fclose(fid);
end
end

function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
