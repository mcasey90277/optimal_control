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
%   Hard rows. (1) A small MULTISTART at the starting N, four cold seeds:
%   'u+' / 'u-' (uniform bump, swirl sense +1 / -1) and 'a+' / 'a-'
%   (alternating amplitudes: bumpScale 1.2 on the odd agents of each ring,
%   0.6 on the even, ring = agents of one |x0| radius in index order). The
%   alternating seeds open the symmetry-broken basins a uniform seed cannot
%   reach (at M = 8 the Z4 "two squares" basin lies 0.12 % below the
%   regular octagon). Each seed is verified with mesh and mirror off; the
%   cheapest whose own G1..G3 pass is the row (columns seedsTried,
%   bestSeed; every seed's J and status is logged). (2) The mesh loop: G4
%   re-solves at 2N WARM from the interpolant (IPOPT warm_start_init_point,
%   mu_init 1e-4). A cheaper 2N solution is ADOPTED
%   (adopt_mesh_solution) if it passes G1..G3 itself (checked without a
%   solve before the swap, so a gated row is never replaced by an ungated
%   one): it becomes the row at 2N and is re-verified at
%   its own N -- a cost drop at 2N may be a different basin, never assumed
%   to be discretisation error. A G4 failure with no cost drop escalates
%   to 2N warm from the interpolant, as before. Both stop at N = 4*opts.N;
%   at the cap G4 is a gate only. (3) The mirror gate (G5) runs once on the
%   final row; a cheaper mirror is adopted as before.
%
%   Barrier rows are the hard basin's twin BY CONSTRUCTION: one solve at
%   the hard row's final N, warm-started from the hard solution (no
%   multistart, no escalation), its mirror warm-started from the hard
%   mirror; gates recorded as they come out. The barrier row's J INCLUDES
%   the penalty integral (trapezoid of the barrier U at the nodes); the
%   log line splits J into its smooth part and the penalty and prints the
%   minimum NODAL gap h against l with the barrier branch it lands in
%   (log for h > l, quadratic for h <= l). Measured (batch 2026-10-23):
%   the penalty is 8.64e-4 = 0.43 % of J at M = 4 (smooth part 0.2006460)
%   and 8.82e-2 = 6.80 % of J at M = 32 (where the barrier row sits 7.4 %
%   above the hard row); min nodal h 8.0e-5 / 1.18e-4 / 8.6e-5 / 6.9e-5 at
%   M = 4 / 8 / 16 / 32, so only M = 8 stays on the log branch.
%
%   gatesRun is 'G1-G4' when the mesh gate ran and 'G1-G3' when it did not
%   (M = 32); gatesOk is the AND of the gates that ran. A skipped mirror
%   gate records mirrorOk false.
%
%   adopt = drone_table1('adopt_mesh_solution') returns the handle of the
%   adoption decision (tests/test_table1_fixture.m exercises it on a
%   synthetic report).
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
%                                                   fixtures); the string
%                                                   'adopt_mesh_solution'
%                                                   returns that decision's
%                                                   function handle instead
%                                                   of a table
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
if ischar(opts) && strcmp(opts, 'adopt_mesh_solution')
    tbl = @adopt_mesh_solution;
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
            % (1) the multistart at the starting N
            seeds = hard_seeds(x0);
            best = [];  bestLab = '';  nTried = 0;
            for ks = 1:numel(seeds)
                lg('M = %d, mode hard: seed %s, solving N = %d (cold)', M, seeds(ks).label, N);
                os = solve_drones_direct(x0, xT, geom, struct('mode', 'hard', 'N', N, ...
                                         'sense', seeds(ks).sense, 'bumpScale', seeds(ks).scale));
                rs = verify_drones_direct(os, struct('mirror', false, 'mesh', false));
                nTried = nTried + 1;
                g13 = rs.G1.ok && rs.G2.ok && rs.G3.ok;
                lg('M = %d, mode hard: seed %s: J = %.7f, %s, %d iterations, G1-G3 %s, wall %.1f s', ...
                   M, seeds(ks).label, os.J, os.ipoptStatus, os.iterations, mat2str(g13), os.wall);
                if g13 && (isempty(best) || os.J < best.J)
                    best = os;  bestLab = seeds(ks).label;
                end
            end
            if isempty(best)
                error('drone_table1:seeds', 'M = %d: no hard seed passed G1-G3', M);
            end
            lg('M = %d, mode hard: best seed %s, J = %.7f', M, bestLab, best.J);
            % (2) the mesh loop: adopt a cheaper 2N solution, else escalate on a G4 failure
            o = best;  Nmax = 4*N;
            while true
                rep = verify_drones_direct(o, struct('mirror', false, 'mesh', mesh));
                if rep.G4.ran
                    lg('M = %d, mode hard: N = %d, G4 2N = %d: J = %.7f -> %.7f (rel %.2e, %s)', ...
                       M, o.N, 2*o.N, o.J, rep.G4.J2N, rep.G4.relChange, rep.G4.status);
                end
                if adopt_mesh_solution(o, rep)
                    if 2*o.N <= Nmax
                        % the cheaper 2N solution replaces a gated row only if it
                        % passes G1-G3 itself (no solve: clearance, bound, defects)
                        repC = verify_drones_direct(rep.G4.out, struct('mirror', false, 'mesh', false));
                        if repC.G1.ok && repC.G2.ok && repC.G3.ok
                            lg('M = %d, mode hard: ADOPTED the cheaper 2N solution (J %.7f < %.7f), row now at N = %d', ...
                               M, rep.G4.out.J, o.J, 2*o.N);
                            o = rep.G4.out;
                            continue
                        end
                        lg('M = %d, mode hard: the cheaper 2N solution fails its own gates (G1 %d, G2 %d, G3 %d): not adopted', ...
                           M, repC.G1.ok, repC.G2.ok, repC.G3.ok);
                    else
                        lg('M = %d, mode hard: the 2N solution is cheaper but N = %d is the cap: G4 is a gate only', M, o.N);
                    end
                end
                if rep.G4.ran && ~rep.G4.ok && o.N < Nmax
                    Nn = 2*o.N;
                    lg('M = %d, mode hard: mesh gate failed (%.2e) with no cost drop, escalating to N = %d (warm)', ...
                       M, rep.G4.relChange, Nn);
                    [Xw, Vw] = hermite_resample(o.X, o.V, o.tN, 2);
                    o = solve_drones_direct(x0, xT, geom, struct('mode', 'hard', 'N', Nn, 'sense', o.sense, ...
                                            'bumpScale', o.bumpScale, 'X0', Xw, 'V0', Vw));
                    continue
                end
                break
            end
            % (3) the mirror once, on the final row (G4 carried over, not re-solved)
            rep = verify_drones_direct(o, struct('mirror', true, 'mesh', mesh, 'G4', rep.G4));
            seedLab = bestLab;
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
            nTried = 0;  seedLab = 'hardWarm';
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
        if rep.G4.ran
            gatesRun = 'G1-G4';  gatesOk = rep.G1.ok && rep.G2.ok && rep.G3.ok && rep.G4.ok;
        else
            gatesRun = 'G1-G3';  gatesOk = rep.G1.ok && rep.G2.ok && rep.G3.ok;
        end
        pp = paperMap(M);
        pb = antipodal_pair_bound(x0, xT, o.T, geom.r);
        rows(end+1, :) = {M, mode, o.N, o.JLB, pb, o.J, 100*(o.J - o.JLB)/o.JLB, pp(1), pp(2), ...
                          rep.G1.minHfine, o.maxDefect, rep.G4.relChange, mirrorRel, mirrorOk, ...
                          o.iterations, o.wall, o.ipoptStatus, nTried, seedLab, gatesRun, gatesOk};
        lg('M = %d, mode %s: N = %d, J = %.7f (J_LB %.4f, pairBound %.7f, OCFormer %.4f), %s, seed %s of %d, %s %s, mirror %s (%.2e), wall %.1f s', ...
           M, mode, o.N, o.J, o.JLB, pb, pp(1), o.ipoptStatus, seedLab, nTried, gatesRun, mat2str(gatesOk), ...
           mat2str(mirrorOk), mirrorRel, o.wall);
        if strcmp(mode, 'barrier')
            hB  = o.T/o.N;
            pen = hB*sum(o.Unodes) - (hB/2)*(o.Unodes(1) + o.Unodes(end));
            if o.minHnodes > geom.l, branch = 'log'; else, branch = 'quadratic'; end
            lg('M = %d, mode barrier: J = %.7f = smooth %.7f + penalty %.3e (%.2f %% of J); min nodal h = %.3e vs l = %.0e: %s branch', ...
               M, o.J, o.J - pen, pen, 100*pen/o.J, o.minHnodes, geom.l, branch);
        end
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
                 'minHfine', 'maxDefect', 'meshRel', 'mirrorRel', 'mirrorOk', 'iters', 'wallSec', 'status', ...
                 'seedsTried', 'bestSeed', 'gatesRun', 'gatesOk'});
if ~isfolder(fileparts(csvOut)), mkdir(fileparts(csvOut)); end
writetable(tbl, csvOut);
lg('wrote %s', csvOut);
disp(tbl);
end

% ------------------------------------------------------------------------
function tf = adopt_mesh_solution(o, rep)
%% Purpose:
%
%   The adoption decision of the mesh loop: true when the mesh gate ran,
%   its 2N re-solve succeeded and is cheaper than the row by more than
%   1e-12. A cost DROP at 2N is not read as discretisation error: the 2N
%   solve may have found a different (symmetry-broken) basin, and the
%   cheaper of two gated local minima is the number of record.
%
%% Inputs:
%
%  o                        struct                  The row's solve
%                                                   (solve_drones_direct)
%
%  rep                      struct                  Its report
%                                                   (verify_drones_direct)
%
%% Outputs:
%
%  tf                       logical                 Adopt rep.G4.out
%
tf = rep.G4.ran && ~isempty(rep.G4.out) && rep.G4.out.success && rep.G4.out.J < o.J - 1e-12;
end

function seeds = hard_seeds(x0)
%% Purpose:
%
%   The four multistart seeds of a hard row: uniform bump, both senses;
%   alternating amplitude (1.2 on the odd agents of each ring, 0.6 on the
%   even; ring = the agents of one |x0| radius, tolerance 1e-9, in index
%   order), both senses.
%
%% Inputs:
%
%  x0                       [3 x M]                 Start positions
%
%% Outputs:
%
%  seeds                    [1 x 4] struct          .label .sense .scale
%
M   = size(x0, 2);
rad = sqrt(sum(x0.^2, 1));
alt = zeros(1, M);
done = false(1, M);
for a = 1:M
    if done(a), continue, end
    ring = find(~done & abs(rad - rad(a)) < 1e-9);
    sc = 0.6*ones(1, numel(ring));
    sc(1:2:end) = 1.2;                          % odd index within the ring
    alt(ring) = sc;
    done(ring) = true;
end
seeds = struct('label', {'u+', 'u-', 'a+', 'a-'}, 'sense', {1, -1, 1, -1}, ...
               'scale', {ones(1, M), ones(1, M), alt, alt});
end

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
