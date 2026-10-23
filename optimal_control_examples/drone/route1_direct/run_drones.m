%% RUN_DRONES  Solve one multi-agent drone problem directly, certify it, plot it.
%
%  Purpose:
%
%   Front door of Route 1 (CasADi/IPOPT direct transcription of the paper's
%   multi-agent drone problem). Choose a layout, the number of drones, the
%   clearance mode and the mesh in section 1; the script solves, runs the
%   five verification gates, prints a summary against the rigorous
%   antipodal-pair bound and (for the 'paper' layout) the paper's Table 1
%   numbers, and writes a four-snapshot PNG and an MP4. The paper's numbers
%   are printed for comparison only and never gate anything.
%
%  Outputs (in results/ beside this script):
%
%   <label>_panels.png     configuration at t = 0, T/3, 2T/3, T
%   <label>.mp4            1280 x 720 movie
%
%  Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% 1. USER SETTINGS ------------------------------------------------------
label  = 'drones_M4';
which  = 'paper';      % 'paper' | 'circle' | 'rings' | 'sphere' | 'headon'   (drone_layouts)
M      = 4;            % number of drones ('headon' ignores it: 2)
r      = [];           % drone radius (m); [] = the layout's default (0.02, 0.01 above M = 32)
mode   = 'hard';       % 'hard' (clearance as a constraint) | 'barrier' (the paper's log barrier in the cost)
N      = 100;          % trapezoid intervals
T      = 10;           % horizon (s)
sense  = +1;           % swirl sense of the initial guess (+1 / -1)
mesh   = true;         % run the mesh gate (re-solve at 2N); costs a second solve
movieOpts = struct('fps', 30, 'nFrames', 300);
%% -------------------------------------------------------------------------

%% 2. set up and solve
here = fileparts(mfilename('fullpath'));
casadi_path;
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
fprintf('paths set (casadi + route1_direct + drone_dynamics)\n');
[x0, xT, geom] = drone_layouts(which, M, r);
M = geom.M;
out = solve_drones_direct(x0, xT, geom, struct('mode', mode, 'N', N, 'T', T, 'sense', sense));
fprintf('solved: %s, status %s, %d iterations, %.1f s\n', out.mode, out.ipoptStatus, out.iterations, out.wall);

%% 3. verify
rep = verify_drones_direct(out, struct('mirror', true, 'mesh', mesh));
fprintf('verified: ok = %s\n', mat2str(rep.ok));

%% 4. summary
pb = antipodal_pair_bound(x0, xT, T, geom.r);
paperTab = [4 0.2003 0.2001; 8 0.4036 0.4046; 16 0.8307 0.8345; 32 1.2531 1.4413; 64 2.448 NaN; 128 5.1424 13.3];
kp = find(paperTab(:, 1) == M, 1);
fprintf('\n%s: layout %s, M = %d, mode %s, N = %d\n', label, which, M, out.mode, out.N);
fprintf('  J = %.6f   J_LB = %.6f   gap %.3f %%\n', out.J, out.JLB, 100*(out.J - out.JLB)/out.JLB);
fprintf('  rigorous antipodal-pair bound = %s\n', numOrNA(pb, '%.6f'));
if strcmp(which, 'paper') && ~isempty(kp)
    fprintf('  paper Table 1: OCFormer %s   SympOCNet %s\n', numOrNA(paperTab(kp, 2), '%.4f'), numOrNA(paperTab(kp, 3), '%.4f'));
else
    fprintf('  paper Table 1: n/a\n');
end
fprintf('  G1 fine-grid clearance  min h = %+.3e   %s\n', rep.G1.minHfine, verdict(rep.G1.ok));
fprintf('  G2 cost vs bound        J - J_LB = %+.3e   %s\n', rep.G2.gap, verdict(rep.G2.ok));
fprintf('  G3 solver               status %s, maxDefect = %.3e   %s\n', rep.G3.status, rep.G3.maxDefect, verdict(rep.G3.ok));
if rep.G4.ran, fprintf('  G4 mesh (2N)            relChange = %.3e   %s\n', rep.G4.relChange, verdict(rep.G4.ok));
else,          fprintf('  G4 mesh (2N)            not run\n');
end
if rep.G5.ran, fprintf('  G5 mirror               relDiff = %.3e   %s\n', rep.G5.relDiff, verdict(rep.G5.ok));
else,          fprintf('  G5 mirror               not run (%s)\n', rep.G5.note);
end
fprintf('  gates ok = %s   (%s)\n', mat2str(rep.ok), rep.why);

%% 5. plot and film
outDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
pngFile = fullfile(outDir, [label '_panels.png']);
fig = plot_drones(out, struct('visible', false, 'file', pngFile));
close(fig);
movie_drones(out, fullfile(outDir, label), movieOpts);
fprintf('  wrote %s\n        %s\n', pngFile, fullfile(outDir, [label '.mp4']));

%% -------------------------------------------------------------------------
function pb = antipodal_pair_bound(x0, xT, T, r)
%% Purpose:
%
%   Rigorous lower bound on the hard-clearance cost, as in drone_table1:
%   each antipodal pair (agent i starts where agent j ends) is the
%   two-agent head-on problem, costing at least ell^2/(4T) with ell = 2
%   sqrt(D^2 - R^2) + R (pi - 2 acos(R/D)), D the pair distance, R = 2r.
%   Each pair is counted once.
%
R  = 2*r;
pb = 0;
for ia = 1:size(x0, 2)
    for ib = ia+1:size(x0, 2)
        if max(abs(x0(:, ib) - xT(:, ia))) < 1e-12
            D   = sqrt(sum((x0(:, ia) - x0(:, ib)).^2));
            ell = 2*sqrt(D^2 - R^2) + R*(pi - 2*acos(R/D));
            pb  = pb + ell^2/(4*T);
        end
    end
end
end

% ------------------------------------------------------------------------
function s = numOrNA(x, fmt)
%% Purpose:
%
%   Format x, or 'n/a' when it is NaN or zero-valued by absence.
%
if isnan(x) || x == 0, s = 'n/a'; else, s = sprintf(fmt, x); end
end

% ------------------------------------------------------------------------
function s = verdict(ok)
%% Purpose:
%
%   PASS / FAIL text.
%
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
