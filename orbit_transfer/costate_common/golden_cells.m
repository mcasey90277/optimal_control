function ok = golden_cells()
%% Purpose:
%
%   GOLDEN-CELL QUALITY REGRESSION for the costate pipeline (principle 7c):
%   fixed benchmark cells with stored reference metrics, where a quality
%   DROP is a failure even when correctness gates still pass. The two
%   subtle bugs this defends against (the Hermite-Simpson midpoint station
%   shift and the missing thrust ratio in continuation mass scaling) both
%   passed every runtime gate because the pipeline is self-healing -- they
%   showed up only in the EFFICIENCY channel (seed quality, iterations).
%
%   Cells:
%     1-3  dro/halo/dpo: ms_tfmin from a flown-then-perturbed seed of one
%          certified 1 N catalog entry per family. Guards the ENGINE
%          (ms_bvp rebinding, Jacobian, guards) + the conjugate test.
%     4    harvest: the full covector path (duals_to_costates ->
%          harvest_ms_seed -> ms_tfmin) on the stored 12x12 reference cell
%          with REAL collocation duals. Guards the HARVEST rules (sign
%          vote, midpoint association, interpolation).
%
%   Checks per cell: solution unchanged (|z - z_ref| < 1e-8), quality
%   unchanged (iters <= iters_ref + 2, normR < 100x normR_ref), conjugate
%   verdict unchanged, and for the harvest cell a unanimous sign vote.
%
%   Run after ANY change to ms_bvp, ms_tfmin, duals_to_costates,
%   harvest_ms_seed, or ms_conjugate_test. Takes ~1 minute.
%
%% Inputs:
%
%   none (references in golden_cells_data.mat beside this file; the
%   harvest cell reads DRO_tulip/direct/results/dsweep_12x12_cells.mat)
%
%% Outputs:
%
%  ok                       logical                 All cells passed
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
G = load(fullfile(here, 'golden_cells_data.mat'));
ok = true;

%% Cells 1-3: engine goldens (flown-then-perturbed seeds):
for kc = 1:numel(G.cells)
    c = G.cells(kc);
    tG = linspace(0, c.z8(8), c.K+1);
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(c.z8(8), [c.rv0(:); 1; c.z8(1:7)], ...
                                       c.Tnd, c.cnd, c.muStar);
    [tu, iu] = unique(tj);
    Y = interp1(tu, yj(iu,1:14), tG, 'pchip')';
    Y(8:14,1) = Y(8:14,1)*(1 + c.pert);
    seed = struct('tf', c.z8(8), 'tGrid', tG, 'Y', Y);
    [z, info] = ms_tfmin(c.rv0, c.rvf, seed, c.Tnd, c.cnd, c.muStar, ...
                         struct('conjTest', true));
    ok = check(c.name, 'converged', info.converged) && ok;
    ok = check(c.name, sprintf('|z-z8| %.1e < 1e-8', norm(z-c.z8)), ...
               norm(z - c.z8) < 1e-8) && ok;
    ok = check(c.name, sprintf('iters %d <= %d', info.iters, c.iters_ref+2), ...
               info.iters <= c.iters_ref + 2) && ok;
    % noise floor at the solver tolerance: a run that converges organically
    % at ~1e-10 is healthy even if the reference was exceptionally deep
    % (review finding, Gemini + GPT 2026-08-08)
    normGate = max(1e-10, 100*c.normR_ref);
    ok = check(c.name, sprintf('normR %.1e < %.1e', info.normR, normGate), ...
               info.normR < normGate) && ok;
    ok = check(c.name, sprintf('conj pass %d == %d', info.conj.pass, ...
               c.conj_pass_ref), info.conj.pass == c.conj_pass_ref) && ok;
end

%% Cell 4: harvest golden (real collocation duals through the full path):
CC = load(fullfile(fileparts(here), 'DRO_tulip', 'direct', 'results', ...
                   'dsweep_12x12_cells.mat'));
cell4 = CC.CELLS{2,5};
o = struct('X', cell4.X(1:7,:), 'lamDef', cell4.lamDef, 'Um', cell4.Um, ...
           'tNodes', cell4.tNodes, 'tf', cell4.tf);
[seed, dg] = harvest_ms_seed(o, 24);
[z, info] = ms_tfmin(cell4.rv0(1:6), cell4.rvf(1:6), seed, ...
                     cell4.Tmax, cell4.c, cell4.muStar, struct());
h = G.harvest;
ok = check('harv', 'converged', info.converged) && ok;
ok = check('harv', sprintf('vote margin %.2f == 1', dg.voteMargin), ...
           dg.voteMargin == 1) && ok;
ok = check('harv', sprintf('lambda_t %.4f ok=%d', dg.lamT, dg.lamTOK), ...
           dg.lamTOK) && ok;
ok = check('harv', sprintf('|z-zref| %.1e < 1e-8', norm(z-h.z_ref)), ...
           norm(z - h.z_ref) < 1e-8) && ok;
ok = check('harv', sprintf('iters %d <= %d', info.iters, h.iters_ref+2), ...
           info.iters <= h.iters_ref + 2) && ok;

if ok, fprintf('GOLDEN CELLS: ALL PASS\n');
else,  fprintf('GOLDEN CELLS: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = check(name, msg, cond)
% CHECK  Report one golden check.  INPUTS: name; msg; cond logical.
% OUTPUTS: ok logical.
ok = logical(cond);
if ok, tag = 'ok  '; else, tag = 'FAIL'; end
fprintf('  [%s] %-4s %s\n', tag, name, msg);
end
