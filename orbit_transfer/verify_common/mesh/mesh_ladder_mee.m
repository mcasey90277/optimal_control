function L = mesh_ladder_mee(matPath, factors, opts)
% MESH_LADDER_MEE  Re-solve a banked earth MEE row at successively finer meshes.
%
% The study's measuring instrument: takes ONE certified row and re-solves the
% same continuous problem (same t_f, same endpoints, same thrust) on meshes
% scaled by `factors`, warm-chained coarse to fine. What changes between levels
% is ONLY the discretization, so the level-to-level movement in a quantity IS
% its discretization error.
%
% WHAT THIS DELIBERATELY DOES NOT DO: apply the certified-quantity guard.
% refresh_duals_mee gates a re-solve on reproducing the row's certified final
% mass, which is correct there (it re-solves the SAME NLP and a mass change
% would mean something broke). Here the mass is EXPECTED to move with the mesh
% -- that movement is the measurement. Guarding on it would reject exactly the
% signal the study exists to collect. Each level is instead gated on being a
% properly converged solution of its OWN NLP: solver success plus a
% machine-tight defect.
%
% WARM-CHAINING IS DELIBERATE, AND IS ALSO THE MAIN RISK. These basins are
% razor-thin -- the 10 N c_tf case flips between 19 and 24 switches on a 2e-5
% relative change in t_f -- so a cold solve at each level could land on a
% different branch and the "convergence" would be comparing unlike things.
% Chaining keeps the levels on one branch. The cost is that a level can inherit
% an unconverged structure from its predecessor and still pass its own gates:
% per-level success proves each level is stationary for its own NLP, NOT that
% all levels lie on the same continuous branch. That is what the plan's
% separate BRANCH CHECK (an independently seeded solve at the finest level) is
% for; it is not this function's job, and a clean ladder here is not evidence
% against warm-chain bias.
%
% A FAILED LEVEL STOPS THE LADDER. Continuing past a bad rung would feed a
% non-solution into an order fit, which cannot detect it -- mesh_order sees
% only numbers. The failed level is recorded with .ok = false and the remaining
% factors are skipped.
%
% INPUTS:
%   matPath - certified row (.mat): results/MEE_M2_*.mat or *_PSR_psr_final.mat
%   factors - node-count multipliers, coarse to fine [1xL], e.g. [1 2 4 8].
%             Four levels give the sliding-slope consistency check; three give
%             a single fragile slope; two support only a stability statement.
%   opts    - struct (optional):
%             .feasTol     max accepted maxDefect per level    [default 1e-8]
%             .maxIterBase IPOPT cap at factor 1               [default: row's]
%             .maxIterScale cap multiplier per unit factor     [default 1]
%                          (cap at level k = maxIterBase * factor^maxIterScale;
%                           finer meshes need more iterations, and under-
%                           iteration collapses the deep-eps tail -- see
%                           process/DEEP_THRUST_LESSONS.md)
%             .printLevel  IPOPT print level                   [default 0]
%             .verbose     per-level progress to stdout        [default true]
%
% OUTPUTS:
%   L - [1xL] struct array, one entry per attempted level:
%       .factor   the multiplier [scalar]
%       .N        number of INTERVALS at this level [scalar]
%       .sigma    the level's node grid [(N+1)x1]
%       .out      casadi_lt_mee result struct (with .model attached)
%       .rep      foc_check report for this level [struct]
%       .wall     solve wall time, seconds [scalar]
%       .ok       true if the level converged and passed its gates [logical]
%       .why      '' when ok, else the reason it failed [char]
%
% Nothing is saved here. The CALLER writes each level under
% results/mesh_study/MESH_<tag>_x<factor>.mat, per the plan's cache isolation
% constraint -- it must be impossible for a study solve to clobber a
% production .mat.
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/verify/refresh_duals_mee.m (the re-solve
%       pattern this mirrors; differs in mesh scaling and in dropping the
%       certified-mass guard, for the reason in the header above).
%   [2] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Task 3.

if nargin < 3, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
feasTol   = d('feasTol', 1e-8);
mIterScl  = d('maxIterScale', 1);
prnt      = d('printLevel', 0);
verbose   = d('verbose', true);

assert(~isempty(factors) && all(factors > 0), 'mesh_ladder_mee:factors', ...
    'factors must be positive');
assert(issorted(factors), 'mesh_ladder_mee:factors', ...
    'factors must run coarse to fine (warm-chaining depends on the order)');

% --- bootstrap the campaign's own paths from the row's location -------------
% The row lives at <module>/results/<tag>.mat, so the module root is two up.
modRoot = fileparts(fileparts(matPath));
sp = fullfile(modRoot, 'setup_paths.m');
assert(isfile(sp), 'mesh_ladder_mee:paths', ...
    'no setup_paths.m at %s -- is %s really a campaign row?', modRoot, matPath);
here = pwd;  cleaner = onCleanup(@() cd(here)); %#ok<NASGU>
cd(modRoot);  setup_paths();
addpath(fullfile(fileparts(fileparts(modRoot)), 'verify_common'));
setup_verify_common();
addpath(fileparts(mfilename('fullpath')));          % this mesh/ folder

saved   = sosc_load_row(matPath);
par     = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
sigma0  = saved.sigma(:);
N0      = numel(sigma0) - 1;
mIterB  = d('maxIterBase', saved.maxIter);

% The warm start for the first level is the row itself.
Xprev = saved.X;  Uprev = saved.U;  dLprev = saved.dL;  sgPrev = sigma0;

L = struct('factor',{},'N',{},'sigma',{},'out',{},'rep',{},'wall',{},'ok',{},'why',{});

for k = 1:numel(factors)
    fk = factors(k);
    Nk = round(fk * N0);
    sgK = local_rescale_grid(sigma0, Nk);

    % Warm start by interpolating the PREVIOUS level onto this grid. pchip,
    % not spline: pchip is shape-preserving, and a spline overshoots at the
    % throttle's bang-bang steps, seeding thr outside [0,1].
    Xk = interp1(sgPrev, Xprev.', sgK, 'pchip').';
    Uk = interp1(sgPrev, Uprev.', sgK, 'pchip').';
    Uk(1:3,:) = local_unitize(Uk(1:3,:));   % beta must stay a unit direction
    Uk(4,:)   = min(max(Uk(4,:), 0), 1);    % and the throttle inside its box

    sopts = struct('par', par, 'mode', 'fixedtf', 'eps', 0, ...
        'tfTarget', saved.tfTarget, 'x0', saved.X(:,1), 'xf', saved.xf, ...
        'maxIter', round(mIterB * fk^mIterScl), 'warmTight', true, ...
        'printLevel', prnt, 'returnModel', true);

    if verbose
        fprintf('mesh_ladder_mee: x%g  N = %d intervals  (maxIter %d) ... ', ...
                fk, Nk, sopts.maxIter);
    end
    t0 = tic;
    out = casadi_lt_mee(sgK, Xk, Uk, dLprev, sopts);
    wall = toc(t0);

    why = '';
    if ~out.success
        why = sprintf('solver did not succeed (%s)', out.ipoptStatus);
    elseif ~(out.maxDefect <= feasTol)
        why = sprintf('maxDefect %.3e exceeds feasTol %.3e', out.maxDefect, feasTol);
    end
    ok = isempty(why);

    rep = struct();
    if ok
        rep = foc_check(out, sgK, foc_manifest('earth_mee'), struct());
    end

    L(end+1) = struct('factor', fk, 'N', Nk, 'sigma', sgK, 'out', out, ...
                      'rep', rep, 'wall', wall, 'ok', ok, 'why', why); %#ok<AGROW>

    if verbose
        if ok
            fprintf('OK  %.1f s  m_f = %.6f kg  defect %.2e  sw %d\n', ...
                    wall, out.X(6,end), out.maxDefect, out.switches);
        else
            fprintf('FAILED  %.1f s  (%s)\n', wall, why);
        end
    end

    if ~ok
        warning('mesh_ladder_mee:levelFailed', ...
            ['level x%g failed (%s); STOPPING the ladder. Continuing would feed ' ...
             'a non-solution into an order fit, which cannot detect it.'], fk, why);
        break
    end

    % chain forward
    Xprev = out.X;  Uprev = out.U;  dLprev = out.dL;  sgPrev = sgK;
end
end

% ---------------------------------------------------------------------------
function sg = local_rescale_grid(sigma0, N)
% LOCAL_RESCALE_GRID  Resample a node grid to N intervals, PRESERVING its
% relative node distribution.
%
% The mesh's character is part of the method being measured, so refinement
% must not quietly convert a Sundman or PSR grid into a uniform one. Treating
% the original sigma as a monotone map from uniform index-fraction to sigma
% and re-sampling that map at the new resolution reproduces the same relative
% clustering at any N. For an already-uniform sigma this returns linspace.
%
% INPUTS:  sigma0 - original grid [(N0+1)x1]; N - target interval count
% OUTPUTS: sg - the rescaled grid [(N+1)x1]
u0 = linspace(0, 1, numel(sigma0)).';
un = linspace(0, 1, N + 1).';
sg = interp1(u0, sigma0(:), un, 'linear');
sg(1) = sigma0(1);  sg(end) = sigma0(end);      % pin the ends exactly
end

% ---------------------------------------------------------------------------
function B = local_unitize(B)
% LOCAL_UNITIZE  Renormalize columns to unit norm, leaving zero columns alone.
% Interpolating a unit direction does not preserve its norm; the NLP's unit
% constraint would otherwise start violated at every node.
% INPUTS: B [3xM]   OUTPUTS: B [3xM]
n = vecnorm(B, 2, 1);
good = n > 0;
B(:, good) = B(:, good) ./ n(good);
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
