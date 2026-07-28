function B = mesh_basin_probe(matPath, opts)
% MESH_BASIN_PROBE  Re-solve a certified row on ITS OWN mesh from several seeds.
%
% WHAT QUESTION THIS ANSWERS, AND WHY IT IS NOT A MESH QUESTION. The mesh study
% asks "is the production mesh fine enough?". This asks a different and, at
% 10 N, a bigger question: "is the certified row even the best solution of its
% OWN discretization?" Those are independent -- a row can be perfectly
% mesh-converged and still sit in an inferior local minimum, because the
% min-fuel bang-bang problem is non-convex and its basins are documented to be
% razor-thin (the 10 N c_tf case flips 19 <-> 24 switches on a 2e-5 relative
% change in t_f).
%
% EVERY SOLVE HERE USES THE ROW'S OWN PRODUCTION MESH. Only the SEED varies.
% That is what isolates basin selection from resolution: if a seed finds a
% better optimum on the identical grid, the improvement cannot be a
% discretization effect.
%
% SEEDS. The productive ones are not random perturbations -- at eps = 0 a
% perturbed seed usually falls straight back into the same basin. They are:
%   'row'        the certified solution itself (the baseline to beat)
%   'downproj:F' the mesh-study level at factor F, interpolated DOWN onto the
%                production grid. This is the seed family that found the better
%                10 N basin: solving on a finer mesh lets the solver escape,
%                and the escaped structure then survives projection back down.
%   'shift:+k'   the row's throttle profile shifted by k nodes, states kept.
%                A cheap, structured nudge of switch placement -- unlike random
%                noise it moves the SWITCH TIMES, which is the coordinate the
%                basins actually differ in.
%
% BEST MEANS HIGHEST FINAL MASS. Minimum fuel maximizes m_f, so a higher mass
% at machine-tight feasibility is strictly better, never a failure. This
% mirrors the one-sided convention in reproduce/verify_row.m and
% certified_guard's better='higher'.
%
% INPUTS:
%   matPath - certified row (.mat) [char]
%   opts    - struct (optional):
%             .factors  down-projection seeds to try [default [2 4 8]]
%             .shifts   throttle-shift seeds, in nodes [default [-1 1]]
%             .meshDir  where the MESH_*_x*.mat levels live
%                       [default <modRoot>/results/mesh_study]
%             .maxIter  IPOPT cap [default: the row's own]
%             .verbose  [default true]
%
% OUTPUTS:
%   B - struct array, one per seed attempted:
%       .seed [char] .ok .m_f_kg .switches .dV_kms .maxDefect .wall
%       .dMass  m_f minus the certified row's m_f (kg; positive = better)
%       .why    '' when ok
%   B(1) is always the 'row' baseline. Sorted output is left to the caller so
%   the baseline's position stays identifiable.
%
% SIDE EFFECTS: none. This function saves nothing -- it is a probe, and the
% certified caches are never touched.
%
% REFERENCES:
%   [1] verify_common/doc/mesh_study_tierA_results.md, Finding 5 (the
%       down-projection test that motivated this).
%   [2] verify/refresh_duals_mee.m (the re-solve pattern; same solver options).

if nargin < 2, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
factors = d('factors', [2 4 8]);
shifts  = d('shifts',  [-1 1]);
verbose = d('verbose', true);

modRoot = fileparts(fileparts(matPath));
here = pwd;  cleaner = onCleanup(@() cd(here)); %#ok<NASGU>
cd(modRoot);  setup_paths();
addpath(fullfile(fileparts(fileparts(modRoot)), 'verify_common'));
setup_verify_common();
meshDir = d('meshDir', fullfile(modRoot, 'results', 'mesh_study'));

saved = sosc_load_row(matPath);
par   = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
sgP   = saved.sigma(:);                       % the PRODUCTION grid, fixed
mIter = d('maxIter', saved.maxIter);
[~, tag] = fileparts(matPath);

sopts = struct('par', par, 'mode', 'fixedtf', 'eps', 0, ...
    'tfTarget', saved.tfTarget, 'x0', saved.X(:,1), 'xf', saved.xf, ...
    'maxIter', mIter, 'warmTight', true, 'printLevel', 0, 'returnModel', false);

% --- assemble the seed list -------------------------------------------------
seeds = {struct('name','row', 'X',saved.X, 'U',saved.U, 'dL',saved.dL)};
for f = factors
    lf = fullfile(meshDir, sprintf('MESH_%s_x%g.mat', tag, f));
    if ~isfile(lf), continue; end
    Lv = load(lf);  lv = Lv.lvl;
    if ~lv.ok, continue; end
    Xd = interp1(lv.sigma(:), lv.out.X.', sgP, 'pchip').';
    Ud = interp1(lv.sigma(:), lv.out.U.', sgP, 'pchip').';
    seeds{end+1} = struct('name', sprintf('downproj:x%g', f), ...
        'X', Xd, 'U', local_fixU(Ud), 'dL', lv.out.dL); %#ok<AGROW>
end
for s = shifts
    U2 = saved.U;
    U2(4,:) = circshift(saved.U(4,:), s);
    seeds{end+1} = struct('name', sprintf('shift:%+d', s), ...
        'X', saved.X, 'U', local_fixU(U2), 'dL', saved.dL); %#ok<AGROW>
end

% --- solve from each seed on the SAME grid ----------------------------------
mfRow = saved.X(6,end) * saved.m0kg;
B = struct('seed',{},'ok',{},'m_f_kg',{},'switches',{},'dV_kms',{}, ...
           'maxDefect',{},'wall',{},'dMass',{},'why',{});
if verbose
    fprintf('\n=== basin probe: %s  (production mesh, N = %d) ===\n', tag, numel(sgP)-1);
end
for k = 1:numel(seeds)
    sd = seeds{k};
    t0 = tic;
    ok = true;  why = '';  o = struct();
    try
        o = casadi_lt_mee(sgP, sd.X, sd.U, sd.dL, sopts);
        if ~o.success
            ok = false;  why = o.ipoptStatus;
        elseif ~(o.maxDefect <= 1e-8)
            ok = false;  why = sprintf('defect %.2e', o.maxDefect);
        end
    catch err
        ok = false;  why = err.message;
    end
    w = toc(t0);
    if ok
        B(end+1) = struct('seed', sd.name, 'ok', true, 'm_f_kg', o.m_f_kg, ...
            'switches', o.switches, 'dV_kms', o.dV_kms, 'maxDefect', o.maxDefect, ...
            'wall', w, 'dMass', o.m_f_kg - mfRow, 'why', ''); %#ok<AGROW>
        if verbose
            fprintf('  %-14s m_f = %12.6f kg  sw = %-4d dV = %.6f  (%+.6f kg)  %.1f s\n', ...
                sd.name, o.m_f_kg, o.switches, o.dV_kms, B(end).dMass, w);
        end
    else
        B(end+1) = struct('seed', sd.name, 'ok', false, 'm_f_kg', NaN, ...
            'switches', NaN, 'dV_kms', NaN, 'maxDefect', NaN, 'wall', w, ...
            'dMass', NaN, 'why', why); %#ok<AGROW>
        if verbose
            fprintf('  %-14s FAILED (%s)  %.1f s\n', sd.name, why, w);
        end
    end
end

if verbose
    okB = B([B.ok]);
    [~, best] = max([okB.m_f_kg]);
    fprintf('  ---\n  BEST: %s  %.6f kg  (%+.6f kg vs the certified row, %+.3e rel)\n', ...
        okB(best).seed, okB(best).m_f_kg, okB(best).dMass, okB(best).dMass/mfRow);
    if okB(best).dMass > 1e-4
        fprintf('  => the certified row is NOT the best optimum of its own mesh.\n');
    else
        fprintf('  => no seed beat the certified row; it survives this probe.\n');
    end
end
end

% ---------------------------------------------------------------------------
function U = local_fixU(U)
% LOCAL_FIXU  Make an interpolated/shifted control admissible as a warm start:
% unit thrust direction, throttle inside its box. Interpolation does not
% preserve either, and a seed that starts outside the feasible box wastes
% iterations shedding an avoidable infeasibility.
% INPUTS: U [4xM]   OUTPUTS: U [4xM]
n = vecnorm(U(1:3,:), 2, 1);
good = n > 0;
U(1:3, good) = U(1:3, good) ./ n(good);
U(4,:) = min(max(U(4,:), 0), 1);
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
