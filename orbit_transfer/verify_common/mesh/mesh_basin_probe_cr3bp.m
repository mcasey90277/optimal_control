function B = mesh_basin_probe_cr3bp(matPath, opts)
% MESH_BASIN_PROBE_CR3BP  Basin probe for the Sundman CR3BP min-fuel campaigns.
%
% The CR3BP analog of mesh_basin_probe (which is MEE/earth-specific). Same
% question, same discipline: re-solve a certified row ON ITS OWN MESH, varying
% only the SEED, so that any improvement is basin selection and cannot be a
% discretization effect. Deliberately a sibling rather than a generalization of
% the earth probe -- the two campaigns share no solver, no row layout, and no
% parameter struct, and merging them would mean a wrapper with two disjoint
% halves.
%
% THE SEED SET DIFFERS FROM THE EARTH PROBE, and the reason is worth stating.
% The earth probe's most productive seeds were down-projections of finer
% mesh-study levels; no mesh ladder has been run for tulip or ELFO, so that
% family is unavailable here. What remains:
%   'row'        the certified solution (the baseline to beat)
%   'shift:+k'   the throttle profile circularly shifted by k nodes. This found
%                a better 10 N basin at k = -1, so it is not a token
%                perturbation. Shifts are given in NODES, and these meshes are
%                ~4000 nodes against the earth rows' 193-1740, so a 1-node
%                shift is a far smaller relative nudge -- hence the wider
%                default ladder.
%   'loose'      the certified seed re-solved with warmTight = false. Not a
%                seed perturbation at all but a different SOLVER PATH: the
%                loose setting uses an adaptive barrier and the default
%                bound_push, which lets IPOPT travel further from the warm
%                start and can leave the basin the tight setting is pinned in.
%
% BEST MEANS HIGHEST FINAL MASS (minimum fuel maximizes m_f), matching the
% one-sided convention in certified_guard(better='higher').
%
% INPUTS:
%   matPath - a certified CR3BP min-fuel .mat carrying sigma, tf, rv0, rvf,
%             tauf0, pSund and out (X [8xN1], U [4xN1]) [char]
%   opts    - struct (optional):
%             .shifts   throttle shifts in nodes [default [-8 -1 1 8]]
%             .doLoose  add the warmTight=false path [default true]
%             .thrustN  thrust for cr3bp_lt_params [default 0.025 N]
%             .m0kg     initial mass  [default 15]
%             .ispS     specific impulse [default 2100]
%             .maxIter  IPOPT cap [default 3000]
%             .verbose  [default true]
%
% OUTPUTS:
%   B - struct array, one per seed: .seed .ok .mf .switches .maxDefect .wall
%       .dMass (mf minus the row's mf, ND mass fraction; positive = better)
%       .why ('' when ok)
%   B(1) is always the 'row' baseline.
%
% SIDE EFFECTS: none. Saves nothing; certified caches are never touched.
%
% REFERENCES:
%   [1] verify_common/doc/mesh_study_tierA_results.md, Finding 6.
%   [2] GTO_tulip/direct/lib/casadi_minfuel_sundman.m (the solver).

if nargin < 2, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
shifts  = d('shifts',  [-8 -1 1 8]);
doLoose = d('doLoose', true);
maxIter = d('maxIter', 3000);
verbose = d('verbose', true);

% paths: the Sundman solver and the CR3BP parameter helpers
ot = fileparts(fileparts(mfilename('fullpath')));      % .../orbit_transfer/verify_common
ot = fileparts(ot);                                    % .../orbit_transfer
addpath(fullfile(ot, 'cr3bp_common'));
addpath(fullfile(ot, 'GTO_tulip', 'direct', 'lib'));

S = load(matPath);
assert(isfield(S,'out') && isfield(S,'sigma'), 'mesh_basin_probe_cr3bp:shape', ...
    '%s lacks out/sigma -- not a Sundman min-fuel row', matPath);
o0    = S.out;
sigma = S.sigma(:);
tf    = S.tf;
rv0   = S.rv0;
rvf   = S.rvf;
tauf0 = local_field(S, 'tauf0', local_field(o0, 'tauf', NaN));
pSund = local_field(S, 'pSund', 1.5);
p     = cr3bp_lt_params(d('thrustN',0.025), d('m0kg',15), d('ispS',2100));
[~, tag] = fileparts(matPath);

sopts = struct('returnModel', false);
mfRow = o0.X(7,end);

% --- seeds ------------------------------------------------------------------
seeds = {struct('name','row', 'X',o0.X, 'U',o0.U, 'tight',true)};
for s = shifts
    U2 = o0.U;
    U2(4,:) = circshift(o0.U(4,:), s);
    seeds{end+1} = struct('name', sprintf('shift:%+d', s), ...
        'X', o0.X, 'U', local_fixU(U2), 'tight', true); %#ok<AGROW>
end
if doLoose
    seeds{end+1} = struct('name','loose', 'X',o0.X, 'U',o0.U, 'tight',false);
end

if verbose
    fprintf('\n=== CR3BP basin probe: %s  (N = %d nodes, tf = %.6f) ===\n', ...
            tag, numel(sigma), tf);
    fprintf('  baseline m_f = %.8f   switches = %d\n', mfRow, local_field(o0,'switches',NaN));
end

B = struct('seed',{},'ok',{},'mf',{},'switches',{},'maxDefect',{},'wall',{}, ...
           'dMass',{},'why',{});
for k = 1:numel(seeds)
    sd = seeds{k};
    t0 = tic;  ok = true;  why = '';  o = struct();
    try
        o = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                sd.X, sd.U, tauf0, pSund, maxIter, 0, sd.tight, sopts);
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
        B(end+1) = struct('seed',sd.name, 'ok',true, 'mf',o.X(7,end), ...
            'switches',o.switches, 'maxDefect',o.maxDefect, 'wall',w, ...
            'dMass', o.X(7,end)-mfRow, 'why',''); %#ok<AGROW>
        if verbose
            fprintf('  %-12s m_f = %.8f  sw = %-4d defect %.1e  (%+.3e)  %.0f s\n', ...
                sd.name, o.X(7,end), o.switches, o.maxDefect, B(end).dMass, w);
        end
    else
        B(end+1) = struct('seed',sd.name, 'ok',false, 'mf',NaN, 'switches',NaN, ...
            'maxDefect',NaN, 'wall',w, 'dMass',NaN, 'why',why); %#ok<AGROW>
        if verbose
            fprintf('  %-12s FAILED (%s)  %.0f s\n', sd.name, why, w);
        end
    end
end

if verbose
    okB = B([B.ok]);
    if isempty(okB)
        fprintf('  ---\n  every seed failed; nothing to compare.\n');
    else
        [~, best] = max([okB.mf]);
        fprintf('  ---\n  BEST: %s  m_f = %.8f  (%+.3e vs the certified row)\n', ...
            okB(best).seed, okB(best).mf, okB(best).dMass);
        if okB(best).dMass > 1e-7
            fprintf('  => the certified row is NOT the best optimum of its own mesh.\n');
        else
            fprintf('  => no seed beat the certified row; it survives this probe.\n');
        end
    end
end
end

% ---------------------------------------------------------------------------
function U = local_fixU(U)
% LOCAL_FIXU  Unit thrust direction, throttle inside its box, so a perturbed
% seed does not start with an avoidable infeasibility.
% INPUTS: U [4xM]   OUTPUTS: U [4xM]
n = vecnorm(U(1:3,:), 2, 1);
good = n > 0;
U(1:3,good) = U(1:3,good) ./ n(good);
U(4,:) = min(max(U(4,:), 0), 1);
end

% ---------------------------------------------------------------------------
function v = local_field(s, f, dflt)
% LOCAL_FIELD  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
