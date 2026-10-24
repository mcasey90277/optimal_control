function rep = verify_drones_direct(out, opts)
%% Purpose:
%
%   The gates on one direct solve, with its own verdict:
%     G1 clearance   min over pairs (and the ball) of the gap on a nSub-times
%                    finer cubic-Hermite reconstruction from (X, V) -- a
%                    nodal constraint says nothing between nodes, and the
%                    interpolant cuts a constraint arc's corner by up to
%                    s^2/(16 r), so the bar is mesh-aware (-1e-2 * 2r) and
%                    the number is REPORTED, not hidden
%     G2 bound       J >= J_LB
%     G3 solver      IPOPT status ok (out.success) and the trapezoid
%                    defects, re-evaluated here from out.X, out.V (not read
%                    from out.maxDefect), < 1e-9
%     G4 mesh        (opts.mesh) re-solve at 2N warm-started from the
%                    interpolant; relative cost change
%     G5 mirror      (opts.mirror) re-solve with the opposite swirl; costs
%                    agree -- only meaningful on a mirror-symmetric layout
%                    (the paper's circles/rings and headon); skipped with
%                    a note otherwise; a G5 that did not run has
%                    .ran = false, .ok = false, .relDiff = NaN and does
%                    not enter rep.ok (a G5 that ran and failed does)
%
%% Inputs:
%
%  out                      struct                  From solve_drones_direct
%
%  opts                     struct (optional)       .mirror [true] .mesh
%                                                   [false] .nSub [10]
%                                                   .tolG1 [1e-2 * 2r]
%                                                   .tolDefect [1e-9]
%                                                   .tolMesh [1e-3]
%                                                   .tolMirror [1e-6]
%                                                   .mirrorX0 [[]] .mirrorV0
%                                                   [[]] warm start for the
%                                                   G5 re-solve (X0/V0 of
%                                                   solve_drones_direct, at
%                                                   the row's N); used only
%                                                   when both are non-empty
%
%% Outputs:
%
%  rep                      struct                  .ok .why .G1 .G2 .G3
%                                                   .G4 .G5 (.out = mirror solve) .notes
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
    casadi_path;
    [x0, xT, geom] = drone_layouts('headon');
    out = solve_drones_direct(x0, xT, geom, struct('N', 50));
    rep = verify_drones_direct(out);
    fprintf('verify_drones_direct demo (headon M = 2, N = 50): ok = %d, why = ''%s''\n', rep.ok, rep.why);
    fprintf('  G1 minHfine = %.3e (bar %.1e)   G2 gap = %.3e   G3 %s, maxDefect = %.1e\n', ...
            rep.G1.minHfine, rep.G1.bar, rep.G2.gap, rep.G3.status, rep.G3.maxDefect);
    return
end
if nargin < 2, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
d = @(f, v) fieldd(opts, f, v);
geom = out.geom;
doMirror = d('mirror', true);
doMesh   = d('mesh', false);
nSub     = d('nSub', 10);
tolG1    = d('tolG1', 1e-2*2*geom.r);
tolDef   = d('tolDefect', 1e-9);
tolMesh  = d('tolMesh', 1e-3);
tolMir   = d('tolMirror', 1e-6);
mirX0    = d('mirrorX0', []);
mirV0    = d('mirrorV0', []);
rep = struct('ok', false, 'why', '', 'notes', {{}});

%% G3 first: a failed solve has no gates worth reading; the defects are
% re-evaluated from the stored X, V
h   = out.T/out.N;
def = out.X(:, 2:end) - out.X(:, 1:end-1) - (h/2)*(out.V(:, 2:end) + out.V(:, 1:end-1));
maxDef = max(abs(def(:)));
rep.G3 = struct('status', out.ipoptStatus, 'maxDefect', maxDef, ...
                'ok', out.success && maxDef < tolDef);

%% G1: the fine grid
[Xf, ~] = hermite_resample(out.X, out.V, out.tN, nSub);
if size(out.pairs, 2) > 0
    [hf, ~, ~, hbf] = clearance_terms(Xf, out.pairs, geom);
    minH = min(hf(:));
else
    [~, ~, ~, hbf] = clearance_terms(Xf, out.pairs, geom);
    minH = Inf;
end
if isempty(hbf), minHb = Inf; else, minHb = min(hbf(:)); end
rep.G1 = struct('minHfine', minH, 'minHballFine', minHb, 'nSub', nSub, 'bar', -tolG1);
switch out.mode
    case 'hard',    rep.G1.ok = minH >= -tolG1 && minHb >= -tolG1;
    case 'barrier', rep.G1.ok = true;   % reported, not gated: the barrier is the constraint
end

%% G2
rep.G2 = struct('gap', out.J - out.JLB, 'ok', out.J >= out.JLB - 1e-12);

%% G4: mesh
rep.G4 = struct('ran', false, 'J2N', NaN, 'relChange', NaN, 'ok', true);
if doMesh && out.success
    [X2, V2] = hermite_resample(out.X, out.V, out.tN, 2);
    o2 = solve_drones_direct(reshape(out.X(:,1), 3, geom.M), reshape(out.X(:,end), 3, geom.M), geom, ...
                             struct('N', 2*out.N, 'T', out.T, 'mode', out.mode, 'sense', out.sense, ...
                                    'bump', out.bump, 'pairs', out.pairs, 'midpoints', out.midpoints, ...
                                    'X0', X2, 'V0', V2));
    rep.G4.ran = true;  rep.G4.J2N = o2.J;
    rep.G4.relChange = abs(o2.J - out.J)/abs(out.J);
    rep.G4.ok = o2.success && rep.G4.relChange < tolMesh;
end

%% G5: mirror
rep.G5 = struct('ran', false, 'Jmirror', NaN, 'relDiff', NaN, 'ok', false, 'note', '', 'status', '', 'success', true, 'out', []);
symmetric = any(strcmp(geom.which, {'paper', 'circle', 'rings', 'headon'}));
if doMirror && out.success
    if symmetric
        % a warm-started row must be mirrored against a warm-started mirror:
        % with opts.mirrorX0/V0 the re-solve starts from them, else cold
        mopts = struct('N', out.N, 'T', out.T, 'mode', out.mode, 'sense', -out.sense, ...
                       'bump', out.bump, 'pairs', out.pairs, 'midpoints', out.midpoints);
        if ~isempty(mirX0) && ~isempty(mirV0)
            mopts.X0 = mirX0;  mopts.V0 = mirV0;
        end
        try
            om = solve_drones_direct(reshape(out.X(:,1), 3, geom.M), reshape(out.X(:,end), 3, geom.M), geom, mopts);
            rep.G5.out = om;  rep.G5.ran = true;  rep.G5.Jmirror = om.J;  rep.G5.status = om.ipoptStatus;  rep.G5.success = om.success;
            rep.G5.relDiff = abs(om.J - out.J)/abs(out.J);
            rep.G5.ok = om.success && rep.G5.relDiff < tolMir;
        catch err
            if ~strcmp(err.identifier, 'solve_drones_direct:guess'), rethrow(err); end
            rep.G5.note = 'barrier mirror needs a warm start (default guess has coincident agents): G5 skipped';
            rep.notes{end+1} = rep.G5.note;
        end
    else
        rep.G5.note = sprintf('layout ''%s'' is not mirror-symmetric: G5 skipped', geom.which);
        rep.notes{end+1} = rep.G5.note;
    end
end

%% Verdict
if ~rep.G3.ok
    rep.why = sprintf('G3: IPOPT %s, max defect %.2e', out.ipoptStatus, maxDef);
elseif ~rep.G1.ok
    rep.why = sprintf('G1: fine-grid clearance %.2e below %.1e', min(minH, minHb), -tolG1);
elseif ~rep.G2.ok
    rep.why = sprintf('G2: J = %.6f below J_LB = %.6f', out.J, out.JLB);
elseif ~rep.G4.ok
    rep.why = sprintf('G4: cost changed %.2e relative at 2N', rep.G4.relChange);
elseif rep.G5.ran && ~rep.G5.ok && ~rep.G5.success
    rep.why = sprintf('G5: mirror re-solve failed (%s)', rep.G5.status);
elseif rep.G5.ran && ~rep.G5.ok
    rep.why = sprintf('G5: mirror cost differs by %.2e relative', rep.G5.relDiff);
end
rep.ok = isempty(rep.why);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
