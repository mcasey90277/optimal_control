function V = verify_with_pumpkyn(T, B, opts)
%% Purpose:
%
%   INDEPENDENT VERIFICATION, made visible: hand our converged costates to
%   pumpkyn's own minimum-time solver and show, component by component, that
%   it does not move them. Then fly ITS answer and check that it too reaches
%   the target.
%
%   Why both halves. Agreement alone is not convergence -- a solver that
%   returned its input on stagnation would report a perfect zero. This one
%   demonstrably does not (measured 2026-09-10: perturbing lambda_0 by
%   1.5x, 3x and 10x moved its answer by 9.3, 37 and 4612), and the test
%   beside this file re-runs that control experiment. Flying its solution
%   closes the question regardless: two independently obtained solutions
%   that both reach the target is a stronger statement than two vectors
%   that happen to match.
%
%   This is NOT an optimality condition. It is a guard against our own
%   solver: a bug in our shooting could produce a self-consistent answer to
%   the wrong problem, and only a second implementation catches that.
%
%% Inputs:
%
%  T                        struct                  certify_root output
%                                                   (.z, .sD, .sA)
%  B                        struct                  arclength_arrival setup
%  opts                     struct (optional)
%   .quiet [false] .tolDz [1e-6] .gateKm [100] .capSec [300] .pool
%
%% Outputs:
%
%  V                        struct                  .z8 .z8Pumpkyn
%                                                   .dComponent [8x1] .dz
%                                                   .moved .converged
%                                                   .flyKm .note
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
tolDz = d('tolDz', 1e-6);  gateKm = d('gateKm', 100);  capSec = d('capSec', 300);
pool = d('pool', gcp('nocreate'));
lStar = B.problem.lStar;

rv0 = B.stateD(T.sD);  rvf = B.stateA(T.sA);
z8 = T.z(:);
V = struct('z8', z8, 'z8Pumpkyn', nan(8,1), 'dComponent', nan(8,1), 'dz', NaN, ...
           'moved', true, 'converged', false, 'flyKm', NaN, 'note', '');

% the foreign solver prints its own fsolve banner; capture it so this
% section reads as one comparison rather than a solver log
[okW, za] = fenced(pool, capSec, @quietTfMin, 1, rv0(1:6)', rvf(1:6)', z8, ...
                   B.Tnd, B.cnd, B.mu);
if ~okW
    V.note = sprintf('pumpkyn.cr3bp.tfMin exceeded its %g s cap', capSec);
elseif ~(isnumeric(za) && numel(za) == 8 && all(isfinite(za(:))))
    V.note = 'pumpkyn.cr3bp.tfMin returned an unusable vector';
else
    V.converged = true;
    V.z8Pumpkyn = za(:);
    V.dComponent = za(:) - z8;
    V.dz = norm(V.dComponent);
    V.moved = ~(V.dz <= tolDz);
    [okF, ~, Ya] = fenced(pool, capSec, @pumpkyn.cr3bp.tfMinProp, 2, za(8), ...
                          [rv0(1:6); 1; za(1:7)], B.Tnd, B.cnd, B.mu);
    if okF, V.flyKm = norm(Ya(end,1:3) - rvf(1:3)')*lStar; end
    V.note = sprintf('|dz| = %.3e (tolerance %.0e); its own flight misses by %.4f km', ...
                     V.dz, tolDz, V.flyKm);
end

if ~d('quiet', false)
    nm = {'lambda_r1','lambda_r2','lambda_r3','lambda_v1','lambda_v2','lambda_v3','lambda_m','t_f'};
    fprintf('\n--- INDEPENDENT VERIFICATION (pumpkyn.cr3bp.tfMin) ---\n');
    fprintf('  our costates handed to a second, independently written solver;\n');
    fprintf('  if ours solve the same problem, it has nothing to change.\n\n');
    fprintf('  %-11s %18s %18s %14s\n', 'component', 'ours', 'pumpkyn''s', 'difference');
    for k = 1:8
        fprintf('  %-11s %18.10g %18.10g %14.2e\n', nm{k}, z8(k), V.z8Pumpkyn(k), V.dComponent(k));
    end
    fprintf('\n  |dz| = %.3e     %s\n', V.dz, tern(~V.moved, 'COSTATES DID NOT MOVE', 'COSTATES MOVED'));
    if isfinite(V.flyKm)
        fprintf('  pumpkyn''s own solution flies to the target within %.4f km  %s\n', ...
                V.flyKm, tern(V.flyKm < gateKm, '(PASS)', '(FAIL)'));
    end
    fprintf('------------------------------------------------------\n');
end
end

function za = quietTfMin(rv0r, rvfr, z8, Tnd, cnd, mu)
% QUIETTFMIN  pumpkyn.cr3bp.tfMin with its console output captured.
% INPUTS: rv0r; rvfr; z8; Tnd; cnd; mu.  OUTPUTS: za [8x1].
evalc('za = pumpkyn.cr3bp.tfMin(rv0r, rvfr, z8, Tnd, cnd, mu);');
end

function varargout = fenced(pool, capSec, fh, nout, varargin)
% FENCED  One external call under a hard cap when a pool is available.
% INPUTS: pool; capSec; fh; nout; varargin.  OUTPUTS: ok, then nout outputs.
varargout = cell(1, nout + 1);
if isempty(pool)
    try
        [varargout{2:nout+1}] = feval(fh, varargin{:});  varargout{1} = true;
    catch
        varargout{1} = false;
    end
    return
end
[varargout{1}, varargout{2:nout+1}] = run_capped(pool, fh, nout, capSec, varargin{:});
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
