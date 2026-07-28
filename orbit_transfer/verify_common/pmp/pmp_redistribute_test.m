function T = pmp_redistribute_test(matPath, opts)
% PMP_REDISTRIBUTE_TEST  Measure a row's discretization error by perturbing its mesh.
%
% The study's only measurement of discretization error that does not rest on an
% extrapolation. Re-solves a certified row at the SAME node count with the nodes
% redistributed, and reads the error off the resulting mass change.
%
% WHY THIS WORKS WHERE THE LADDER DID NOT. A mesh ladder changes the node COUNT,
% which changes the discrete problem enough that the solver can leave its
% solution branch -- and then the two masses being differenced belong to
% different extremals. Here the node count is FIXED and only the placement
% moves, which is a far milder perturbation. The switch count is checked
% afterwards, and a change in it VOIDS the measurement rather than being
% absorbed into it.
%
% THE REDISTRIBUTION IS DERIVED, NOT CHOSEN. The measured residual obeys
% E ~ C*h^3 per interval (verify_common/doc/pmp_residual_results.md, Finding 3:
% median R_x / h^3 is constant to 5% across three rows whose revolutions,
% switches and node counts each vary ninefold). Minimizing sum h^3 subject to
% fixed total time and fixed node count gives h constant -- a mesh uniform in
% physical TIME. That is the redistribution used.
%
% ITS KNOWN LIMITATION, measured rather than assumed: at 10 N the worst
% interval got 8x WORSE under this redistribution even as the total improved,
% so C is NOT phase-independent and uniform-in-time is not the OPTIMAL mesh.
% It does not need to be. It only needs to be a controlled perturbation of
% known size, which it is.
%
% READING THE ERROR OFF. With two same-branch solutions whose sum h^3 are S1
% and S2 and whose masses differ by dm, and taking E = k*sum h^3,
%
%     k = dm / (S1 - S2),     E_production = k*S1
%
% This is two points and one assumption. The assumption is supported by an
% INDEPENDENT measurement (the R_x/h^3 constancy above), not fitted here, which
% is the difference between this and the Richardson attempts that were
% withdrawn earlier in the study.
%
% INPUTS:
%   matPath - certified row (.mat) [char]
%   opts    - struct (optional):
%             .blend    how far to move toward the uniform-in-time grid, in
%                       [0,1] [default 1.0 = all the way]. Lower values are
%                       milder perturbations, for rows whose branch the full
%                       move does not survive.
%             .maxIter  IPOPT cap for the re-solve [default 3000]
%             .verbose  [default true]
%
% OUTPUTS:
%   T - struct:
%       .ok            false if the switch count changed (measurement VOID)
%       .sw0 .sw1      switch counts before and after
%       .mf0 .mf1      final masses, kg
%       .dm            mf1 - mf0, kg
%       .sumH3_0 .sumH3_1
%       .k             kg per unit sum h^3
%       .errProd       inferred discretization error of the PRODUCTION mesh, kg
%       .errRedist     inferred error of the redistributed mesh, kg
%       .relErr        .errProd / .mf0
%       .hRatio0 .hRatio1   step-size range before and after
%       .why           '' when ok
%
% REFERENCES:
%   [1] pmp_residual_mee.m (produces the residual map this rests on).
%   [2] verify_common/doc/pmp_residual_results.md, Steps 1 and 2.

if nargin < 2, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
verbose = d('verbose', true);

here = fileparts(mfilename('fullpath'));
ot   = fileparts(fileparts(here));
addpath(here);  addpath(fullfile(ot,'verify_common','mesh'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
modRoot = fileparts(fileparts(matPath));
cwd0 = pwd;  cleaner = onCleanup(@() cd(cwd0)); %#ok<NASGU>
cd(modRoot);  setup_paths();
addpath(fullfile(ot,'verify_common'));  setup_verify_common();

[~, tag] = fileparts(matPath);
R0 = pmp_residual_mee(matPath, struct('verbose', false));
saved = sosc_load_row(matPath);
h0  = R0.hPhys;
sg0 = R0.sigma(:);  t0 = R0.X(7,:).';
N1  = numel(sg0);

% --- uniform-in-time grid at the SAME node count ---------------------------
% A BLEND toward uniform-in-time, not necessarily all the way. The full move
% is the optimal one for sum h^3, but it is also a large perturbation, and the
% deeper rows are basin-fragile enough that it knocks them onto another branch
% (measured: 76 -> 75 switches at 2.5 N, 171 -> 169 at 1 N, both VOID). A
% partial blend is a smaller perturbation that may hold the branch while still
% changing sum h^3 enough to read the error off. Any blend works for the
% measurement -- it only has to be a controlled change of known size.
blend = d('blend', 1.0);
tU  = linspace(t0(1), t0(end), N1).';
sgU = interp1(t0, sg0, tU, 'pchip');
sgN = (1-blend)*sg0 + blend*sgU;
sgN(1) = sg0(1);  sgN(end) = sg0(end);
sgN = sort(sgN);
T_blend = blend;

Xn = interp1(sg0, R0.X.', sgN, 'pchip').';
Un = interp1(sg0, R0.U.', sgN, 'pchip').';
nn = vecnorm(Un(1:3,:),2,1);  Un(1:3,:) = Un(1:3,:) ./ max(nn, eps);
Un(4,:) = min(max(Un(4,:),0),1);

par = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
sopts = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
   'x0',saved.X(:,1),'xf',saved.xf,'maxIter',d('maxIter',3000),'warmTight',true, ...
   'printLevel',0,'returnModel',false);

if verbose
    fprintf('\n=== redistribute test: %s (N = %d, uniform in time) ===\n', tag, N1-1);
end
tic;  o = casadi_lt_mee(sgN, Xn, Un, saved.dL, sopts);  w = toc;

h1 = diff(o.X(7,:));
T = struct();
T.sw0 = numel(R0.switchTimes);  T.sw1 = o.switches;
T.mf0 = R0.X(6,end) * saved.m0kg;
T.mf1 = o.m_f_kg;
T.dm  = T.mf1 - T.mf0;
T.sumH3_0 = sum(h0.^3);  T.sumH3_1 = sum(h1.^3);
T.hRatio0 = max(h0)/min(h0);  T.hRatio1 = max(h1)/min(h1);
T.wall = w;  T.ipoptStatus = o.ipoptStatus;  T.maxDefect = o.maxDefect;
T.blend = T_blend;

T.why = '';
if ~o.success || o.maxDefect > 1e-8
    T.ok = false;  T.why = sprintf('re-solve failed (%s, defect %.2e)', o.ipoptStatus, o.maxDefect);
elseif T.sw1 ~= T.sw0
    % A topology change makes the two masses belong to different extremals, so
    % differencing them measures the branch, not the mesh. VOID, not a result.
    T.ok = false;
    T.why = sprintf('switch count changed %d -> %d: different branch, measurement VOID', ...
                    T.sw0, T.sw1);
else
    T.ok = true;
end

if T.ok
    T.k        = T.dm / (T.sumH3_0 - T.sumH3_1);
    T.errProd  = T.k * T.sumH3_0;
    T.errRedist= T.k * T.sumH3_1;
    T.relErr   = abs(T.errProd) / T.mf0;
else
    T.k = NaN;  T.errProd = NaN;  T.errRedist = NaN;  T.relErr = NaN;
end

if verbose
    fprintf('  %s  %.1f s   h ratio %.1f -> %.1f   sum h^3 %.4g -> %.4g (%.2fx)\n', ...
        o.ipoptStatus, w, T.hRatio0, T.hRatio1, T.sumH3_0, T.sumH3_1, T.sumH3_0/T.sumH3_1);
    fprintf('  switches %d -> %d   m_f %.6f -> %.6f  (%+.6f kg)\n', ...
        T.sw0, T.sw1, T.mf0, T.mf1, T.dm);
    if T.ok
        fprintf('  => production-mesh error %.4f kg  (%.2e relative); redistributed %.4f kg\n', ...
            T.errProd, T.relErr, T.errRedist);
    else
        fprintf('  => VOID: %s\n', T.why);
    end
end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
