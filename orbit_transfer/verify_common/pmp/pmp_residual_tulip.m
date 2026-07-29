function R = pmp_residual_tulip(matPath, opts)
% PMP_RESIDUAL_TULIP  Continuous PMP residuals for a certified GTO->tulip row.
%
% The Sundman-CR3BP sibling of pmp_residual_mee. Same question, same
% discipline: does the converged direct solution, read as a CONTINUOUS object,
% satisfy the continuous Pontryagin necessary conditions, and to what residual?
%
% TWO STRUCTURAL DIFFERENCES FROM THE MEE ADAPTER, both in this campaign's
% favour:
%
%   1. AUTONOMOUS. The Sundman clock is kappa = r1^p with r1 a function of the
%      STATE, so the regularized dynamics carry no explicit dependence on the
%      independent variable tau. The Hamiltonian H_tau is therefore CONSTANT
%      along an extremal -- a genuine conserved quantity, and a stronger check
%      than the MEE campaign can offer, where the Gauss equations carry sigma
%      explicitly through L = pi + sigma*DeltaL and only lam_t is conserved.
%      Both are reported here.
%   2. THE DYNAMICS MUST BE REBUILT, not reused. casadi_minfuel_sundman
%      constructs its right-hand side INLINE and exports no reusable function,
%      unlike lt_mee_rhs. The rebuild below mirrors that block line for line,
%      INCLUDING the 1e-12 guards inside the distance norms -- omit them and
%      the two right-hand sides differ at perigee, where it matters most.
%
%      This is exactly the situation the self-validating gate exists for. It
%      recomputes the transcription's own trapezoidal defects from the rebuilt
%      dynamics and compares them against the solver's reported maxDefect. If
%      the rebuild is wrong the gate fails and the run aborts, because a
%      residual computed with the wrong dynamics is worse than no residual --
%      it looks like evidence.
%
% DUALS REQUIRE A RE-SOLVE. The certified artifact stores no multipliers, so
% this warm re-solves at the saved primal with returnModel/lamDef attached,
% mirroring run_foc_tulip. At a converged point that costs ~0 iterations.
%
% STATE LAYOUT differs from MEE and is easy to get wrong: x = [r(3); v(3); m; t]
% so MASS is row 7 and TIME is row 8, where the MEE state has mass at 6 and
% time at 7.
%
% INPUTS:
%   matPath - certified row (.mat) carrying out, sigma, rv0, rvf, tauf0, pSund
%             [char]. Default: the flagship.
%   opts    - struct (optional):
%             .maxIter  cap for the dual re-solve [default 800]
%             .relTol .absTol  integrator tolerances [default 1e-10, 1e-12]
%             .verbose  [default true]
%
% OUTPUTS:
%   R - struct: .Rx .Rlam [1xN] per-interval residuals; .dx [8xN] signed;
%       .defectCell [1xN]; .isSwitchCell; .hPhys; .rEarth (Earth distance at
%       interval midpoints, for reading residuals against orbit phase);
%       .Htau .HtauCoV (the conserved Hamiltonian and its variation);
%       .lamTimeCoV; .lamSign; .selfDefect .gatePass; .switchTimes; .nArcs;
%       .sigma .X .U .lam .lamDef
%
% REFERENCES:
%   [1] pmp_residual.m (the engine); pmp_residual_mee.m (the sibling adapter).
%   [2] GTO_tulip/direct/lib/casadi_minfuel_sundman.m (the dynamics mirrored).
%   [3] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, PHASE 3.

if nargin < 2, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
verbose = d('verbose', true);

here = fileparts(mfilename('fullpath'));
ot   = fileparts(fileparts(here));
addpath(here);
addpath(fullfile(ot,'verify_common','mesh'));
addpath(fullfile(ot,'cr3bp_common'));
addpath(fullfile(ot,'GTO_tulip','direct','lib'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
if nargin < 1 || isempty(matPath)
    matPath = fullfile(ot,'GTO_tulip','direct','lib','sundman_minfuel_certified.mat');
end
import casadi.*

S      = load(matPath);
sigma  = S.sigma(:).';
N1     = numel(sigma);
pSund  = 1.5;  if isfield(S,'pSund') && ~isempty(S.pSund), pSund = S.pSund; end
tauf0  = S.tauf0;
tf     = S.out.X(8,end);
p      = cr3bp_lt_params(0.025, 15, 2100);

% --- warm re-solve for the duals (the artifact stores none) -----------------
out = casadi_minfuel_sundman(sigma(:), tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
        S.out.X, S.out.U, tauf0, pSund, d('maxIter',800), 0, true, ...
        struct('returnModel', true));
assert(out.success && out.maxDefect <= 1e-8, 'pmp_residual_tulip:resolve', ...
    'dual re-solve failed (%s, defect %.2e)', out.ipoptStatus, out.maxDefect);
X = out.X;  U = out.U;

% --- rebuild the transcription's dynamics, MIRRORING the solver block -------
% Guards (1e-12) are part of the definition, not decoration: without them the
% two right-hand sides diverge at perigee, which is where the residual lives.
xs = MX.sym('x',8);  us = MX.sym('u',4);
r = xs(1:3);  v = xs(4:6);  m = xs(7);  al = us(1:3);  s = us(4);
dd = [r(1)+p.muStar; r(2); r(3)];
rr = [r(1)-1+p.muStar; r(2); r(3)];
r1 = sqrt(dd.'*dd + 1e-12);
d3 = (dd.'*dd + 1e-12)^1.5;  r3 = (rr.'*rr + 1e-12)^1.5;
gr = [r(1); r(2); 0] - (1-p.muStar)*dd/d3 - p.muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (s*p.Tmax/m)*al;
mdot  = -(p.Tmax/p.c)*s;
kappa = r1^pSund;
% THE tauf0 FACTOR IS NOT OPTIONAL. The collocation runs on sigma in [0,1],
% not on tau, and the solver's own defect line reads
%   D = X(:,2:end) - X(:,1:end-1) - tauf*(dsig/2).*(F(:,1:end-1)+F(:,2:end))
% i.e. dX/dsigma = tauf * dX/dtau. Omitting it (as a first version of this file
% did) understates the trapezoid term by a factor of tauf0 = 151.68, so the
% recomputed defect came out at 1.628 against the solver's 1.976e-14 -- caught
% by the self-gate, which refused the run rather than reporting the numbers
% that followed. The objective integrand carries the same factor.
fFun = Function('f', {xs,us}, {tauf0*kappa*[v; accel; mdot; 1]});
LFun = Function('L', {xs,us}, {tauf0*s*kappa});      % eps = 0 fuel integrand

% --- nodal costates from the interval defect duals --------------------------
lam0 = local_dual_to_costate(out.lamDef, sigma);

% --- switch structure --------------------------------------------------------
thr  = U(4,:);  burn = thr > 0.5;  tN = X(8,:);
sw   = mesh_switch_times(sigma, tN, burn, [], thr);
R.switchTimes  = sw.tSw;
R.isSwitchCell = false(1, N1-1);
if sw.n > 0, R.isSwitchCell(sw.bracket(:,1)) = true; end

% --- costate sign, measured on a probe arc rather than assumed --------------
probe = 1:min(40, N1);
rb = zeros(1,2);
for q = 1:2
    sgn = 3 - 2*q;
    Rp = pmp_residual(sigma(probe), X(:,probe), U(:,probe), sgn*lam0(:,probe), ...
        fFun, LFun, struct('variant','interp','verbose',false,'defectTol',Inf, ...
            'relTol', d('relTol',1e-10), 'absTol', d('absTol',1e-12)));
    rb(q) = median(Rp.Rlam_interp);
end
[~,qb] = min(rb);  lamSign = 3 - 2*qb;
lam = lamSign * lam0;
R.lamSign = lamSign;
if verbose
    fprintf('costate sign probe: +1 -> %.3e, -1 -> %.3e  => using %+d\n', rb(1), rb(2), lamSign);
end

% --- arc-by-arc residuals ----------------------------------------------------
edges = unique([1, find(diff(burn) ~= 0) + 1, N1]);
R.Rx = nan(1,N1-1);  R.Rlam = nan(1,N1-1);
R.defectCell = nan(1,N1-1);  R.dx = nan(8,N1-1);
selfD = 0;  Hall = nan(1,N1);
for a = 1:numel(edges)-1
    i1 = edges(a);  i2 = edges(a+1);
    if i2 - i1 < 1, continue; end
    idx = i1:i2;
    Ra = pmp_residual(sigma(idx), X(:,idx), U(:,idx), lam(:,idx), fFun, LFun, ...
        struct('variant','interp','verbose',false,'defectTol',Inf, ...
            'relTol', d('relTol',1e-10), 'absTol', d('absTol',1e-12)));
    R.Rx(i1:i2-1)   = Ra.Rx_interp;
    R.Rlam(i1:i2-1) = Ra.Rlam_interp;
    R.defectCell(i1:i2-1) = Ra.defectCell;
    R.dx(:,i1:i2-1) = Ra.dx;
    Hall(idx) = Ra.H;
    selfD = max(selfD, Ra.selfDefect);
end
R.selfDefect = selfD;  R.gatePass = selfD <= 1e-8;
R.nArcs = numel(edges)-1;

% --- conserved quantities ----------------------------------------------------
% H_tau is constant here because the regularized dynamics are AUTONOMOUS in
% tau. lam_t is constant because t enters only through its own equation.
R.Htau     = Hall;
R.HtauCoV  = std(Hall(isfinite(Hall))) / max(abs(mean(Hall(isfinite(Hall)))), realmin);
lt = lam(8,:);
R.lamTimeCoV = std(lt) / max(abs(mean(lt)), realmin);

lamMag = vecnorm(lam,2,1);
R.RlamRel = R.Rlam ./ max(0.5*(lamMag(1:end-1)+lamMag(2:end)), realmin);
R.hPhys   = diff(tN);
rE = vecnorm(X(1:3,:) - [-p.muStar;0;0], 2, 1);
R.rEarth  = 0.5*(rE(1:end-1) + rE(2:end));
R.sigma = sigma;  R.X = X;  R.U = U;  R.lam = lam;  R.lamDef = lamSign*out.lamDef;

if verbose
    [~,tag] = fileparts(matPath);
    fprintf('\n=== pmp_residual_tulip: %s ===\n', tag);
    fprintf('  nodes %d, arcs %d, switches %d\n', N1, R.nArcs, sw.n);
    fprintf('  self-gate defect : %.3e   %s   (solver reported %.3e)\n', ...
        R.selfDefect, local_tf(R.gatePass,'PASS','FAIL'), out.maxDefect);
    fprintf('  H_tau CoV        : %.3e   <- CONSTANT on an extremal (autonomous)\n', R.HtauCoV);
    fprintf('  lam_t CoV        : %.3e\n', R.lamTimeCoV);
    fprintf('  R_x   : median %.3e   max %.3e\n', median(R.Rx,'omitnan'), max(R.Rx));
    fprintf('  R_lam : median %.3e   max %.3e   (relative median %.3e)\n', ...
        median(R.Rlam,'omitnan'), max(R.Rlam), median(R.RlamRel,'omitnan'));
    fprintf('  ratio R_x / discrete defect : %.2e\n', ...
        median(R.Rx,'omitnan') / max(median(R.defectCell,'omitnan'), realmin));
end
end

% ---------------------------------------------------------------------------
function lam = local_dual_to_costate(LamDef, sigma)
% LOCAL_DUAL_TO_COSTATE  Interval defect duals -> nodal costates.
%
% The trapezoid defect of interval k couples nodes k and k+1, so LamDef(:,k) is
% a property of the INTERVAL, not of a node. The standard discrete-adjoint
% result is that the interior nodal costate is the STEP-WEIGHTED AVERAGE of its
% two adjacent interval duals, reducing to the plain average on a uniform mesh;
% the endpoints are one-sided. This mirrors mee_dual_to_costate, which the MEE
% campaign adopted after the same correctness argument.
%
% INPUTS:  LamDef [nx x N]; sigma [1x(N+1)]   OUTPUTS: lam [nx x (N+1)]
sigma = sigma(:).';
h  = diff(sigma);
nx = size(LamDef,1);  N = size(LamDef,2);
lam = zeros(nx, N+1);
lam(:,1)   = LamDef(:,1);
lam(:,N+1) = LamDef(:,N);
for k = 2:N
    w1 = h(k-1);  w2 = h(k);
    lam(:,k) = (w1*LamDef(:,k-1) + w2*LamDef(:,k)) / max(w1+w2, realmin);
end
end

% ---------------------------------------------------------------------------
function s = local_tf(c,a,b)
% LOCAL_TF  Pick a label by condition. INPUTS: c,a,b  OUTPUTS: s [char]
if c, s = a; else, s = b; end
end

% ---------------------------------------------------------------------------
function v = local_default(s,f,dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
