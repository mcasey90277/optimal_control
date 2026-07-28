function R = pmp_residual_mee(matPath, opts)
% PMP_RESIDUAL_MEE  Continuous PMP residuals for a certified earth MEE row.
%
% Campaign adapter for pmp_residual: rebuilds the L-domain dynamics as CasADi
% Functions from the transcription's OWN right-hand side, recovers nodal
% costates from the NLP duals, splits the trajectory at its bang-bang switches,
% and evaluates the two-point-BVP residuals arc by arc.
%
% WHY ARC BY ARC AND NOT ON ONE GRID. At a bang-bang switch the control is
% genuinely double-valued -- thr = 1 as the burn arc closes, thr = 0 as the
% coast arc opens -- so a single shared node cannot carry both, and an interval
% spanning the switch interpolates between them. The toy validation measured
% exactly this: a shared node produced a spurious 1.77e-02 residual that the
% self-gate correctly refused. Splitting at switches makes the control
% continuous on each piece; the intervals that STRADDLE a switch are reported
% separately, because their residual is the physically interesting quantity
% rather than an artifact.
%
% THE DYNAMICS ARE THE TRANSCRIPTION'S, NOT A RE-DERIVATION. lt_mee_rhs is
% called directly and dH/dx is taken from it by AD. The check is independent in
% INTEGRATION only; re-deriving the Gauss equations here would test the
% author's algebra rather than the mesh.
%
% NON-AUTONOMOUS IN SIGMA. The node longitude is L = pi + sigma*DeltaL, so the
% Gauss equations depend explicitly on the independent variable and the
% Hamiltonian in sigma is NOT expected to be constant. The conserved quantity
% for this transcription is the TIME COSTATE: t enters the dynamics only
% through its own trivial equation, so lam_t is constant along an extremal and
% equals -H_t. That, not H_sigma, is the Hamiltonian-type diagnostic here, and
% it is reported as lamTimeCoV.
%
% COSTATE SIGN. The adjoint equation is NOT invariant under lam -> -lam,
% because the running-cost term dL/dx does not flip while (df/dx)'*lam does.
% Both signs are therefore evaluated and the one with the smaller costate
% residual is reported, together with BOTH residuals, so the choice is visible
% rather than silent. This is a measurement of which convention the duals carry,
% not a fit.
%
% INPUTS:
%   matPath - certified row (.mat), e.g. results/MEE_M2_10N.mat [char]
%   opts    - struct (optional):
%             .maxArcs   cap on arcs evaluated, for a quick look [default Inf]
%             .relTol .absTol  integrator tolerances [default 1e-10, 1e-12]
%             .verbose   [default true]
%
% OUTPUTS:
%   R - struct:
%       .Rx .Rlam        [1xN] per-interval residuals over the whole trajectory
%                        (NaN on intervals that straddle a switch and were not
%                        evaluated within an arc)
%       .isSwitchCell    [1xN] logical, intervals containing a switch
%       .hPhys           [1xN] physical time step per interval
%       .sigma .X .U .lam  the solution the residuals belong to
%       .lamSign         +1 or -1, the convention selected
%       .lamTimeCoV      coefficient of variation of the time costate
%       .selfDefect .gatePass   the self-validating gate
%       .switchTimes     physical switch times
%       .RlamRel         [1xN] costate residual RELATIVE to local |lam|. The
%                        raw .Rlam is an absolute norm and costates are
%                        sensitivities that can be enormous, so only this
%                        normalized form is interpretable.
%       .rOrb .rCell     orbital radius at nodes / interval midpoints, so the
%                        residual can be read against orbit PHASE
%       .revs .nodesPerRev .switchesPerRev
%       .nArcs
%
% REFERENCES:
%   [1] pmp_residual.m (the engine); mesh_switch_times.m (sub-grid switches).
%   [2] earth_elliptic_to_geo/direct/core/lt_mee_rhs.m (the dynamics reused).
%   [3] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, PHASE 3.

if nargin < 2, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
verbose = d('verbose', true);

here = fileparts(mfilename('fullpath'));
ot   = fileparts(fileparts(here));                    % .../orbit_transfer
addpath(here);
addpath(fullfile(ot, 'verify_common', 'mesh'));
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
modRoot = fileparts(fileparts(matPath));
cwd0 = pwd;  cleaner = onCleanup(@() cd(cwd0)); %#ok<NASGU>
cd(modRoot);  setup_paths();
addpath(fullfile(ot, 'verify_common'));  setup_verify_common();

import casadi.*

% --- refresh the row to obtain corrected duals alongside the primal ----------
[out, par, sigma] = refresh_duals_mee(matPath, struct('returnModel', false));
X = out.X;  U = out.U;  dL = out.dL;
sigma = sigma(:).';
N1 = numel(sigma);
lam0 = mee_dual_to_costate(out.lamDef, sigma(:));      % [7 x N1]

% --- the transcription's dynamics, as sigma-parameterized CasADi Functions ---
% dX/dsigma = DeltaL * dX/dL, with par.L set from sigma exactly as
% casadi_lt_mee does (L_k = pi + sigma_k*DeltaL).
xs = MX.sym('x',7);  us = MX.sym('u',4);  ss = MX.sym('s',1);
parS = par;  parS.L = pi + ss*dL;
[dXdL, Ldot] = lt_mee_rhs(xs, us, parS);
fFun = Function('f', {xs,us,ss}, {dL * dXdL});
% eps = 0 fuel objective: dJ/dsigma = (DeltaL/Ldot)*thr
LFun = Function('L', {xs,us,ss}, {(dL / Ldot) * us(4)});

% --- switch structure --------------------------------------------------------
thr  = U(4,:);
burn = thr > 0.5;
tN   = X(7,:);
sw   = mesh_switch_times(sigma, tN, burn, [], thr);
R.switchTimes  = sw.tSw;
R.isSwitchCell = false(1, N1-1);
if sw.n > 0, R.isSwitchCell(sw.bracket(:,1)) = true; end

% --- costate sign: measured, not assumed ------------------------------------
% Evaluate a short probe arc under both conventions and keep the better.
probe = 1:min(20, N1);
rBoth = zeros(1,2);
for q = 1:2
    sgn = 3 - 2*q;                                    % +1 then -1
    Rp = pmp_residual(sigma(probe), X(:,probe), U(:,probe), sgn*lam0(:,probe), ...
        fFun, LFun, struct('variant','interp','verbose',false,'defectTol',Inf, ...
                           'relTol', d('relTol',1e-10), 'absTol', d('absTol',1e-12)));
    rBoth(q) = median(Rp.Rlam_interp);
end
[~, qBest] = min(rBoth);
lamSign = 3 - 2*qBest;
lam = lamSign * lam0;
R.lamSign = lamSign;
if verbose
    fprintf('costate sign probe: +1 -> %.3e, -1 -> %.3e  => using %+d\n', ...
            rBoth(1), rBoth(2), lamSign);
end

% --- time-costate constancy (this transcription's conserved quantity) -------
lt = lam(7,:);
R.lamTimeCoV = std(lt) / max(abs(mean(lt)), realmin);

% --- arc-by-arc residuals ----------------------------------------------------
edges = [1, find(diff(burn) ~= 0) + 1, N1];
edges = unique(edges);
nArcs = numel(edges) - 1;
if isfinite(d('maxArcs', Inf)), nArcs = min(nArcs, d('maxArcs', Inf)); end
R.Rx   = nan(1, N1-1);
R.Rlam = nan(1, N1-1);
R.defectCell = nan(1, N1-1);
selfD  = 0;
for a = 1:nArcs
    i1 = edges(a);  i2 = edges(a+1);
    if i2 - i1 < 1, continue; end
    idx = i1:i2;
    Ra = pmp_residual(sigma(idx), X(:,idx), U(:,idx), lam(:,idx), fFun, LFun, ...
        struct('variant','interp','verbose',false,'defectTol',Inf, ...
               'relTol', d('relTol',1e-10), 'absTol', d('absTol',1e-12)));
    R.Rx(i1:i2-1)   = Ra.Rx_interp;
    R.Rlam(i1:i2-1) = Ra.Rlam_interp;
    R.defectCell(i1:i2-1) = Ra.defectCell;
    selfD = max(selfD, Ra.selfDefect);
end
% --- normalization and orbital context --------------------------------------
% R_lam is an ABSOLUTE norm and costates are sensitivities, which can be
% enormous -- a residual of 13.9 against a costate of 1e5 is relatively tiny.
% Interpreting the raw number is therefore meaningless; normalize by the local
% costate magnitude (mean of the interval's two endpoints).
lamMag = vecnorm(lam, 2, 1);
lamCell = 0.5*(lamMag(1:end-1) + lamMag(2:end));
R.lamMag  = lamMag;
R.RlamRel = R.Rlam ./ max(lamCell, realmin);

% Orbital radius from the MEE state, so residuals can be read against orbit
% PHASE rather than only against time. In MEE with p = X(1),
%   r = p / (1 + ex*cos(L) + ey*sin(L)),   L = pi + sigma*DeltaL.
Lnode = pi + sigma * dL;
R.rOrb = X(1,:) ./ (1 + X(2,:).*cos(Lnode) + X(3,:).*sin(Lnode));
R.rCell = 0.5*(R.rOrb(1:end-1) + R.rOrb(2:end));
R.revs        = dL / (2*pi);
R.nodesPerRev = (N1-1) / R.revs;
R.switchesPerRev = sw.n / R.revs;

R.selfDefect = selfD;
R.gatePass   = selfD <= 1e-8;
R.nArcs = nArcs;
R.sigma = sigma;  R.X = X;  R.U = U;  R.lam = lam;
R.hPhys = diff(tN);

if verbose
    fprintf('\n=== pmp_residual_mee: %s ===\n', matPath);
    fprintf('  nodes %d, arcs %d, switches %d\n', N1, nArcs, sw.n);
    fprintf('  self-gate defect  : %.3e   %s   (transcription reported %.3e)\n', ...
        R.selfDefect, local_tf(R.gatePass,'PASS','FAIL'), out.maxDefect);
    fprintf('  lam_t CoV         : %.3e   (constant on an extremal)\n', R.lamTimeCoV);
    fprintf('  R_x   over arcs   : median %.3e   max %.3e\n', ...
        median(R.Rx,'omitnan'), max(R.Rx));
    fprintf('  R_lam over arcs   : median %.3e   max %.3e\n', ...
        median(R.Rlam,'omitnan'), max(R.Rlam));
    fprintf('  ratio R_x / nodal defect : %.2e   <- how much bigger the CONTINUOUS\n', ...
        median(R.Rx,'omitnan') / max(out.maxDefect, realmin));
    fprintf('                              residual is than the DISCRETE one\n');
end
end

% ---------------------------------------------------------------------------
function s = local_tf(c, a, b)
% LOCAL_TF  Pick a label by condition. INPUTS: c,a,b  OUTPUTS: s [char]
if c, s = a; else, s = b; end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
