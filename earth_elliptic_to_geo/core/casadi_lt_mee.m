function out = casadi_lt_mee(sigma, X0, U0, dL0, opts)
% CASADI_LT_MEE  L-domain MEE-collocated 2-body low-thrust NLP (CasADi+IPOPT).
%
% Trapezoidal collocation in the true-longitude-fraction variable sigma, with
% the total true-longitude span DeltaL as a single scalar opti.variable (the
% L-domain analog of casadi_lt_2body's fixed-tauf/cScale-slack trick: here
% DeltaL multiplies every defect and the objective, exactly like tau_f would,
% but as ONE scalar column rather than a per-node slack state -- confirmed at
% the smoke-test wall-time budget to stay sparse, not dense-KKT).
%
% State per node x = [P; ex; ey; hx; hy; m; t] (7), control u = [beta(3); thr]
% (4) with beta a unit RTN thrust direction and thr in [0,1]. Node longitude
% L_k = pi + sigma_k*DeltaL (x0's L is pi, the paper's apogee start); dynamics
% at each node come from lt_mee_rhs (d/dL, MX-safe) with par.L set per node.
% Time is state row 7 (t(1)=0 pinned via the initial-state constraint); the
% terminal is expressed directly in elements as opts.xf (default GEO: P=1,
% ex=ey=hx=hy=0), leaving L free -- DeltaL is the transfer's free DOF.
%
% Modes:
%   'mintime' - thr == 1 (all-burn restriction) at every node; objective
%               min t(end); DeltaL (hence transfer time) found by the solve.
%   'fixedtf' - constraint t(end) = opts.tfTarget; Bertrand-Epenoy objective
%               J(eps) = Sum (dsig/2)*(w_k+w_{k+1}), w_k = (DeltaL/Ldot_k)*
%               (thr_k - eps*thr_k*(1-thr_k))  (dt = (DeltaL/Ldot)*dsigma).
%               eps=1 energy (smooth), eps=0 fuel (bang-bang).
%
% INPUTS:
%   sigma - uniform node parameter, 0->1 [(N+1)x1]
%   X0    - warm-start MEE states [P;ex;ey;hx;hy;m;t] [7x(N+1)]
%   U0    - warm-start controls [beta(3);thr] [4x(N+1)]
%   dL0   - warm-start total true-longitude span [scalar]
%   opts  - struct: .par (kepler_lt_params struct; .LdotMin default 1e-3 if
%           absent), .mode 'mintime'|'fixedtf', .eps (fixedtf homotopy param),
%           .tfTarget (fixedtf target transfer time), .x0 [7x1 initial MEE
%           state, t=0], .xf [5x1 terminal target [P;ex;ey;hx;hy], default
%           [1;0;0;0;0] = GEO], .maxIter, .warmTight, .printLevel, .selftest
%           (true -> skip the solve and return early with .xf resolved, for
%           option-parse unit tests -- Task 1), .returnModel (default false;
%           true -> ADDITIONALLY attach out.model with the live opti/X/U/dL
%           handles + a constraint registry creg and variable registry vreg,
%           for SOSC KKT assembly -- Task 2. Purely additive: with the flag
%           absent/false, X/U/dL are byte-identical and out.model is absent)
%
% OUTPUTS:
%   out - struct: .X [7x(N+1)] .U [4x(N+1)] .dL (converged DeltaL) .success
%         .ipoptStatus .maxDefect .maxUnit .termErr (distance to opts.xf)
%         .tfErr (|t(end)-tfTarget|, fixedtf only; NaN in mintime -- Task 7c
%         diagnostic) .mf .m_f_kg .dV_kms .tf .switches .edge .lamDef [7xN]
%         .LdotMin .incDeg (terminal inclination, deg -- should be ~0 for the
%         h=0 equatorial default target). .model (opts.returnModel only) -
%         struct('opti',opti,'X',X,'U',U,'dL',dL,'creg',creg,'vreg',vreg);
%         creg is a struct array, one per subject_to group, with fields
%         label[char] kind['eq'|'ineqLo'|'ineqHi'] rows[1xk] (row range in
%         opti.g, partitioning 1:size(opti.g,1) exactly once) bound[scalar|[]]
%         node[1xk|[]]; vreg indexes into opti.x (Xrows/Urows/nNode).
%
% A "casadi_lt_mee:boundSaturation" WARNING fires if t(end), P, or m sits
% within 1e-6 of its box bound at the solution -- maxDefect/termErr alone
% can both read machine precision while a box bound is quietly the real
% blocker (Task 7c, 2026-07-17: the old flat t<=300 ceiling did exactly
% this at the 1 N fuel rung).
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/casadi_lt_2body.m (Cartesian template this
%       mirrors: same IPOPT option regimes, out-struct reporting style).
%   [2] earth_elliptic_to_geo/lt_mee_rhs.m (L-domain Gauss dynamics, Task 1).
%   [3] earth_elliptic_to_geo/mee_seed.m (warm-start generator, Task 2).
%   [4] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004 (problem statement).
cp = getenv('CASADI_PATH');
if isempty(cp), cp = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cp);
par = opts.par;
if ~isfield(par, 'LdotMin') || isempty(par.LdotMin), par.LdotMin = 1e-3; end
if ~isfield(par, 'LdotFloor') || isempty(par.LdotFloor), par.LdotFloor = 1e-6; end
d = @(f,v) optdef(opts, f, v);
mode      = d('mode', 'fixedtf');
epsv      = d('eps', 0);
tfTarget  = d('tfTarget', []);
maxIter   = d('maxIter', 1500);
warmTight = d('warmTight', false);
printLvl  = d('printLevel', 0);
xf        = d('xf', [1;0;0;0;0]);
returnModel = d('returnModel', false);
assert(numel(xf)==5, 'casadi_lt_mee: opts.xf must be 5x1 [P;ex;ey;hx;hy]');
if d('selftest', false), out = struct('xf', xf(:)); return; end
assert(isfield(opts, 'x0') && ~isempty(opts.x0), ...
    'casadi_lt_mee requires opts.x0 (7x1 initial MEE state, t=0)');

N    = numel(sigma) - 1;
sg   = sigma(:);
dsig = diff(sg).';                                  % [1xN]

opti = casadi.Opti();
% SOSC registry (Task 2, additive only): records the row range each
% subject_to group occupies in opti.g, purely for later KKT assembly --
% zero effect on the solve itself (bracketing reads size(opti.g,1), never
% writes it).
creg = struct('label',{},'kind',{},'rows',{},'bound',{},'node',{});
addc = @(lab,kind,r0,bnd,nd) struct('label',lab,'kind',kind,'rows',r0:size(opti.g,1),'bound',bnd,'node',nd);
X  = opti.variable(7, N+1);
U  = opti.variable(4, N+1);
dL = opti.variable();                               % scalar DOF (single column)
m  = X(6,:);  t = X(7,:);  beta = U(1:3,:);  thr = U(4,:);
Lk = pi + sg*dL;                                     % [(N+1)x1] MX, x0's L = pi

% node dynamics dXdL_k = d/dL, Ldot_k = dL/dt (par.L set per node; Lk(k) is an
% MX expression depending on dL, which lt_mee_rhs handles fine since it only
% does arithmetic with par.L -- no norm/abs/max/if on state-dependent terms)
dXdLc = cell(1, N+1);  Ldotc = cell(1, N+1);
for k = 1:N+1
    parK = par;  parK.L = Lk(k);
    [dXdLc{k}, Ldotc{k}] = lt_mee_rhs(X(:,k), U(:,k), parK);
end
dXdL = [dXdLc{:}];      % [7x(N+1)] MX
Ldot = [Ldotc{:}];      % [1x(N+1)] MX

% collocation defects in sigma (KEEP HANDLES for the duals); dX/dsigma =
% DeltaL*dXdL, so DeltaL scales the defect exactly like tau_f would
r0 = size(opti.g,1)+1;
conDef = cell(1, N);
for k = 1:N
    conDef{k} = X(:,k+1) - X(:,k) - (dsig(k)/2)*dL*(dXdL(:,k) + dXdL(:,k+1)) == 0;
    opti.subject_to(conDef{k});
end
if returnModel, creg(end+1) = addc('defect','eq',r0,0,1:N); end

% Ldot degeneracy guard at every node
r0 = size(opti.g,1)+1;
for k = 1:N+1
    opti.subject_to(Ldot(k) >= par.LdotMin);
end
if returnModel, creg(end+1) = addc('ldotGuard','ineqLo',r0,par.LdotMin,1:N+1); end

% control cone + throttle (NEVER chain a<=x<=b -- MATLAB gotcha)
r0 = size(opti.g,1)+1;
for k = 1:N+1
    opti.subject_to(beta(1,k)^2 + beta(2,k)^2 + beta(3,k)^2 == 1);
end
if returnModel, creg(end+1) = addc('betaNorm','eq',r0,1,1:N+1); end
if strcmp(mode, 'mintime')
    r0 = size(opti.g,1)+1;
    opti.subject_to(thr == 1);
    if returnModel, creg(end+1) = addc('thrEq','eq',r0,1,1:N+1); end
else
    r0 = size(opti.g,1)+1;
    opti.subject_to(thr(:) >= 0);
    if returnModel, creg(end+1) = addc('thrLo','ineqLo',r0,0,1:N+1); end
    r0 = size(opti.g,1)+1;
    opti.subject_to(thr(:) <= 1);
    if returnModel, creg(end+1) = addc('thrHi','ineqHi',r0,1,1:N+1); end
end

% generous boxes (review lesson: bounds only block divergence); non-square/
% row-vector slices flattened with (:) per the Cartesian file's convention
r0 = size(opti.g,1)+1;  opti.subject_to(X(1,:).' >= 0.05);
if returnModel, creg(end+1) = addc('boxP_lo','ineqLo',r0,0.05,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(1,:).' <= 3);      % P
if returnModel, creg(end+1) = addc('boxP_hi','ineqHi',r0,3,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(2,:).' >= -1.5);
if returnModel, creg(end+1) = addc('boxEx_lo','ineqLo',r0,-1.5,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(2,:).' <= 1.5);    % ex
if returnModel, creg(end+1) = addc('boxEx_hi','ineqHi',r0,1.5,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(3,:).' >= -1.5);
if returnModel, creg(end+1) = addc('boxEy_lo','ineqLo',r0,-1.5,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(3,:).' <= 1.5);    % ey
if returnModel, creg(end+1) = addc('boxEy_hi','ineqHi',r0,1.5,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(4,:).' >= -2);
if returnModel, creg(end+1) = addc('boxHx_lo','ineqLo',r0,-2,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(4,:).' <= 2);      % hx
if returnModel, creg(end+1) = addc('boxHx_hi','ineqHi',r0,2,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(5,:).' >= -2);
if returnModel, creg(end+1) = addc('boxHy_lo','ineqLo',r0,-2,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(X(5,:).' <= 2);      % hy
if returnModel, creg(end+1) = addc('boxHy_hi','ineqHi',r0,2,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(m(:) >= 0.3);
if returnModel, creg(end+1) = addc('boxM_lo','ineqLo',r0,0.3,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(m(:) <= 1.001);
if returnModel, creg(end+1) = addc('boxM_hi','ineqHi',r0,1.001,1:N+1); end
% t upper bound: MUST scale with the problem's own timescale, not a fixed
% constant. Task 7c incident (2026-07-17): this was hardcoded to a flat 300
% (copied verbatim from casadi_lt_2body.m's "generous box"), which silently
% became a structural infeasibility at the 1 N fuel rung -- 'fixedtf' pins
% t(end)==tfTarget=335.71 ND (1.5x the certified 1 N min-time anchor), so
% t<=300 made the BC constraint and the box bound mutually unsatisfiable
% before IPOPT ever started iterating. Every resulting
% "Infeasible_Problem_Detected" showed t(end) parked at ~300.25 (the box
% ceiling) with IPOPT's own "Overall NLP error" reading ~35.46 (=tfTarget-
% 300.25) even while the re-derived maxDefect/termErr read machine
% precision -- a blind spot this file's diagnostics did not surface (see
% the bound-saturation WARNING below, added for the same reason). Fix:
% derive the ceiling from tfTarget (fixedtf) or the seed's own time span
% (mintime, where tfTarget is empty and t(end) is a free DOF found by the
% solve) so the box can never outrun the BC it is meant to merely bound.
if strcmp(mode, 'fixedtf')
    tUB = max(300, 2*tfTarget);
else
    tUB = max(300, 3*X0(7,end));   % mintime: 3x the warm-start seed's t-span
end
r0 = size(opti.g,1)+1;  opti.subject_to(t(:) >= 0);
if returnModel, creg(end+1) = addc('tBox_lo','ineqLo',r0,0,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(t(:) <= tUB);
if returnModel, creg(end+1) = addc('tBox_hi','ineqHi',r0,tUB,1:N+1); end
r0 = size(opti.g,1)+1;  opti.subject_to(beta(:) >= -1.01);
if returnModel, creg(end+1) = addc('betaBox_lo','ineqLo',r0,-1.01,[]); end
r0 = size(opti.g,1)+1;  opti.subject_to(beta(:) <= 1.01);
if returnModel, creg(end+1) = addc('betaBox_hi','ineqHi',r0,1.01,[]); end
r0 = size(opti.g,1)+1;  opti.subject_to(dL >= 0.1);
if returnModel, creg(end+1) = addc('dLbox_lo','ineqLo',r0,0.1,[]); end
% dL upper bound is RUNG-ADAPTIVE (external review, GPT-5.6, 2026-07-19): a fixed
% dL<=2000 (~318 rev) made the deep rungs structurally INFEASIBLE -- 0.2 N needs
% DeltaL~2168 (~345 rev) and 0.1 N ~4335 (~690 rev), both above 2000. Derive a
% generous ceiling from the warm-start C-law estimate dL0 (bounds "only block
% divergence"); the max(2000,..) floor leaves the shallow rungs unchanged (10 N
% dL0~46 -> 2000).
dLub = max(2000, 5*dL0);
r0 = size(opti.g,1)+1;  opti.subject_to(dL <= dLub);
if returnModel, creg(end+1) = addc('dLbox_hi','ineqHi',r0,dLub,[]); end

% boundary conditions
r0 = size(opti.g,1)+1;
opti.subject_to(X(:,1) == opts.x0(:));
if returnModel, creg(end+1) = addc('initBC','eq',r0,0,1); end
% terminal target in elements (default GEO [1;0;0;0;0]); L free (DeltaL is DOF).
% Prograde automatic for the h=0 equatorial default; a custom xf is the
% caller's responsibility (see run_gergaud scope note).
r0 = size(opti.g,1)+1;
for kt = 1:5
    opti.subject_to(X(kt,end) == xf(kt));
end
if returnModel, creg(end+1) = addc('termBC','eq',r0,0,N+1); end

% objective + t_f handling
if strcmp(mode, 'mintime')
    opti.minimize(t(end));
else
    assert(~isempty(tfTarget), 'fixedtf mode requires opts.tfTarget');
    r0 = size(opti.g,1)+1;
    opti.subject_to(t(end) == tfTarget);
    if returnModel, creg(end+1) = addc('tfPin','eq',r0,tfTarget,N+1); end
    w = (dL ./ fmax(Ldot, par.LdotFloor)) .* (thr - epsv*thr.*(1 - thr)); % dt/dsigma; Ldot guarded (see LdotFloor)
    opti.minimize(sum((dsig/2) .* (w(1:N) + w(2:N+1))));
end

% warm start + IPOPT. In 'mintime' mode thr is pinned ==1 at every node (a
% hard equality, not a box bound in CasADi's canonicalization -- confirmed by
% a minimal isolated Opti probe: opti.g() lists it as a general constraint,
% not lbx=ubx); seeding the warm start with the seed's own constant thr
% (e.g. mee_seed's thr=0.5) leaves a uniform ~0.5 primal infeasibility on
% that constraint from iteration 0, which this problem's conditioning does
% not shed quickly (empirically stuck >0.4 through 60 iterations). Priming
% the throttle warm start consistently with the pin removes that avoidable
% infeasibility outright (residual ~1e-13 at iter 0); state warm start X0
% is left untouched (it seeds the physical trajectory, mode-independent).
U0w = U0;
if strcmp(mode, 'mintime'), U0w(4,:) = 1; end
opti.set_initial(X, X0);
opti.set_initial(U, U0w);
opti.set_initial(dL, dL0);
ip = struct('max_iter', maxIter, 'tol', 1e-9, 'constr_viol_tol', 1e-10, ...
            'print_level', printLvl, 'mu_strategy', 'adaptive', ...
            'linear_solver', 'mumps', 'mumps_pivot_order', 0, ...
            'nlp_scaling_method', 'gradient-based');   % IPOPT default, set
            % explicitly and recorded on purpose (Campaign-B review lesson:
            % don't inherit a silent default); applies in both regimes below,
            % since warmTight only overrides the mu_/warm_start_ fields.
            % mumps_pivot_order=0 forces MUMPS' analysis-phase ordering
            % (ICNTL(7)) to AMD instead of the automatic choice, which at
            % N~193 nodes (Task 4's cross-formulation gate problem size)
            % selects METIS and hits a hard crash (abort(), not a catchable
            % MATLAB exception) inside this machine's bundled METIS build --
            % confirmed via crash-dump backtrace: MumpsSolverInterface::
            % SymbolicFactorization -> dmumps_ana_driver -> mumps_metis_
            % nodend_mixedto32 -> METIS_NodeND -> __CompressGraph/__GKfree.
            % Reproduced identically across two fresh MATLAB relaunches (not
            % the sporadic ~1/10 MEX-init crash) and independent of tfTarget
            % aggressiveness, so it is a genuine library bug for this
            % sparsity pattern, not a physics/transcription issue. AMD is
            % the standard robust fallback ordering and does not change any
            % problem formulation; verified this does not regress the N=75
            % Task-3 smoke test (both modes still pass).
if warmTight
    ip.mu_strategy = 'monotone';  ip.mu_init = 1e-4;
    ip.warm_start_init_point = 'yes';
    ip.warm_start_bound_push = 1e-9;  ip.warm_start_mult_bound_push = 1e-9;
end
opti.solver('ipopt', struct('print_time', printLvl > 0), ip);
success = true;
try
    sol = opti.solve();
catch
    sol = opti.debug;  success = false;
end
st = opti.stats();
status = st.return_status;
success = success && any(strcmp(status, {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));

% extraction + numeric re-check of the defects
Xs  = sol.value(X);   Us = sol.value(U);   dLs = sol.value(dL);
dmax = 0;  fn = zeros(7, N+1);  LdotN = zeros(1, N+1);
for k = 1:N+1
    parK = par;  parK.L = pi + sg(k)*dLs;
    [fn(:,k), LdotN(k)] = lt_mee_rhs(Xs(:,k), Us(:,k), parK);
end
for k = 1:N
    dk = Xs(:,k+1) - Xs(:,k) - (dsig(k)/2)*dLs*(fn(:,k) + fn(:,k+1));
    dmax = max(dmax, max(abs(dk)));
end
lamDef = nan(7, N);
try
    for k = 1:N, lamDef(:,k) = sol.value(opti.dual(conDef{k})); end
catch
end
ss = Us(4,:);
burn = ss > 0.5;
termErr = norm(Xs(1:5,end) - xf(:));
mf = Xs(6,end);
incDeg = 2*atand(sqrt(Xs(4,end)^2 + Xs(5,end)^2));
if strcmp(mode, 'fixedtf')
    tfErr = abs(Xs(7,end) - tfTarget);
else
    tfErr = nan;
end

% Bound-saturation diagnostic (Task 7c fix-round, 2026-07-17): maxDefect and
% termErr can BOTH read machine precision while a box bound quietly sits
% saturated at the solution -- that exact combination (t(end) pinned at the
% t<=300 ceiling, defect/termErr fine) is what cost a full task cycle on the
% 1 N fuel rung before the true root cause (the box, not the physics) was
% found. Cheap check, always run: flag t(end), P, or m within 1e-6 of
% either box edge. t(1)=0 and m(1)=1 are excluded (those are the pinned
% initial BCs, expected to sit at/near their bounds, not a symptom).
tolSat = 1e-6;
satMsgs = {};
if abs(Xs(7,end) - tUB) < tolSat || abs(Xs(7,end) - 0) < tolSat
    satMsgs{end+1} = sprintf('t(end)=%.6g saturates its box [0,%.6g]', Xs(7,end), tUB);
end
Pvals = Xs(1,:);
if any(abs(Pvals - 0.05) < tolSat) || any(abs(Pvals - 3) < tolSat)
    satMsgs{end+1} = 'P saturates its box [0.05,3] at >=1 node';
end
mvals = Xs(6,2:end);   % exclude the pinned m(1)=1 initial BC
if any(abs(mvals - 0.3) < tolSat) || any(abs(mvals - 1.001) < tolSat)
    satMsgs{end+1} = 'm saturates its box [0.3,1.001] at >=1 node';
end
if ~isempty(satMsgs)
    warning('casadi_lt_mee:boundSaturation', ...
        ['casadi_lt_mee: box-bound saturation at the reported solution ' ...
         '(maxDefect/termErr alone do NOT rule this out) -- %s'], ...
        strjoin(satMsgs, '; '));
end

out = struct('X', Xs, 'U', Us, 'dL', dLs, 'success', success, ...
    'ipoptStatus', status, 'maxDefect', dmax, ...
    'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), 'termErr', termErr, ...
    'tfErr', tfErr, ...
    'mf', mf, 'm_f_kg', par.m0kg*mf, 'dV_kms', par.c*log(1/mf)*par.VU_kms, ...
    'tf', Xs(7,end), 'switches', sum(abs(diff(burn))), ...
    'edge', mean(ss > 0.95 | ss < 0.05), 'lamDef', lamDef, ...
    'LdotMin', min(LdotN), 'incDeg', incDeg);

if returnModel
    vreg = struct('Xrows',1:7,'Urows',1:4,'nNode',N+1);
    out.model = struct('opti',opti,'X',X,'U',U,'dL',dL,'creg',creg,'vreg',vreg);
end
end
