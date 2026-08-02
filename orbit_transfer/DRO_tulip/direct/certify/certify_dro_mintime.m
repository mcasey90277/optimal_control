function C = certify_dro_mintime(o, p, Tmax, c, opts)
% CERTIFY_DRO_MINTIME  Pass/fail gate for a direct DRO->tulip min-time solution.
%
% WHY THIS EXISTS. Every solve in this campaign reported a defect of 1e-14 and
% every one of them was wrong. The defect is a statement about the
% discretization; it is not evidence that the trajectory is physical. This gate
% is the standard a solve must meet before its t_f is quoted as an answer.
%
% THE TWO GATES THAT MATTER, and both are independent of the solver:
%
%   G1  ACCURACY. The worst continuous-time local error must be below a
%       position tolerance stated in KILOMETRES, not in solver units. Default
%       1 km. Measured by re-integrating each interval at 1e-11/1e-13 -- see
%       dro_residual.
%
%   G2  AGREEMENT. This is the one problem in the catalog where an independently
%       converged INDIRECT solution exists (pumpkyn.cr3bp.tfMin, t_f = 4.015243).
%       A working direct solver must reproduce it. Default tolerance 1e-4
%       relative, i.e. agreement to four digits.
%
% G2 is the sharper test and it is the reason DRO->tulip was chosen to start
% with. A direct method that cannot reproduce a known answer on the one problem
% where the answer is known should not be trusted on the problems where it is
% not.
%
% THE SUPPORTING CHECKS (G3-G6) are the usual transcription hygiene: defect,
% control unit norm, terminal boundary error, throttle saturation, the lifted-
% time spread, and -- when a floor was imposed -- whether the trajectory honours
% it BETWEEN nodes and not merely at them. They are necessary, not sufficient,
% and they are reported separately so that a solve which passes all of them and
% still fails G1 is visible for what it is.
%
% INPUTS:
%   o     - solution struct from casadi_mintime_dro
%   p     - params from dro_tulip_endpoints (.muStar .lStar .tStar)
%   Tmax  - ND thrust acceleration at unit mass fraction [scalar]
%   c     - ND exhaust velocity [scalar]
%   opts  - struct (optional):
%           .posTolKm  [1.0]      G1 tolerance, km
%           .tfRef     [4.0152430] G2 reference t_f, ND. [] skips G2.
%           .tfRelTol  [1e-4]     G2 relative tolerance
%           .defectTol [1e-9]     G3
%           .rMoonKm   [1737.4]
%           .verbose   [true]
%
% OUTPUTS:
%   C - struct: .pass (G1 && G2), .passAll (every gate), .gates (struct array
%       with .id .name .value .tol .pass .units), .resid (from dro_residual),
%       .worstAltKm, .tfErrRel
%
% REFERENCES:
%   [1] orbit_transfer/DRO_tulip/FINDINGS.md
%   [2] orbit_transfer/OPTIMALITY_CERTIFICATION.md -- the repo-wide register
%       this gate reports into.

if nargin < 5, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
posTolKm  = d('posTolKm', 1.0);
tfRef     = d('tfRef', 4.0152430);
tfRelTol  = d('tfRelTol', 1e-4);
defectTol = d('defectTol', 1e-9);
rMoonKm   = d('rMoonKm', 1737.4);
verbose   = d('verbose', true);

lStar = p.lStar;  mu = p.muStar;

%% the measurement
R = dro_residual(o, mu, Tmax, c);
worstKm = R.RxMax * lStar;

r2 = vecnorm(o.X(1:3,:) - [1-mu;0;0], 2, 1);
altKm = r2*lStar - rMoonKm;
worstAltKm = altKm(min(R.kWorst+1, numel(altKm)));

%% the gates
G = struct('id',{},'name',{},'value',{},'tol',{},'pass',{},'units',{});
add = @(id,nm,val,tol,ps,un) struct('id',id,'name',nm,'value',val,'tol',tol,'pass',ps,'units',un);

G(end+1) = add('G1','continuous accuracy (worst interval)', worstKm, posTolKm, ...
               worstKm <= posTolKm, 'km');

if isempty(tfRef)
    tfErrRel = NaN;
else
    tfErrRel = abs(o.tf - tfRef)/tfRef;
    G(end+1) = add('G2','agreement with the indirect reference t_f', tfErrRel, ...
                   tfRelTol, tfErrRel <= tfRelTol, 'rel');
end

G(end+1) = add('G3','NLP defect', o.maxDefect, defectTol, ...
               o.maxDefect <= defectTol, '-');
G(end+1) = add('G4','control unit-norm error', o.maxUnit, 1e-8, ...
               o.maxUnit <= 1e-8, '-');
G(end+1) = add('G5','terminal boundary error', o.termErr, 1e-9, ...
               o.termErr <= 1e-9, '-');
G(end+1) = add('G6','throttle saturation (min u; min-time expects 1)', ...
               o.thrMin, 1-1e-6, o.thrMin >= 1-1e-6, '-');

% the floor, if one was imposed, checked BETWEEN nodes as well as at them
if isfield(o,'minAltKm') && ~isnan(o.minAltKm)
    trueMinKm = local_true_min_alt(o, mu, Tmax, c, lStar, rMoonKm);
    G(end+1) = add('G7','altitude floor honoured between nodes', trueMinKm, ...
                   o.minAltKm, trueMinKm >= o.minAltKm, 'km');
end

pass    = all([G(strcmp({G.id},'G1')).pass]);
if ~isempty(tfRef), pass = pass && G(strcmp({G.id},'G2')).pass; end
passAll = all([G.pass]);

C = struct('pass',pass, 'passAll',passAll, 'gates',G, 'resid',R, ...
           'worstAltKm',worstAltKm, 'tfErrRel',tfErrRel, 'worstKm',worstKm);

%% report
if verbose
    fprintf('\n===== CERTIFICATION: DRO->tulip min-time (%s, N = %d) =====\n', ...
        R.scheme, numel(o.s)-1);
    fprintf('%-4s %-48s %12s %12s  %s\n','','gate','value','tolerance','');
    for k = 1:numel(G)
        fprintf('%-4s %-48s %12.4e %12.4e  %s\n', G(k).id, G(k).name, ...
            G(k).value, G(k).tol, local_mark(G(k).pass));
    end
    fprintf('  worst interval sits at %.0f km lunar altitude\n', worstAltKm);
    fprintf('  t_f = %.6f ND = %.3f days\n', o.tf, o.tf*p.tStar/86400);
    if pass
        fprintf('  VERDICT: CERTIFIED (G1 accuracy and G2 agreement both pass)\n');
    else
        fprintf('  VERDICT: NOT CERTIFIED -- this t_f is not an answer\n');
    end
end
end

% ---------------------------------------------------------------------------
function amin = local_true_min_alt(o, mu, Tmax, c, lStar, rMoonKm)
% LOCAL_TRUE_MIN_ALT  Minimum lunar altitude of the PROPAGATED trajectory, not
% of the nodes. A collocation floor binds at nodes only; this is what checks it.
% INPUTS: o; mu; Tmax; c; lStar; rMoonKm   OUTPUTS: amin [km]
t = o.s(:).' * o.tf;
odeo = odeset('RelTol',1e-11,'AbsTol',1e-13);
amin = inf;
hasMid = isfield(o,'Um') && ~isempty(o.Um);
for k = 1:numel(t)-1
    a = t(k);  b = t(k+1);  h = max(b-a,eps);
    ua = o.U(:,k);  ub = o.U(:,k+1);
    if hasMid
        um = o.Um(:,k);
        uOf = @(tt) local_q(ua, um, ub, (tt-a)/h);
    else
        uOf = @(tt) ua + ((tt-a)/h)*(ub-ua);
    end
    [~,Z] = ode113(@(tt,z) local_f(z, uOf(tt), mu, Tmax, c), linspace(a,b,9), o.X(:,k), odeo);
    rr = vecnorm(Z(:,1:3) - [1-mu 0 0], 2, 2);
    amin = min(amin, min(rr)*lStar - rMoonKm);
end
end

% ---------------------------------------------------------------------------
function u = local_q(ua, um, ub, w)
% LOCAL_Q  Lagrange quadratic through ua, um, ub at w = 0, 1/2, 1.
% INPUTS: ua, um, ub [4x1]; w   OUTPUTS: u [4x1]
u = (2*w-1).*(w-1).*ua + 4*w.*(1-w).*um + w.*(2*w-1).*ub;
end

% ---------------------------------------------------------------------------
function dz = local_f(z, u, mu, Tmax, c)
% LOCAL_F  CR3BP with thrust; mirrors dro_residual/local_rhs.
% INPUTS: z [7x1]; u [4x1]; mu; Tmax; c   OUTPUTS: dz [7x1]
r = z(1:3);  v = z(4:6);  m = z(7);
al = u(1:3);  al = al/max(norm(al),eps);  th = min(max(u(4),0),1);
dd = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
gr = [r(1) - (1-mu)*(r(1)+mu)/dd^3 - mu*(r(1)-1+mu)/rr^3;
      r(2) - (1-mu)*r(2)/dd^3      - mu*r(2)/rr^3;
           - (1-mu)*r(3)/dd^3      - mu*r(3)/rr^3];
dz = [v; gr + [2*v(2); -2*v(1); 0] + (th*Tmax/m)*al; -(Tmax/c)*th];
end

% ---------------------------------------------------------------------------
function s = local_mark(tf)
% LOCAL_MARK  'PASS'/'FAIL' text for the report.
% INPUTS: tf (logical)   OUTPUTS: s [char]
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
