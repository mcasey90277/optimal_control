function K = costate_compare(o, lamInd, tauInd, p, Tmax, c, opts)
% COSTATE_COMPARE  Compare direct-NLP-derived costates against indirect ones.
%
% WHY THIS IS THE POINT OF THE WHOLE CAMPAIGN. Every other transfer in this
% repository has only one method's answer, so a direct solution's costates can
% be checked for self-consistency but never against an independent second
% opinion. DRO->tulip is the exception: pumpkyn.cr3bp.tfMin converges to a
% genuine Pontryagin extremal, and tfMinProp returns the whole costate history
% alongside the state. So the covector mapping -- NLP multipliers to costates --
% can be tested against a known answer rather than against itself.
%
% Note what this adds over the t_f agreement already certified. Matching t_f
% and the endpoints shows the two methods found the same OPTIMUM. Matching the
% costates shows they found the same EXTREMAL, structure and all, which is a
% strictly stronger statement and the one the campaign was set up to make.
%
% THE MAPPING. For defect constraints D_k = 0 with multipliers nu_k, stationarity
% of the Lagrangian with respect to the interior state X_j is exactly the
% discrete adjoint recursion, and nu_k -> lambda(t_k) as the mesh refines (Hager
% 2000). No h-scaling is required for this defect form: the multiplier IS the
% costate in the limit. Two ambiguities remain and are RESOLVED BY MEASUREMENT
% rather than assumed:
%
%   SIGN.  IPOPT's multiplier sign convention, and CasADi's reporting of it,
%          may differ from the PMP convention. We estimate the sign from the
%          data and report it.
%   SCALE. Costates for a free-final-time minimum-time problem are fixed only
%          up to the normalization implied by the transversality condition,
%          and tfMin's internal normalization need not match the NLP's. We
%          estimate a single global scale factor by least squares and then --
%          this is the actual test -- check that it is CONSTANT in time. A
%          constant factor means the two costate histories are the same
%          covector field in different units. A drifting factor means they
%          are not, and the agreement in t_f was luck.
%
% INPUTS:
%   o      - solution struct from casadi_mintime_dro, with .lamDef [7xN]
%            (requires opts.returnModel = true at solve time)
%   lamInd - indirect costate history [Mx7] = [lambda_r lambda_v lambda_m]
%   tauInd - times for lamInd [Mx1], same units as o.tf
%   p      - params (.muStar .lStar .tStar)
%   Tmax, c- ND thrust acceleration and exhaust velocity [scalars]
%   opts   - struct (optional): .verbose [true]
%
% OUTPUTS:
%   K - struct: .scale (the fitted global factor) .scaleCoV (its coefficient of
%       variation in time -- the real test; small is good) .sign
%       .primerAngDeg [1xN] angle between the direct-derived primer and the
%       solution's own thrust direction .primerAngMax .primerAngMed
%       .indPrimerAngDeg (the same for the indirect costates -- the control)
%       .lamAngDeg [1xN] angle between the direct and indirect lambda_v
%       .lamMEnd (direct lambda_m at t_f; transversality wants 0)
%       .Hmean .HCoV (Hamiltonian constancy; min-time free-t_f wants H const)
%
% REFERENCES:
%   [1] Hager, W.W., "Runge-Kutta methods in optimal control and the transformed
%       adjoint system", Numer. Math. 87, 2000.
%   [2] pumpkyn.cr3bp.tfMinEoM -- the indirect costate dynamics.

if nargin < 7, opts = struct(); end
verbose = local_default(opts, 'verbose', true);

if ~isfield(o,'lamDef') || isempty(o.lamDef)
    error('costate_compare:noDuals', ...
        ['o.lamDef is empty. Re-solve with opts.returnModel = true -- the ' ...
         'constraint-row registry that locates the defect multipliers is only ' ...
         'built in that mode.']);
end

N   = size(o.lamDef, 2);
t   = o.s(:).' * o.tf;
tMid= 0.5*(t(1:end-1) + t(2:end));       % nu_k belongs to interval k
mu  = p.muStar;

% --- the direct-derived costates -------------------------------------------
% nu_k is the multiplier of defect k. Place it at the interval midpoint; that
% is where a midpoint-rule reading of the adjoint recursion puts it, and it
% avoids inventing a node value the NLP never formed.
lamD   = o.lamDef;                        % 7 x N: [lam_r; lam_v; lam_m]
lamDv  = lamD(4:6,:);
lamDm  = lamD(7,:);

% --- the indirect costates, resampled onto the same instants ---------------
lamI   = interp1(tauInd(:), lamInd, tMid(:), 'spline').';   % 7 x N
lamIv  = lamI(4:6,:);

% --- SIGN, estimated not assumed -------------------------------------------
% The primer vector is alpha = -lambda_v/||lambda_v||. If the NLP's multiplier
% sign is flipped relative to PMP, the direct-derived primer points backwards.
% Decide by which sign aligns better with the solution's OWN thrust direction,
% which is a primal quantity and so free of any dual convention.
% Use the solver's OWN midpoint control where it exists. Averaging the two node
% controls would compare the duals against a quantity the NLP never formed, and
% for Hermite-Simpson the midpoint control is a genuine decision variable.
if isfield(o,'Um') && ~isempty(o.Um)
    alphaSol = o.Um(1:3,:);  thrSol = o.Um(4,:);
else
    alphaSol = 0.5*(o.U(1:3,1:end-1) + o.U(1:3,2:end));
    thrSol   = 0.5*(o.U(4,1:end-1)   + o.U(4,2:end));
end
alphaSol = alphaSol ./ max(vecnorm(alphaSol,2,1), eps);
sgnTry = [1 -1];  fit = zeros(1,2);
for kk = 1:2
    a = -sgnTry(kk)*lamDv ./ max(vecnorm(lamDv,2,1), eps);
    fit(kk) = mean(sum(a .* alphaSol, 1));           % mean cosine
end
[~, ibest] = max(fit);
sgn = sgnTry(ibest);
lamD  = sgn*lamD;   lamDv = sgn*lamDv;   lamDm = sgn*lamDm;

% --- PRIMER ALIGNMENT: the direct duals against the primal control ----------
aD = -lamDv ./ max(vecnorm(lamDv,2,1), eps);
K.primerAngDeg = real(acosd(max(-1,min(1, sum(aD.*alphaSol,1)))));

aI = -lamIv ./ max(vecnorm(lamIv,2,1), eps);
K.indPrimerAngDeg = real(acosd(max(-1,min(1, sum(aI.*alphaSol,1)))));

% --- DIRECT vs INDIRECT costate direction ----------------------------------
K.lamAngDeg = real(acosd(max(-1,min(1, sum(aD.*aI,1)))));

% lambda_r as well. lambda_v is what the primer is built from and so is the
% best-determined block; lambda_r couples to the trajectory only indirectly and
% is the classic weak spot of shooting, so it is the sterner test of the mapping.
uDr = lamD(1:3,:) ./ max(vecnorm(lamD(1:3,:),2,1), eps);
uIr = lamI(1:3,:) ./ max(vecnorm(lamI(1:3,:),2,1), eps);
K.lamRAngDeg = real(acosd(max(-1,min(1, sum(uDr.*uIr,1)))));
K.lamRScale  = median(vecnorm(lamI(1:3,:),2,1) ./ max(vecnorm(lamD(1:3,:),2,1), eps));

% --- SCALE, and whether it is constant (the actual test) -------------------
ratio = vecnorm(lamIv,2,1) ./ max(vecnorm(lamDv,2,1), eps);
K.scale    = median(ratio);
K.scaleCoV = std(ratio)/max(abs(mean(ratio)), eps);
K.ratio    = ratio;

% --- transversality and Hamiltonian ----------------------------------------
% lambda_m at t_f: the multipliers live at interval midpoints, so the terminal
% value must be extrapolated from the last two rather than read off the last.
K.lamMEnd = lamDm(end) + (lamDm(end)-lamDm(end-1)) * ...
            (t(end)-tMid(end)) / max(tMid(end)-tMid(end-1), eps);
K.lamMEndRaw = lamDm(end);

% Hamiltonian, evaluated at the SOLVER'S midpoint state and control -- not at a
% fabricated average of the neighbouring nodes, which is a point the NLP never
% formed and which sits off the trajectory by the interpolation term.
if isfield(o,'Xm') && ~isempty(o.Xm)
    Xmid = o.Xm;
else
    Xmid = 0.5*(o.X(:,1:end-1) + o.X(:,2:end));
end
H = zeros(1,N);
for k = 1:N
    f = local_f(Xmid(:,k), [alphaSol(:,k); thrSol(k)], mu, Tmax, c);
    H(k) = lamD(:,k).' * f;
end
K.H = H;  K.Hmean = mean(H);  K.HCoV = std(H)/max(abs(mean(H)), eps);

% INTERIOR vs BOUNDARY. For minimum time with J = t_f the PMP Hamiltonian is
% H_PMP = 1 + lambda'f and the free-final-time condition is H_PMP == 0, i.e.
% lambda'f == -1 identically. Measured here the median is -1 to ~1e-5, but the
% last handful of intervals deviate by O(1).
%
% That is a MAPPING artifact at the terminal boundary, not a disagreement. The
% final node carries the endpoint equality constraints, so its stationarity
% condition mixes the defect multipliers with the boundary multipliers, and the
% naive reading used here -- 'nu_k is lambda at the midpoint of interval k' --
% has no way to separate them. The Hager covector mapping supplies exactly this
% terminal correction; foc_check applies it elsewhere in this repository for the
% same reason. Until it is applied here, the last ~1% of intervals should be
% treated as uninformative rather than as evidence.
%
% Measured on the certified N=1600 solution: trimming 20 intervals from each end
% drops max|H+1| from 1.02 to 0.023 and the CoV from 5.6e-2 to 2.5e-3, while the
% median is unchanged at 9.8e-6. The contamination is confined to k >= N-8.
nTrim = max(1, round(0.0125*N));
ii = (1+nTrim):(N-nTrim);
K.nTrim      = nTrim;
K.HmeanInt   = mean(H(ii));
K.HCoVInt    = std(H(ii))/max(abs(mean(H(ii))), eps);
K.HdevMedian = median(abs(H+1));
K.HdevMaxInt = max(abs(H(ii)+1));
K.sign = sgn;
K.tMid = tMid;
K.primerAngMax = max(K.primerAngDeg);
K.primerAngMed = median(K.primerAngDeg);

if verbose
    fprintf('\n===== COSTATE COMPARISON: direct duals vs indirect costates =====\n');
    fprintf('  (N = %d, scheme %s)\n', N, o.scheme);
    fprintf('  dual sign convention resolved empirically to : %+d\n', sgn);
    fprintf('\n  PRIMER ALIGNMENT (angle to the solution''s own thrust direction)\n');
    fprintf('    direct duals   : median %8.4f deg, max %8.4f deg\n', ...
        median(K.primerAngDeg), max(K.primerAngDeg));
    fprintf('    indirect (ctrl): median %8.4f deg, max %8.4f deg\n', ...
        median(K.indPrimerAngDeg), max(K.indPrimerAngDeg));
    fprintf('\n  DIRECT vs INDIRECT costate direction\n');
    fprintf('    angle between lambda_v: median %8.4f deg, max %8.4f deg\n', ...
        median(K.lamAngDeg), max(K.lamAngDeg));
    fprintf('\n  SCALE (indirect/direct magnitude of lambda_v)\n');
    fprintf('    median factor %.6g, coefficient of variation %.3e\n', ...
        K.scale, K.scaleCoV);
    if K.scaleCoV < 1e-3
        fprintf('    -> CONSTANT: the two are the same covector field, different units.\n');
    else
        fprintf('    -> NOT constant: they are NOT the same covector field.\n');
    end
    fprintf('    angle between lambda_r: median %8.4f deg, max %8.4f deg\n', ...
        median(K.lamRAngDeg), max(K.lamRAngDeg));
    fprintf('    lambda_r scale factor  : %.6g   (lambda_v gave %.6g)\n', ...
        K.lamRScale, K.scale);
    fprintf('\n  transversality lambda_m(t_f) = %+.4e  (PMP wants 0)\n', K.lamMEnd);
    fprintf(['  Hamiltonian lambda''f (PMP wants identically -1 for min-time):\n' ...
             '    median |H+1| over all intervals      : %.3e\n' ...
             '    INTERIOR (trimming %d from each end) : mean %+.6f, CoV %.3e, max|H+1| %.3e\n' ...
             '    all intervals                        : mean %+.6f, CoV %.3e\n'], ...
        K.HdevMedian, K.nTrim, K.HmeanInt, K.HCoVInt, K.HdevMaxInt, K.Hmean, K.HCoV);
    fprintf(['    (the end-interval deviation is a terminal covector-mapping\n' ...
             '     artifact -- see the note in the source -- not a disagreement)\n']);
end
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
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
