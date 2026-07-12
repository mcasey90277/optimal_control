function res = ifs_tf_arclength(fac0, facEnd, ds, opts)
% IFS_TF_ARCLENGTH  Pseudo-arclength continuation of the IFS k=0 all-burn branch
% across the min-time fold, up in transfer time to the first switch birth.
%
% Naive t_f-stepping cannot cross the min-time point: d(terminal)/d(lambda0) is
% singular there (vertical tangent of the branch). Pseudo-arclength (Keller)
% treats the transfer-time factor as an UNKNOWN and parameterizes the branch by
% arclength: the extended Jacobian [dR/dlambda0 | dR/dfactor] (8x9) stays full
% rank through the fold, so the corrector is nonsingular. State x =
% [lambda0(8); factor(1)]; residual F = [R(x); tangent.(x-x_prev) - ds]; solved
% by a predictor (tangent step) + Newton corrector in a scaled metric.
%
% The Sundman length tauf is held fixed within each corrector (so the t=tf
% residual is meaningful and the R-block is square 8x8), and updated between
% accepted steps to the physically consistent value.
%
% INPUTS:
%   fac0   - starting factor (anchor), e.g. 1.00 [scalar]
%   facEnd - factor to march toward, e.g. 1.15 [scalar]
%   ds     - arclength step in the scaled metric, e.g. 0.05 [scalar]
%   opts   - (optional) struct: maxSteps [200], corrTol [1e-8], corrIter [20],
%            odeOpts, verbose [true]
% OUTPUTS:
%   res - struct array per accepted point: factor, tf, tf_days, tauf, resNorm,
%         maxS, minS, dsUsed, lam0
%
% STATUS (2026-07-11): the predictor/corrector machinery is correct (correctors
% converge; the extended 8x9 Jacobian is full rank, so the tangent is a clean
% 1-D null vector). BUT the min-time anchor is pathologically degenerate: the
% all-burn branch is effectively VERTICAL in (t_f, lambda0) and carries a
% near-null costate GAUGE (smallest scaled singular value ~8e-4). The loose
% corrector (forced loose by the ~1e-6 min-time seed accuracy) drifts along that
% gauge at fixed t_f, so a march never advances `factor` and `max S` wanders
% non-monotonically -- the reported switch "birth" at factor=1.0000 is a gauge
% artifact, NOT physical. Next lever: regularize the gauge (phase/pinning
% condition on the near-null lambda0 direction) or anchor above the fold. Full
% record: RESULTS_RUNG01_RUNG2.md (section Rung 2b).
%
% REFERENCES:
%   [1] Keller, "Numerical solution of bifurcation and nonlinear eigenvalue
%       problems," 1977 (pseudo-arclength).
%   [2] ifs/RESULTS_RUNG01_RUNG2.md (the fold finding this crosses).

here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
if nargin < 4, opts = struct(); end
if ~isfield(opts,'maxSteps'),   opts.maxSteps = 200;  end
if ~isfield(opts,'corrTol'),    opts.corrTol = 1e-7;  end
if ~isfield(opts,'corrAcceptR'),opts.corrAcceptR = 1e-5; end  % accept this R-floor near the fold
if ~isfield(opts,'corrIter'),   opts.corrIter = 25;   end
if ~isfield(opts,'dsMax'),      opts.dsMax = 0.6;     end
if ~isfield(opts,'verbose'),  opts.verbose = true;  end
if ~isfield(opts,'odeOpts'),  opts.odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-14); end
tfMin = 6.290694;  tStarDays = 382981.289129055/86400;

% ---- anchor -----------------------------------------------------------
[Za, proba] = ifs_seed_mintime(fac0, opts.odeOpts);
outa = ifs_solve2(Za, proba, struct('tolR',1e-8,'maxIter',80,'verbose',false));
lam0 = outa.Z;  tauf = proba.tauf;
x = [lam0(:); fac0];
xscale = [max(abs(lam0(:)), 1e-2); 1];         % arclength metric scale
proto = proba;

res = mkrec(x, tauf, norm(Rk0(x, proto, tfMin, tauf)), proto, tStarDays, ds);
if opts.verbose
    fprintf('%7s %9s %10s %11s %8s %10s\n','factor','tf(day)','||R||','maxS','corrOK','ds');
    printrec(res(end));
end

% initial tangent: prefer increasing factor
t = branchTangent(x, proto, tfMin, tauf, xscale, [zeros(8,1);1]);

step = 0;
while x(9) < facEnd - 1e-9 && step < opts.maxSteps
    step = step + 1;
    xs      = x ./ xscale;
    xs_pred = xs + ds * t;                     % predictor (scaled space)

    [xs_new, ok, rn, facTry] = corrector(xs_pred, xs, t, ds, xscale, proto, tfMin, tauf, opts);
    if ~ok
        if opts.verbose
            fprintf('  [corr fail] step %d ds=%.4f -> facTry=%.5f rn=%.2e; halving ds\n', ...
                    step, ds, facTry, rn);
        end
        ds = ds/2;                             % corrector failed: shorten and retry
        if ds < 1e-4
            fprintf('  [stop] arclength step underflow near factor=%.4f\n', x(9));
            break;
        end
        continue;
    end

    xnew = xscale .* xs_new;
    % update tangent (oriented), tauf (physical), and record
    told = t;
    t = branchTangent(xnew, proto, tfMin, tauf, xscale, told);
    taufNew = updateTauf(xnew, proto, tfMin, opts.odeOpts, tauf);
    r = mkrec(xnew, taufNew, rn, proto, tStarDays, ds);
    res(end+1) = r; %#ok<AGROW>
    if opts.verbose, printrec(r); end

    x = xnew;  tauf = taufNew;
    if rn < opts.corrAcceptR, ds = min(ds*1.5, opts.dsMax); end   % grow on easy steps

    if r.maxS > 0
        fprintf(['\n*** FIRST SWITCH BIRTH: max S crossed 0 at factor ~ %.4f ' ...
                 '(t_f=%.3f d, maxS=%.3e). The all-burn (k=0) branch ends here ' ...
                 '-- lower edge of the transition band from the indirect side. ***\n'], ...
                r.factor, r.tf_days, r.maxS);
        break;
    end
end
save(fullfile(here,'tf_arclength_results.mat'),'res');
end

% ======================================================================
function R = Rk0(x, proto, tfMin, tauf)
% k=0 rendezvous residual [8x1] at state x=[lam0;factor] with fixed tauf.
lam0 = x(1:8);  factor = x(9);
prob = proto;  prob.k = 0;  prob.uArc = 1;  prob.tauParam = 'direct';
prob.tf = factor*tfMin;  prob.tauf = tauf;
R = ifs_residual(lam0, prob);
end

% ----------------------------------------------------------------------
function J = dRdxs(xs, xscale, proto, tfMin, tauf)
% Complex-step Jacobian dR/dxs (8x9) in scaled coordinates xs (x = xscale.*xs).
h = 1e-30;  n = numel(xs);  J = zeros(8, n);
for j = 1:n
    xsp = xs;  xsp(j) = xsp(j) + 1i*h;
    Rp = Rk0(xscale.*xsp, proto, tfMin, tauf);
    J(:, j) = imag(Rp)/h;
end
end

% ----------------------------------------------------------------------
function t = branchTangent(x, proto, tfMin, tauf, xscale, prefer)
% Unit tangent (scaled space) = null vector of dR/dxs, oriented so t.prefer>0.
xs = x ./ xscale;
Js = dRdxs(xs, xscale, proto, tfMin, tauf);     % 8x9 (full rank 8)
[~,~,V] = svd(Js);                              % FULL svd: V is 9x9
t = V(:, end);                                  % 9th col = null vector = branch tangent
if t.'*prefer(:) < 0, t = -t; end
t = t / norm(t);
end

% ----------------------------------------------------------------------
function [xs, ok, rn, facTry] = corrector(xs, xs_prev, t, ds, xscale, proto, tfMin, tauf, opts)
% Newton corrector on G(xs) = [R(xscale.*xs); t.(xs-xs_prev) - ds] = 0 (9x9).
% Near the min-time fold the R-block cannot beat its own accuracy floor (the
% min-time seed is ~1e-6-accurate), so acceptance is on ||R|| reaching
% corrAcceptR with the arclength constraint satisfied -- not on ||G||<corrTol.
ok = false;  rn = NaN;
for it = 1:opts.corrIter
    R  = Rk0(xscale.*xs, proto, tfMin, tauf);
    Na = t.'*(xs - xs_prev) - ds;
    G  = [R; Na];  gn = norm(G);  rn = norm(R);
    if gn < opts.corrTol, ok = true;  facTry = xscale(9)*xs(9);  return; end
    Js = dRdxs(xs, xscale, proto, tfMin, tauf);  % 8x9
    Jf = [Js; t.'];                              % 9x9
    dxs = -(Jf \ G);
    a = 1;  accepted = false;
    for ls = 1:24
        xt = xs + a*dxs;
        Rt = Rk0(xscale.*xt, proto, tfMin, tauf);
        Gt = [Rt; t.'*(xt - xs_prev) - ds];
        if norm(Gt) < (1 - 1e-4*a)*gn, accepted = true; break; end
        a = a/2;
    end
    if ~accepted
        % stalled: accept if R hit its achievable floor and arclength holds
        ok = (rn < opts.corrAcceptR) && (abs(Na) < 1e-3);
        facTry = xscale(9)*xs(9);
        return;
    end
    xs = xs + a*dxs;
end
R  = Rk0(xscale.*xs, proto, tfMin, tauf);  rn = norm(R);
Na = t.'*(xs - xs_prev) - ds;
ok = (rn < opts.corrAcceptR) && (abs(Na) < 1e-3);
facTry = xscale(9)*xs(9);
end

% ----------------------------------------------------------------------
function tauf = updateTauf(x, proto, tfMin, odeOpts, taufOld)
% Physical Sundman length: integrate hard-burn from lam0 to t=tf.
lam0 = x(1:8);  tf = x(9)*tfMin;
y0 = [proto.rv0(:);1;0;lam0(:)];
ev = odeset(odeOpts,'Events',@(s,y) tfEv(s,y,tf));
sol = ode113(@(s,y) ifs_eom(s,y,proto.Tmax,proto.c,proto.muStar,proto.pSund,1),[0 400],y0,ev);
if isempty(sol.xe), tauf = taufOld; else, tauf = sol.xe(end); end
end

% ----------------------------------------------------------------------
function [mx,mn] = maxSarc(lam0, proto, tauf)
[~,Y] = ode113(@(s,y) ifs_eom(s,y,proto.Tmax,proto.c,proto.muStar,proto.pSund,1), ...
               [0 tauf], [proto.rv0(:);1;0;lam0(:)], proto.odeOpts);
S = 1 - sqrt(sum(Y(:,12:14).^2,2)).*proto.c./Y(:,7) - Y(:,15);
mx = max(S);  mn = min(S);
end

function r = mkrec(x, tauf, rn, proto, tStarDays, ds)
[mx,mn] = maxSarc(x(1:8), proto, tauf);
r = struct('factor',x(9),'tf',x(9)*6.290694,'tf_days',x(9)*6.290694*tStarDays, ...
    'tauf',tauf,'resNorm',rn,'maxS',mx,'minS',mn,'dsUsed',ds,'lam0',x(1:8));
end

function printrec(r)
fprintf('%7.4f %9.3f %10.3e %11.4e %8s %10.4f\n', ...
        r.factor, r.tf_days, r.resNorm, r.maxS, 'ok', r.dsUsed);
end

function [v,it,d] = tfEv(~,y,tf)
v = y(8)-tf; it = 1; d = 1;
end
