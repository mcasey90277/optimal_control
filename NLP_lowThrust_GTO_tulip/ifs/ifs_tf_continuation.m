function res = ifs_tf_continuation(fac0, facEnd, dfac, opts)
% IFS_TF_CONTINUATION  March the IFS min-fuel solution up in transfer time from
% the min-time all-burn anchor, tracking the switching function toward the first
% switch birth.
%
% Rung 2 (t_f-continuation). Anchored at the k=0 all-burn min-time solution
% (IFS_SEED_MINTIME), each step raises the transfer-time factor by dfac and
% re-solves the k=0 rendezvous shooting problem warm-started from the previous
% converged initial costate. The min-fuel optimum stays all-burn (k=0) while the
% switching function S<0 everywhere; the maximum of S over the arc rises toward
% 0 as t_f grows. When max S crosses 0 a coast becomes beneficial -- the FIRST
% SWITCH BIRTH -- which locates the lower edge of the campaign's open 1.01-1.11x
% transition band from the indirect side. (Switch INSERTION past the birth is a
% follow-on; this driver detects and localizes it.)
%
% INPUTS:
%   fac0   - starting t_f factor (anchor), e.g. 1.00 [scalar]
%   facEnd - final t_f factor to attempt, e.g. 1.15 [scalar]
%   dfac   - factor step, e.g. 0.01 [scalar]
%   opts   - (optional) struct: tolR [1e-8], maxIter [80], warmTol accept
%            threshold for a usable warm step [1e-5], odeOpts
% OUTPUTS:
%   res - struct array per factor: factor, tf, tf_days, tauf, resNorm, success,
%         maxS (max switching fn over the arc), minS, lam0
%
% STATUS: naive parameter-stepping in t_f does NOT cross the min-time fold
% (d(terminal)/d(lambda0) is singular there): even a 0.1% step fails to converge.
% See RESULTS_RUNG01_RUNG2.md. The fold-aware sibling is IFS_TF_ARCLENGTH (which
% crosses the corrector issue but exposes a deeper min-time costate-gauge
% degeneracy). This driver is kept as the diagnostic that localized the fold.
%
% REFERENCES: ifs/PLAN_OF_ATTACK.md (Rung 2); ms_band/sms_seed_mintime.m.

here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
if nargin < 4, opts = struct(); end
if ~isfield(opts,'tolR'),    opts.tolR = 1e-8;    end
if ~isfield(opts,'maxIter'), opts.maxIter = 80;   end
if ~isfield(opts,'warmTol'), opts.warmTol = 1e-5; end
if ~isfield(opts,'odeOpts'), opts.odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-14); end
solveOpts = struct('tolR',opts.tolR,'maxIter',opts.maxIter,'verbose',false);
tStarDays = 382981.289129055/86400;

% --- anchor ------------------------------------------------------------
[Z, prob, meta] = ifs_seed_mintime(fac0, opts.odeOpts);
out = ifs_solve2(Z, prob, solveOpts);
lam0 = out.Z;
res = mkrec(fac0, prob, out, lam0, tStarDays);
fprintf('%7s %9s %10s %11s %11s %8s\n','factor','tf(day)','||R||','maxS','minS','ok');
printrec(res(end));

fac = fac0;
while fac < facEnd - 1e-9
    facNew = min(fac + dfac, facEnd);
    probN  = buildK0Prob(facNew, lam0, prob, opts.odeOpts);
    outN   = ifs_solve2(lam0, probN, solveOpts);
    r      = mkrec(facNew, probN, outN, outN.Z, tStarDays);
    res(end+1) = r; %#ok<AGROW>
    printrec(r);

    if outN.resNorm > opts.warmTol
        fprintf('  [warn] step did not converge tight (||R||=%.2e > %.0e); continuation seed degrading.\n', ...
                outN.resNorm, opts.warmTol);
    end
    lam0 = outN.Z;  prob = probN;  fac = facNew;

    if r.maxS > 0
        fprintf(['\n*** FIRST SWITCH BIRTH: max S crossed 0 at factor ~ %.4f ' ...
                 '(t_f=%.3f d, maxS=%.3e). A coast becomes optimal here -- the ' ...
                 'all-burn (k=0) branch ends; this is the lower edge of the ' ...
                 'transition band from the indirect side. ***\n'], ...
                r.factor, r.tf_days, r.maxS);
        break;
    end
end
save(fullfile(here,'tf_continuation_results.mat'),'res');
end

% ======================================================================
function prob = buildK0Prob(factor, lam0, protoProb, odeOpts)
% k=0 problem at a new factor: same endpoints/params, new tf, and tauf found by
% integrating the hard-burn EOM from the warm lam0 to the event t=tf.
tfMin = 6.290694;
tf = factor*tfMin;
y0 = [protoProb.rv0(:); 1; 0; lam0(:)];
ev = odeset(odeOpts,'Events',@(s,y) tf_event(y, tf));
sol = ode113(@(s,y) ifs_eom(s,y,protoProb.Tmax,protoProb.c,protoProb.muStar, ...
                            protoProb.pSund,1), [0 400], y0, ev);
assert(~isempty(sol.xe),'buildK0Prob: t never reached tf (factor=%.3f)',factor);
prob = protoProb;  prob.tf = tf;  prob.tauf = sol.xe(end);
end

% ----------------------------------------------------------------------
function [mx, mn] = maxSarc(lam0, prob)
[~,Y] = ode113(@(s,y) ifs_eom(s,y,prob.Tmax,prob.c,prob.muStar,prob.pSund,1), ...
               [0 prob.tauf], [prob.rv0(:);1;0;lam0(:)], prob.odeOpts);
S = 1 - sqrt(sum(Y(:,12:14).^2,2)).*prob.c./Y(:,7) - Y(:,15);
mx = max(S);  mn = min(S);
end

% ----------------------------------------------------------------------
function r = mkrec(factor, prob, out, lam0, tStarDays)
[mx,mn] = maxSarc(lam0, prob);
r = struct('factor',factor,'tf',prob.tf,'tf_days',prob.tf*tStarDays, ...
    'tauf',prob.tauf,'resNorm',out.resNorm,'success',out.success, ...
    'maxS',mx,'minS',mn,'lam0',lam0(:));
end

function printrec(r)
fprintf('%7.4f %9.3f %10.3e %11.4e %11.4e %8d\n', ...
        r.factor, r.tf_days, r.resNorm, r.maxS, r.minS, r.success);
end

function [val,isterm,dir] = tf_event(y, tf)
val = y(8) - tf;  isterm = 1;  dir = 1;
end
