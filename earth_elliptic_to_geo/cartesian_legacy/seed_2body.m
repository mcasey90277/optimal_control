function [sigma, X0, U0, tauf0, info] = seed_2body(par, rv0, opts)
% SEED_2BODY  Dynamically-exact tangential-thrust warm start on uniform-tau nodes.
%
% Propagates s = sbar, alpha = vhat from rv0 (ode113, tight tol), computes the
% Sundman clock tau(t) = int dt / r^pSund, and samples the DENSE ode solution at
% N+1 uniform-tau nodes (deval => defect-free to ODE tolerance; no downsampling
% interpolation — the campaign's no-resample lesson). If opts.targetLf is given,
% bisects sbar so the arrival unwrapped equatorial longitude lands within pi/2
% of it (winding is monotone-decreasing in sbar: more thrust climbs sooner ->
% slower angular rate -> less longitude wound).
%
% INPUTS:
%   par  - kepler_lt_params struct
%   rv0  - initial inertial state [6x1]
%   opts - .sbar [scalar], .tDur [scalar ND, or [] => stop at GEO energy -mu/2],
%          .N [segments], .targetLf [optional, rad unwrapped from L0]
%
% OUTPUTS:
%   sigma - [(N+1)x1] uniform 0->1;  X0 [9x(N+1)] = [r;v;m;t;cScale=1];
%   U0    - [4x(N+1)] = [vhat; sbar];  tauf0 - total tau length;
%   info  - .sbar .Larr (arrival unwrapped longitude) .tEnd .mEnd
%
% REFERENCES: [1] DESIGN.md secs 3-4. [2] sundman_minfuel/sundman_seed_map.m.
sbar = opts.sbar;  N = opts.N;
if isfield(opts,'targetLf') && ~isempty(opts.targetLf)
    lo = 0.25;  hi = 1.0;                     % Larr(lo) > Larr(hi)
    for kb = 1:14
        mid  = 0.5*(lo+hi);
        Lm   = propagate(par, rv0, mid, opts.tDur).Larr;
        if abs(Lm - opts.targetLf) < pi/2, sbar = mid; break; end
        if Lm > opts.targetLf, lo = mid; else, hi = mid; end
        sbar = 0.5*(lo+hi);
    end
end
S = propagate(par, rv0, sbar, opts.tDur);
% Sundman map on a fine grid, then uniform-tau node times
tt  = linspace(0, S.tEnd, 20*N+1);
xx  = deval(S.sol, tt);
rr  = sqrt(sum(xx(1:3,:).^2, 1));
tau = cumtrapz(tt, rr.^(-par.pSund));         % dtau = dt / r^p  (cScale=1)
tauf0 = tau(end);
tN  = interp1(tau, tt, linspace(0, tauf0, N+1));
XN  = deval(S.sol, tN);                       % 8 x (N+1), exact to ODE tol
vN  = XN(4:6,:);  vn = max(sqrt(sum(vN.^2,1)), 1e-9);
sigma = linspace(0, 1, N+1).';
X0  = [XN; ones(1, N+1)];
U0  = [vN ./ vn; sbar*ones(1, N+1)];
info = struct('sbar', sbar, 'Larr', S.Larr, 'tEnd', S.tEnd, 'mEnd', XN(7,end));
end

% ---------------------------------------------------------------------------
function S = propagate(par, rv0, sbar, tDur)
% Tangential constant-throttle propagation; empty tDur => stop at GEO energy.
odef = @(t,x) lt2b_rhs_time(x, [x(4:6)/max(norm(x(4:6)),1e-9); sbar], par);
oo = odeset('RelTol',1e-11, 'AbsTol',1e-12);
if isempty(tDur)
    oo  = odeset(oo, 'Events', @geoEnergyEvent);
    sol = ode113(odef, [0 500], [rv0(:); 1; 0], oo);
else
    sol = ode113(odef, [0 tDur], [rv0(:); 1; 0], oo);
end
tf_ = sol.x(end);
tq  = linspace(0, tf_, 4000);
xq  = deval(sol, tq);
Lun = unwrap(atan2(xq(2,:), xq(1,:)));
S = struct('sol', sol, 'tEnd', tf_, 'Larr', Lun(end) + (pi - Lun(1)));
end

function [val, isterm, dir_] = geoEnergyEvent(~, x)
% Stop when two-body energy reaches the GEO value -mu/(2a), a=1, mu=1.
val = 0.5*(x(4)^2+x(5)^2+x(6)^2) - 1/norm(x(1:3)) - (-0.5);
isterm = 1;  dir_ = 1;
end
