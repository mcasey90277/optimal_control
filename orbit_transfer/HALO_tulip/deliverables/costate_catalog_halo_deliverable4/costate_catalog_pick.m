function [tf_nd, z8, info] = costate_catalog_pick(cat_, tauDRO, Np, ...
                                 depPhaseDays, arrPhaseDays, thrustN, warnFlag)
%% Purpose:
%
%   Looks up the multi-orbit costate CATALOG by five coordinates -- DRO
%   period, tulip petal count, departure phase, arrival phase, thrust --
%   and returns the minimum flight time (interpolated between the thrust
%   rungs that bracket the request) plus a converged z8 seed from the
%   nearest stored rung.
%
%   HONESTY CONTRACT (the warning flag): whenever what is RETURNED differs
%   from what was REQUESTED, a WARNING states exactly what you are getting:
%     - a different sheet (nearest DRO period and/or petal count),
%     - a snapped phase pair (requests between grid points),
%     - the nearest SOLVED pair (the requested cell is unsolved/red),
%     - a seed from a different thrust rung.
%   The info output always carries requested-vs-delivered, warnings on or
%   off.
%
%  ASSUMPTIONS / NOTES:
%
% • Phases are days past each orbit's reference point (the family getters),
%   modulo the period. The DRO period is the sheet's tauDRO (ND).
% • For thrusts between rungs, fly the returned z8's own rung or hand it to
%   tfMin at your thrust and let it converge the small correction.
% • pm = +1 tulips: z-mirror the pm = -1 solution (see cat_.derive).
%
%% Inputs:
%
%  cat_                     struct                  A compact costate catalog
%                                                   (any family: dro_tulip,
%                                                   halo_tulip, ...)
%
%  tauDRO                   double                  Requested DRO period (ND)
%
%  Np                       double                  Requested petal count
%
%  depPhaseDays             double                  Departure phase (days)
%
%  arrPhaseDays             double                  Arrival phase (days)
%
%  thrustN                  double                  Requested thrust (N)
%
%  warnFlag                 logical                 Print warnings when the
%                                                   delivery differs from
%                                                   the request (default
%                                                   true)
%
%% Outputs:
%
%  tf_nd                    double                  Minimum flight time (ND)
%                                                   at the requested thrust,
%                                                   linear between rungs
%
%  z8                       [8 x 1]                 [lambda0; tf] seed from
%                                                   the nearest stored rung
%                                                   (tfMin-ready)
%
%  info                     struct                  .requested / .delivered
%                                                   (tauDRO, Np, depDays,
%                                                   arrDays, thrustN each),
%                                                   .warned [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: an off-grid request, warnings on. Runs on whichever compact
   %catalog sits in the current folder:
     F = dir('costate_catalog_*.mat');
     assert(~isempty(F), ...
            'demo: no costate_catalog_*.mat in the current folder');
       L = load(F(1).name);
      fn = fieldnames(L);
    cat_ = L.(fn{1});
   fprintf('demo catalog: %s\n', F(1).name);
      sh = cat_.sheets(1);
[tf,z,inf_] = costate_catalog_pick(cat_, sh.tauDRO*0.93, sh.Np+1, ...
                                   3.3, 11.0, 4.0);
   fprintf('tf = %.4f ND;  delivered sheet tau=%.2f Np=%d\n', ...
           tf, inf_.delivered.tauDRO, inf_.delivered.Np);
   return;
end
if ~exist('warnFlag','var'), warnFlag = true; end
    tStar = cat_.constants.tStar_s;
   warned = false;

%% Nearest sheet (log-distance in period; petal count integer-nearest):
     taus = [cat_.sheets.tauDRO];
      nps = [cat_.sheets.Np];
     dist = (log([cat_.sheets.tauDRO]) - log(tauDRO)).^2 ...
          + 0.5*(nps - Np).^2;                   % petal mismatch weighted
   [~, ks] = min(dist);
       sh = cat_.sheets(ks);
if warnFlag && (abs(sh.tauDRO - tauDRO) > 1e-9 || sh.Np ~= Np)
    fprintf(['WARNING (costate_catalog_pick): no sheet at DRO period %.3f ND / %d petals;\n', ...
             '  you are getting the NEAREST SHEET: period %.2f ND (%.2f d), Np = %d.\n'], ...
        tauDRO, Np, sh.tauDRO, sh.tauDRO*tStar/86400, sh.Np);
    warned = true;
end

%% Phase pair on that sheet's torus:
       Pd = sh.tauDRO * tStar/86400;             % DRO period (days)
       Pa = sh.period_tulip_nd * tStar/86400;    % tulip period (days)
       fd = mod(depPhaseDays/Pd, 1);
       fa = mod(arrPhaseDays/Pa, 1);
       dd = abs(sh.sD_frac - fd);  dd = min(dd, 1-dd);
       da = abs(sh.sA_frac - fa);  da = min(da, 1-da);
   [~,iD] = min(dd);
   [~,iA] = min(da);
if warnFlag && (dd(iD)*Pd > 1e-6 || da(iA)*Pa > 1e-6)
    fprintf(['WARNING (costate_catalog_pick): requested phases (dep %.4f d, arr %.4f d) are\n', ...
             '  OFF-GRID on this sheet; returning the nearest grid pair: dep %.4f d, arr %.4f d.\n'], ...
        depPhaseDays, arrPhaseDays, sh.sD_frac(iD)*Pd, sh.sA_frac(iA)*Pa);
    warned = true;
end

%% Unsolved cell? Deliver the nearest SOLVED pair, loudly:
      hav = squeeze(sh.has_solution(iD,iA,:))';
if ~any(hav)
    [okD, okA] = find(any(sh.has_solution,3));
    d2 = min(abs(sh.sD_frac(okD)'-fd),1-abs(sh.sD_frac(okD)'-fd)).^2 ...
       + min(abs(sh.sA_frac(okA)'-fa),1-abs(sh.sA_frac(okA)'-fa)).^2;
    [~, kb] = min(d2);
    iD = okD(kb);  iA = okA(kb);
    hav = squeeze(sh.has_solution(iD,iA,:))';
    if warnFlag
        fprintf(['WARNING (costate_catalog_pick): the requested phase pair is UNSOLVED on this\n', ...
                 '  sheet. You are getting the nearest SOLVED pair: dep %.4f d, arr %.4f d.\n'], ...
            sh.sD_frac(iD)*Pd, sh.sA_frac(iA)*Pa);
        warned = true;
    end
end

%% Flight time: linear between bracketing rungs; seed: nearest stored rung:
     rung = cat_.rungs_N;
   rAvail = rung(hav);
       tf = squeeze(sh.tf_nd(iD,iA,:))';
  tfAvail = tf(hav);
if thrustN <= min(rAvail)
    tf_nd = tfAvail(rAvail == min(rAvail));
elseif thrustN >= max(rAvail)
    tf_nd = tfAvail(rAvail == max(rAvail));
else
       lo = max(rAvail(rAvail <= thrustN));
       hi = min(rAvail(rAvail >= thrustN));
    if hi == lo, tf_nd = tfAvail(rAvail == lo);
    else
        tf_nd = tfAvail(rAvail==lo) + ...
            (tfAvail(rAvail==hi) - tfAvail(rAvail==lo))*(thrustN-lo)/(hi-lo);
    end
end
   [~,kn] = min(abs(rAvail - thrustN));
    kNear = find(rung == rAvail(kn), 1);
       z8 = sh.z8(:, sh.entry_index(iD,iA,kNear));
if warnFlag && abs(rung(kNear) - thrustN) > 1e-9
    fprintf(['WARNING (costate_catalog_pick): no stored rung at %.3f N here; the z8 seed is\n', ...
             '  the %.2f N entry. The interpolated tf above refers to your requested thrust.\n'], ...
        thrustN, rung(kNear));
    warned = true;
end

info = struct( ...
    'requested', struct('tauDRO',tauDRO,'Np',Np,'depDays',depPhaseDays, ...
                        'arrDays',arrPhaseDays,'thrustN',thrustN), ...
    'delivered', struct('tauDRO',sh.tauDRO,'Np',sh.Np, ...
                        'depDays',sh.sD_frac(iD)*Pd,'arrDays',sh.sA_frac(iA)*Pa, ...
                        'tfThrustN',min(max(thrustN,min(rAvail)),max(rAvail)), ...
                        'seedThrustN',rung(kNear)), ...
    'warned', warned);
end
