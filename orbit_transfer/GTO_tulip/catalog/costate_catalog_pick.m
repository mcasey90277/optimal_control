function [tf_nd, z8, info] = costate_catalog_pick(cat_, tauDep, arrKey, ...
                                 depPhaseDays, arrPhaseDays, thrustN, warnFlag)
%% Purpose:
%
%   Looks up the multi-orbit costate CATALOG by five coordinates --
%   departure key, arrival key, departure phase, arrival phase, thrust --
%   and returns the minimum flight time (interpolated between the thrust
%   rungs that bracket the request) plus a converged z8 seed from the
%   nearest stored rung.
%
%   GTO-CATALOG COPY (fix round 1, 2026-08-29 review): this file is
%   costate_catalog_gto_tulip's OWN copy of the shared picker -- every
%   deliverable ships its own, and this one is NOT a byte-identical clone
%   of the period-keyed copies (DRO/HALO/DPO/HALO_HALO). Two GTO-specific
%   differences from the generic template, both required because the GTO
%   departure "family" is a pseudo-family with NO orbital period (see
%   costate_common/get_family_orbit.m's 'gto' case):
%
%   (1) DEPARTURE KEY IS orientDeg (DEGREES), NOT A PERIOD. sheets(k).
%       tauDRO / .tau_dep hold orientDeg in {0,90,180,270} for this
%       catalog (the review-mandated sheet-key override, Stage B task 6).
%       orientation degrees are LINEARLY spaced, so the nearest-sheet
%       metric on the departure axis here is PLAIN LINEAR (ABSOLUTE/
%       SQUARED) DIFFERENCE, unlike the period-keyed catalogs' LOG-
%       DISTANCE convention (log tau_dep). Log-distance is zero-UNSAFE:
%       log(0) = -Inf, so orientDeg = 0 produced a NaN self-distance in
%       the generic picker and a SILENT WRONG-SHEET return for exactly
%       that (a full quarter of the sheet space). The Np / tau_arr axis
%       is untouched -- tulip petal counts and halo/dpo periods are never
%       zero, so the generic log/linear choice there is unaffected.
%
%   (2) DEPARTURE-PHASE DAY CONVERSION uses the sheet's TRUE GTO Kepler
%       period, reconstructed from dep_params.sma_km -- NOT sh.tauDRO
%       (which is orientDeg here, not a period; using it as one would
%       divide by zero at orientDeg = 0 and use nonsense units at every
%       other orientation). The Kepler-period formula below is copied
%       verbatim from get_family_orbit.m's 'gto' case (self-contained by
%       design: deliverable pickers do not depend on the source tree).
%
%   ARRIVAL KEY (schema v2, family-aware): for TULIP-arrival catalogs (this
%   one) the third argument is the petal count Np. For catalogs whose
%   arrival family is not tulip, it is the ARRIVAL PERIOD (ND) and sheets
%   are selected by nearest (departure-axis metric above, log tau_arr).
%   The mode is read off the catalog itself.
%
%   HONESTY CONTRACT (the warning flag): whenever what is RETURNED differs
%   from what was REQUESTED, a WARNING states exactly what you are getting:
%     - a different sheet (nearest orientation and/or petal count),
%     - a snapped phase pair (requests between grid points),
%     - the nearest SOLVED pair (the requested cell is unsolved/red),
%     - a seed from a different thrust rung.
%   The info output always carries requested-vs-delivered, warnings on or
%   off.
%
%  ASSUMPTIONS / NOTES:
%
% • Phases are days past each orbit's reference point (the family getters),
%   modulo the period. The departure "period" for phase-snapping purposes
%   is the sheet's TRUE GTO Kepler period (see (2) above), not tauDRO.
% • For thrusts between rungs, fly the returned z8's own rung or hand it to
%   tfMin at your thrust and let it converge the small correction.
% • pm = +1 tulips: z-mirror the pm = -1 solution (see cat_.derive).
% • This copy assumes dep_family = 'gto' and dep_params.sma_km present on
%   every sheet (true for costate_catalog_gto_tulip since Stage B task 6);
%   it is not the family-agnostic picker.
%
%% Inputs:
%
%  cat_                     struct                  costate_catalog_gto_tulip
%                                                   (or a struct shaped like
%                                                   it)
%
%  tauDep                   double                  Requested departure
%                                                   orientation (DEGREES)
%
%  arrKey                   double                  Petal count Np (tulip
%                                                   arrivals) or arrival
%                                                   period (ND) otherwise
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
%  M. Casey                                                   (c) 08/29/2026
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
   if isfield(sh,'arr_family') && ~strcmpi(sh.arr_family,'tulip')
       a2 = sh.tau_arr*1.07;          % arrival-period keying
   else
       a2 = sh.Np + 1;                % petal-count keying
   end
[tf,~,inf_] = costate_catalog_pick(cat_, sh.tauDRO + 5, a2, ...
                                   3.3, 11.0, 4.0);
   fprintf('tf = %.4f ND;  delivered sheet orient=%g Np arrKey=%.4g\n', ...
           tf, inf_.delivered.tauDRO, inf_.delivered.arrKey);
   return;
end
if ~exist('warnFlag','var'), warnFlag = true; end
    tStar = cat_.constants.tStar_s;
   warned = false;

%% Nearest sheet -- keying mode read off the catalog (schema v2). The
%  DEPARTURE axis is orientDeg (deg) for this catalog: LINEAR (squared)
%  difference, zero-safe -- see header note (1). The arrival axis keeps
%  whatever the generic template does (log tau_arr for period-keyed
%  arrivals; petal-count weighting for tulip arrivals), unchanged:
arrByPeriod = isfield(cat_.sheets, 'arr_family') && ...
              ~strcmpi(cat_.sheets(1).arr_family, 'tulip');
if arrByPeriod
    % departure axis: LINEAR (deg), zero-safe; arrival axis: log tau_arr
     dist = ([cat_.sheets.tau_dep] - tauDep).^2 ...
          + (log([cat_.sheets.tau_arr]) - log(arrKey)).^2;
   [~, ks] = min(dist);
       sh = cat_.sheets(ks);
    if warnFlag && (abs(sh.tau_dep - tauDep) > 1e-9 || ...
                    abs(sh.tau_arr - arrKey) > 1e-9)
        fprintf(['WARNING (costate_catalog_pick): no sheet at departure %.3f deg / arrival %.3f ND;\n', ...
                 '  you are getting the NEAREST SHEET: departure %.4g deg, arrival %.4g ND.\n'], ...
            tauDep, arrKey, sh.tau_dep, sh.tau_arr);
        warned = true;
    end
else
    % Petal-count keying. Departure axis: LINEAR (deg), zero-safe --
    % NOT log(tauDRO), which is -Inf at orientDeg = 0 (a quarter of this
    % catalog's key space) and silently mis-selects the sheet:
      nps = [cat_.sheets.Np];
     dist = ([cat_.sheets.tauDRO] - tauDep).^2 ...
          + 0.5*(nps - arrKey).^2;               % petal mismatch weighted
   [~, ks] = min(dist);
       sh = cat_.sheets(ks);
    if warnFlag && (abs(sh.tauDRO - tauDep) > 1e-9 || sh.Np ~= arrKey)
        fprintf(['WARNING (costate_catalog_pick): no sheet at orientation %.3f deg / %d petals;\n', ...
                 '  you are getting the NEAREST SHEET: orientation %.2f deg, Np = %d.\n'], ...
            tauDep, arrKey, sh.tauDRO, sh.Np);
        warned = true;
    end
end

%% Phase pair on that sheet's torus. Departure "period" for phase-snapping
%  is the sheet's TRUE GTO Kepler period, reconstructed from
%  dep_params.sma_km -- NOT sh.tauDRO (orientDeg here, not a period; the
%  generic template's Pd = sh.tauDRO*tStar/86400 divides by zero at
%  orientDeg = 0 and is nonsense-units at every other orientation).
%  Formula copied verbatim from get_family_orbit.m's 'gto' case:
       muE  = 6.67384e-20*(1 - cat_.constants.muStar)*(5.9736E24 + 7.35E22);
       Tkep = 2*pi*sqrt((sh.dep_params.sma_km)^3/muE);      % Kepler period, s
       Pd = Tkep/86400;                          % TRUE departure period (days)
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

if arrByPeriod, dArr = sh.tau_arr; else, dArr = sh.Np; end
info = struct( ...
    'requested', struct('tauDRO',tauDep,'Np',arrKey,'arrKey',arrKey, ...
                        'depDays',depPhaseDays, ...
                        'arrDays',arrPhaseDays,'thrustN',thrustN), ...
    'delivered', struct('tauDRO',sh.tauDRO,'Np',sh.Np,'arrKey',dArr, ...
                        'depDays',sh.sD_frac(iD)*Pd,'arrDays',sh.sA_frac(iA)*Pa, ...
                        'tfThrustN',min(max(thrustN,min(rAvail)),max(rAvail)), ...
                        'seedThrustN',rung(kNear)), ...
    'warned', warned);
end
