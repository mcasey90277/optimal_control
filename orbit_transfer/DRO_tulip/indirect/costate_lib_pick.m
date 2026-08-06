function [tf_days, e, bracket] = costate_lib_pick(lib, depPhaseDays, arrPhaseDays, thrustN, warnFlag)
%% Purpose:
%
%   Selects a costate-library entry by orbit phasing AND thrust level, and
%   returns the minimum flight time interpolated linearly between the two
%   thrust rungs that bracket the requested thrust. Phases are periodic, so
%   the phase pair is matched on the torus (with wraparound), each axis in
%   units of its own orbit period.
%
%   The returned entry holds a converged minimum-time PMP solution at the
%   NEAREST available rung: its field z8 = [lambda_r(3); lambda_v(3);
%   lambda_m; tf] can be handed directly to pumpkyn.cr3bp.tfMin (it is
%   accepted unchanged) or flown with pumpkyn.cr3bp.tfMinProp. For a thrust
%   between rungs, use z8 as the initial guess at the requested thrust and
%   let tfMin converge the small correction.
%
%  ASSUMPTIONS / NOTES:
%
% • Phases are measured in days past each orbit's reference point (the state
%   returned by getDRO / getTulip), modulo the orbit period.
% • Interpolation is linear in thrust between bracketing rungs; the rung set
%   was chosen so that stays accurate (denser where t_f(T) steepens, below
%   2 N). Outside the rung range the nearest rung's value is returned.
% • lib.grid.has_solution(i,j,k) is the availability map over
%   (departure phase x arrival phase x thrust).
%
%% Inputs:
%
%  lib                      struct                  costate_lib_dro_tulip_v2
%                                                   loaded from its .mat
%
%  depPhaseDays             double                  Departure phase along
%                                                   the DRO, days past the
%                                                   reference point
%
%  arrPhaseDays             double                  Arrival phase along the
%                                                   tulip, days past the
%                                                   reference point
%
%  thrustN                  double                  Requested thrust (N),
%                                                   nominally within the
%                                                   library's rung range
%
%  warnFlag                 logical                 true (default): print a
%                                                   WARNING whenever what is
%                                                   RETURNED differs from
%                                                   what was REQUESTED --
%                                                   off-grid phases snapped
%                                                   to the nearest cell, an
%                                                   unsolved (red/grey) cell
%                                                   answered by the nearest
%                                                   SOLVED pair, or a thrust
%                                                   served from a different
%                                                   rung. The message states
%                                                   the phase pair actually
%                                                   delivered.
%
%% Outputs:
%
%  tf_days                  double                  Minimum flight time at
%                                                   the requested thrust,
%                                                   linearly interpolated
%                                                   between rungs (days)
%
%  e                        struct                  Library entry at the
%                                                   nearest rung:
%                                                   -------
%                                                   lambda0        [7 x 1]
%                                                   tf_nd, tf_days
%                                                   z8             [8 x 1]
%                                                   thrust_N, isp_s, m0_kg
%                                                   departure_phase_*
%                                                   arrival_phase_*
%
%  bracket                  struct                  .thrust_N [1 x 2] and
%                                                   .tf_days  [1 x 2] of the
%                                                   rungs used (NaN where a
%                                                   rung is unavailable)
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: mid-band thrust at an arbitrary phasing:
       L = load('costate_lib_dro_tulip_v2.mat');
     lib = L.costate_lib_dro_tulip_v2;
[tf,e,br] = costate_lib_pick(lib, 1.0, 10.0, 4.0);
   fprintf('4.0 N: tf = %.3f days (rungs %.1f/%.1f N -> %.3f/%.3f d)\n', ...
           tf, br.thrust_N(1), br.thrust_N(2), br.tf_days(1), br.tf_days(2));
   fprintf('nearest-rung seed at %.1f N, z8 = [%s]\n', e.thrust_N, ...
           sprintf('%.4f ', e.z8));
   return;
end

if ~exist('warnFlag','var'), warnFlag = true; end

%% Locate the phase pair on the torus:
      Pd = lib.departure_params.period_days;
      Pa = lib.arrival_params.period_days;
      fd = mod(depPhaseDays/Pd, 1);
      fa = mod(arrPhaseDays/Pa, 1);
     sDg = lib.grid.departure_phase_frac;
     sAg = lib.grid.arrival_phase_frac;
      dd = abs(sDg - fd);  dd = min(dd, 1-dd);
      da = abs(sAg - fa);  da = min(da, 1-da);
   [~,iD] = min(dd);
   [~,iA] = min(da);

%% Off-grid snap? Report how far the returned cell is from the request:
      Pd = lib.departure_params.period_days;     %#ok<NASGU> (already set)
   dSnap = min(abs(sDg(iD)-fd), 1-abs(sDg(iD)-fd)) * lib.departure_params.period_days;
   aSnap = min(abs(sAg(iA)-fa), 1-abs(sAg(iA)-fa)) * lib.arrival_params.period_days;
if warnFlag && (dSnap > 1e-6 || aSnap > 1e-6)
    fprintf(['WARNING (costate_lib_pick): requested phases (dep %.4f d, arr %.4f d) are OFF-GRID;\n', ...
             '  returning the nearest grid pair (dep %.4f d, arr %.4f d) -- %.2f h and %.2f h away.\n'], ...
        depPhaseDays, arrPhaseDays, ...
        sDg(iD)*lib.departure_params.period_days, sAg(iA)*lib.arrival_params.period_days, ...
        24*dSnap, 24*aSnap);
end

%% Which rungs have a solution at this phase pair:
    rung = lib.grid.thrust_N;
     hav = squeeze(lib.grid.has_solution(iD,iA,:))';
if ~any(hav)
    % this grid pair is red/grey everywhere: answer from the nearest pair
    % that IS solved somewhere, and say so loudly
    [okD, okA] = find(any(lib.grid.has_solution,3));
    dd2 = min(abs(sDg(okD)'-fd),1-abs(sDg(okD)'-fd)).^2 ...
        + min(abs(sAg(okA)'-fa),1-abs(sAg(okA)'-fa)).^2;
    [~, kb] = min(dd2);
    iD = okD(kb);  iA = okA(kb);
    hav = squeeze(lib.grid.has_solution(iD,iA,:))';
    if warnFlag
        fprintf(['WARNING (costate_lib_pick): the requested phase pair has NO SOLUTION in this\n', ...
                 '  library (unsolved cell). You are getting the nearest SOLVED pair instead:\n', ...
                 '  dep %.4f d, arr %.4f d.\n'], ...
            sDg(iD)*lib.departure_params.period_days, ...
            sAg(iA)*lib.arrival_params.period_days);
    end
end
   rAvail = rung(hav);
  tfAvail = squeeze(lib.grid.tf_days(iD,iA,:))';
  tfAvail = tfAvail(hav);

%% Linear interpolation in thrust between bracketing rungs:
if thrustN <= min(rAvail)
    tf_days = tfAvail(rAvail == min(rAvail));
    bracket = struct('thrust_N',[min(rAvail) NaN], 'tf_days',[tf_days NaN]);
elseif thrustN >= max(rAvail)
    tf_days = tfAvail(rAvail == max(rAvail));
    bracket = struct('thrust_N',[max(rAvail) NaN], 'tf_days',[tf_days NaN]);
else
       lo = max(rAvail(rAvail <= thrustN));
       hi = min(rAvail(rAvail >= thrustN));
     tfLo = tfAvail(rAvail == lo);
     tfHi = tfAvail(rAvail == hi);
    if hi == lo
        tf_days = tfLo;
    else
        tf_days = tfLo + (tfHi - tfLo)*(thrustN - lo)/(hi - lo);
    end
    bracket = struct('thrust_N',[lo hi], 'tf_days',[tfLo tfHi]);
end

%% Nearest available rung supplies the seed (its entry carries deltaV_kms,
%  m_final_kg and propellant_kg for that rung):
    [~,kn] = min(abs(rAvail - thrustN));
      kNear = find(rung == rAvail(kn), 1);
         ei = lib.grid.entry_index(iD,iA,kNear);
          e = lib.entries(ei);
if warnFlag && abs(e.thrust_N - thrustN) > 1e-9
    fprintf(['WARNING (costate_lib_pick): no stored rung at %.3f N for this pair;\n', ...
             '  the returned SEED is the %.2f N entry (t_f there %.4f d). The interpolated\n', ...
             '  flight time above still refers to your requested thrust.\n'], ...
        thrustN, e.thrust_N, e.tf_days);
end

end
