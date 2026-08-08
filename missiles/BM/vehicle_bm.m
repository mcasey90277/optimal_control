function veh = vehicle_bm()
%% Purpose:
%
%  Ballistic-missile re-entry vehicle parameters for the run_ballistic entry
%  script. Starts from the library defaults so that anything not specific to
%  this airframe has exactly one home, and overrides the four values that are
%  specific to it.
%
%% Note -- this is NOT the glide vehicle, and the differences are the point:
%
%  coorbital.util.vehicleDefaults describes a lifting waverider: a large
%  reference area, a healthy lift coefficient, and L/D = 2.5. Flown ballistic
%  that vehicle does not fall, it glides, and the "ballistic" range would be
%  set by lift rather than by the Keplerian arc the boost put it on. The
%  re-entry body below is the opposite article -- a blunt, essentially
%  non-lifting cone with a high ballistic coefficient:
%
%     beta = mass/(CD*Sref) = 900/(0.25*0.385) = 9351 kg/m^2
%
%  which is what makes the trajectory genuinely ballistic and makes the
%  Keplerian cross-check in run_ballistic a meaningful comparison rather than
%  a comparison against a completely different flight regime.
%
%% Note -- L/D = 0.02 is not a typo:
%
%  There is no CD field to set. coorbital.aero.constLD derives CD = CL/LD, so
%  drag is fixed by the two values below and cannot fall out of sync with
%  them. Expressing "almost no lift, plenty of drag" in that parameterisation
%  therefore requires a SMALL CL over a SMALL L/D, and the ratio, not either
%  value alone, is what sets CD:
%
%     CD = CL/LD = 0.005/0.02 = 0.25,   lift/drag acceleration ratio = 0.02
%
%  CL cannot simply be set to zero, because that would drive CD to zero with
%  it and leave a vehicle that does not feel the atmosphere at all. The two
%  percent of residual lift is deliberate: a real cone at small incidence
%  carries a little, and run_ballistic separates its contribution from drag's
%  in the cross-check rather than letting the two hide inside one number.
%
%  The mass is deliberately left at the 900 kg of vehicleDefaults. The mass
%  bookkeeping documented in coorbital.util.boosterDefaults -- 32400 kg at
%  liftoff, 2400 kg at burnout, 900 kg after separation -- is quoted against
%  that payload, and changing it here would silently invalidate every number
%  in that header.
%
%% Note -- alphaMaxDeg lives HERE because it is a property of the VEHICLE:
%
%  The clamp on the magnitude of the angle of attack is a CONTROL-AUTHORITY
%  limit: how far off the relative wind this airframe may be flown before its
%  structure, its actuators or its aerothermal margins say no. It is not a
%  targeting degree of freedom, and it must not be raised because a desired
%  trajectory would otherwise be out of reach.
%
%  BM/run_ballistic and BM/run_ballistic_target both fly this vehicle, and
%  until 08/08/2026 they carried the limit as a USER PARAMETERS entry of their
%  own -- 6 deg in one and 12 deg in the other. Two different control
%  authorities for one airframe made the two scripts' performance
%  non-comparable, and the 12 deg had been chosen to bring a demonstration
%  target inside the depressed branch, which is exactly the reasoning a
%  vehicle limit may not be set by. Both scripts now read the value below, and
%  both still accept an explicit alphaMax override for a deliberate
%  sensitivity study.
%
%  6 deg IS A PLACEHOLDER AWAITING A QUALIFICATION BASIS, like every other
%  number in this file. It is the value BM/run_ballistic has always flown, and
%  it was not chosen to make anything reachable. What the choice costs is
%  written down rather than left to be rediscovered:
%
%    at 6 deg  maximum range 5211.525 km; the depressed branch spans roughly
%              4708 to 5212 km, and its max-range loft angle sits near
%              -42.9 deg, so a loft bracket must reach below that to find it;
%    at 12 deg maximum range falls by about 156 km, and the depressed branch
%              reaches down to about 1684 km.
%
%  Raising the limit therefore BUYS depressed-branch reach and PAYS about
%  156 km of maximum range. Neither number is a reason to move it: the
%  qualified value is whatever the airframe is cleared for, and that number
%  does not exist yet.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  veh              Struct                      Vehicle parameters; see
%                                               coorbital.util.vehicleDefaults
%                                               mass       (kg)
%                                               Sref       (m^2)
%                                               CL         (-)
%                                               LD         (-)
%                                               noseRadius (m)
%                                               alphaMaxDeg(deg) clamp on the
%                                                          MAGNITUDE of the
%                                                          angle of attack
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  alphaMaxDeg moved here from the two scripts   08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               veh = coorbital.util.vehicleDefaults();

%% Re-entry-body overrides (PLACEHOLDER values, open-literature magnitudes for
%% a generic conical RV; replace before any result leaves this machine):
          veh.mass = 900;              %kg, PLACEHOLDER, matches the payload mass
                                       %   assumed by coorbital.util.boosterDefaults
          veh.Sref = 0.385;            %m^2, PLACEHOLDER, pi/4*(0.7 m)^2 base area
            veh.CL = 0.005;            %PLACEHOLDER, near-zero-incidence cone lift
            veh.LD = 0.02;             %PLACEHOLDER, NOT a typo -- see the header note;
                                       %   gives CD = CL/LD = 0.25
    veh.noseRadius = 0.10;             %m, PLACEHOLDER, blunted nose; unused until
                                       %   Sutton-Graves heating arrives
   veh.alphaMaxDeg = 6;                %deg, PLACEHOLDER control-authority limit on
                                       %   |angle of attack| [1 .. 15], AWAITING A
                                       %   QUALIFICATION BASIS. Read by BM/run_ballistic
                                       %   and BM/run_ballistic_target alike -- one
                                       %   vehicle, one limit. See the header note for
                                       %   what the value costs and why it is not a
                                       %   targeting degree of freedom
end
