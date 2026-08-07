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
%
%% Revision History:
%  Michael Casey                                                08/07/2026
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
end
