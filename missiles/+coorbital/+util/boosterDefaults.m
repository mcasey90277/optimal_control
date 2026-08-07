function bst = boosterDefaults()
%% Purpose:
%
%  Default single-stage booster parameters for lofting the generic glide
%  vehicle of coorbital.util.vehicleDefaults. PLACEHOLDER VALUES: they are
%  open-literature, order-of-magnitude figures for a large solid-propellant
%  first stage and are not any specific motor. Replace before any result
%  leaves this machine.
%
%  Mass bookkeeping (docs/plan_2026-08-07_boost_descent_chain.md, "Mass
%  bookkeeping") -- the propagated state mass is ALWAYS the total mass
%  currently carried:
%
%     m at liftoff       = veh.mass + bst.massDry + bst.massProp
%     m at burnout       = veh.mass + bst.massDry
%     m after separation = veh.mass
%     burn time          = bst.massProp / mdot
%
%  There is deliberately NO massWet field. A wet mass would have to mean
%  either booster-only or stack-including-payload, and either reading is
%  defensible, so the field would silently invite the one ambiguity the
%  bookkeeping rule above exists to remove. Add the two masses that are here.
%
%  The aerodynamic fields describe the BOOSTED STACK, not the glide vehicle:
%  a slender body of revolution flying near zero incidence, which is a much
%  poorer lifting shape than the waverider it carries. They exist so the boost
%  phase can fly with coorbital.aero.constLD, the same aerodynamic model the
%  glide uses. As in vehicleDefaults there is no CD field -- constLD derives
%  CD = CL/LD, so drag cannot fall out of sync with the pair below.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  bst              Struct                      Booster parameters:
%                                               massDry   (kg) spent booster
%                                                  structure, jettisoned at
%                                                  separation; EXCLUDES the
%                                                  payload
%                                               massProp  (kg) usable
%                                                  propellant
%                                               thrustVac (N) vacuum thrust,
%                                                  constant through the burn
%                                               Isp       (s) VACUUM specific
%                                                  impulse; mdot is derived
%                                                  from it, see the note below
%                                               Aexit     (m^2) nozzle exit
%                                                  area, sets the ambient
%                                                  back-pressure thrust debit
%                                               Sref      (m^2) stack
%                                                  reference area
%                                               CL        (-) stack lift
%                                                  coefficient
%                                               LD        (-) stack
%                                                  lift-to-drag ratio
%
%% Note -- Isp is the VACUUM value, and that matters:
%
%  coorbital.prop.constThrust derives the propellant mass flow as
%  mdot = thrustVac/(Isp*g0). That identity holds only when the thrust and the
%  specific impulse are quoted at the SAME ambient condition. Both are quoted
%  in vacuum here. Substituting a sea-level Isp against thrustVac would
%  overstate the flow by the sea-level thrust debit and burn the motor out
%  early, which reads as a trajectory error rather than a units error.
%
%% Note -- the resulting performance, computed from the values below:
%
%  With the 900 kg vehicleDefaults payload, and g0 = 9.80665 m/s^2:
%
%     liftoff mass    = 900 + 1500 + 30000        = 32400    kg
%     burnout mass    = 900 + 1500                =  2400    kg
%     mass ratio                                  =    13.5
%     structural coeff = 1500/(1500 + 30000)      =     0.0476  (4.76 %)
%     exhaust speed   = Isp*g0 = 260*9.80665      =  2549.729 m/s
%     mass flow       = 950000/2549.729           =   372.589 kg/s
%     burn time       = 30000/372.5886            =    80.518 s
%     liftoff T/W     = 950000/(32400*9.80665)    =     2.990
%     ideal delta-V   = 2549.729*log(13.5)        =  6636.2  m/s
%
%  The ideal delta-V is the loss-free Tsiolkovsky figure. What a real ascent
%  keeps of it depends on the PITCH PROGRAM flown, so the losses below are a
%  measurement of one trajectory and not a property of this booster.
%
%  Measured along the shipped BM/run_ballistic ascent -- 89 deg at the pad
%  pitching to 34 deg at burnout -- by integrating each term of the dV
%  equation over the flown boost phase:
%
%     gravity        int gr sin(gamma) dt              =  697.23 m/s
%     drag           int D/m dt                        =   48.37 m/s
%     back pressure  int (Tvac - T)/m dt               =   98.85 m/s
%     steering       int T(1 - cos(alpha))/m dt        =   33.51 m/s
%                                                        --------
%     total losses                                     =  877.96 m/s
%
%     burnout speed  = 10 + 6636.15 - 877.96           = 5768.20 m/s
%
%  which closes against the propagated burnout speed to 0.01 m/s. (The leading
%  10 m/s is the nonzero speed run_ballistic must start from, the equations
%  being singular at V = 0.)
%
%  An earlier version of this note quoted 500-700 m/s. That is too low for a
%  lofted ascent of this duration: gravity loss ALONE is 697 m/s here, because
%  a trajectory steep enough to reach a 1600 km apogee spends most of the burn
%  fighting gravity nearly head-on. A flatter program trades gravity loss for
%  drag loss and lands lower. Expect 700-900 m/s for anything lofted, and
%  measure it rather than assuming it.
%
%% Note -- the structural coefficient is OPTIMISTIC, and it is load-bearing:
%
%  massDry/(massDry + massProp) = 1500/31500 = 4.76 %. Real large solid
%  boosters run 10-15 % -- cases, nozzles, insulation, TVC actuators and
%  skirts do not scale away -- so this stage is light by roughly a factor of
%  two to three. The value is DELIBERATELY LEFT ALONE: it is a marked
%  placeholder, and changing it would move every headline number in
%  docs/README.md and every pinned literal in tests/test_constThrust.m for no
%  gain in what this library is demonstrating, which is the machinery and not
%  the motor. But it must be read as optimistic, because the ideal delta-V
%  above is what pays for the entry interface, and a light structure buys most
%  of it.
%
%  What a realistic structure would cost, at an 11.2 % coefficient with the
%  same 30000 kg propellant load and the same 900 kg payload:
%
%     massDry         = 30000*0.112/0.888          =  3783.78 kg
%     liftoff mass    = 900 + 3783.78 + 30000      = 34683.78 kg
%     burnout mass    = 900 + 3783.78              =  4683.78 kg
%     ideal delta-V   = 2549.729*log(7.4051)       =  5104.98 m/s
%
%  a fall of 1531 m/s from 6636 m/s. Re-splitting the SAME 31500 kg stack at
%  11.2 % instead of enlarging it gives 5074.50 m/s, so it is about 5.1 km/s
%  either way -- the counterfactual is insensitive to which mass is held.
%
%  So: with the losses measured above, this booster delivers a burnout speed
%  near 6 km/s, and the 60 km, 6 km/s entry interface in HGV/run_glide is
%  consistent with it. That alignment is an ARTEFACT OF THE OPTIMISTIC MASS
%  FRACTION, not an expectation to design against. A stage with a realistic
%  structure and this propellant load lands roughly 1.5 km/s short of the
%  ideal figure quoted here, and a burnout near 6 km/s would then need either
%  more propellant, a higher Isp, or a second stage.
%
%  Burn time and liftoff T/W are not independent: for a given propellant mass
%  fraction f, burn time = Isp*f/(T/W). Holding T/W in the 2-3 band with a
%  solid-class Isp therefore pins the burn near 80 s, and no choice of masses
%  shortens it without either raising T/W above the band or giving up the
%  delta-V needed to reach entry speed.
%
%% References:
%   [1] Sutton, G.P., Biblarz, O., "Rocket Propulsion Elements," 9th ed.,
%       Wiley, 2017, Ch. 2-3 and Ch. 12. (Isp, thrust coefficient, and
%       representative HTPB solid-motor performance.)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Masses -- see the bookkeeping block above; there is no massWet:
       bst.massDry = 1500;            %kg, PLACEHOLDER, spent structure only
      bst.massProp = 30000;           %kg, PLACEHOLDER, usable propellant

%% Propulsion -- thrustVac and Isp are BOTH vacuum values:
     bst.thrustVac = 950000;          %N, PLACEHOLDER, constant through the burn
           bst.Isp = 260;             %s, PLACEHOLDER, VACUUM specific impulse;
                                      %   typical HTPB solid VACUUM range
                                      %   250-270 s
         bst.Aexit = 1.25;            %m^2, PLACEHOLDER, nozzle exit area;
                                      %   1.26 m exit diameter

%% Boosted-stack aerodynamics for coorbital.aero.constLD (PLACEHOLDER):
          bst.Sref = 1.77;            %m^2, PLACEHOLDER, pi/4*(1.5 m)^2 body
                                      %   cross-section
            bst.CL = 0.05;            %PLACEHOLDER, near-zero-incidence body lift
            bst.LD = 0.25;            %PLACEHOLDER, gives CD = CL/LD = 0.20
end
