function test_constThrust()
%% Purpose:
%
%  Verify the constant-mass-flow propulsion model and the booster parameter
%  set it is driven with: vacuum thrust, the ambient back-pressure debit, a
%  propellant flow that is fixed by the choked-nozzle assumption and therefore
%  independent of ambient pressure, and the clamped branch, in which thrust
%  and flow must go to zero TOGETHER.
%
%  Every expected value below is a literal worked out by hand from the
%  documented placeholders in coorbital.util.boosterDefaults, with the
%  arithmetic shown. None of them is obtained by re-running the function's own
%  formula on the function's own outputs, which would only prove the code
%  agrees with itself.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               bst = coorbital.util.boosterDefaults();
               veh = coorbital.util.vehicleDefaults();

%% Required booster fields exist and are positive:
          required = {'massDry','massProp','thrustVac','Isp','Aexit', ...
                      'Sref','CL','LD'};
    for k = 1:numel(required)
        assert(isfield(bst,required{k}),'boosterDefaults missing %s',required{k});
        assert(bst.(required{k}) > 0,'%s must be positive',required{k});
    end

%% There must be NO massWet field -- the mass bookkeeping rule forbids it
%% because it cannot be read unambiguously (booster only, or stack including
%% payload?). See docs/plan_2026-08-07_boost_descent_chain.md:
    assert(~isfield(bst,'massWet'), ...
        'boosterDefaults must not define massWet; add massDry and massProp instead');

%% Anchor the placeholder values. Every hand-computed literal further down was
%% derived from exactly these numbers, so if the parameter set is retuned this
%% is the assertion that fires first and says why the rest did:
    assert(bst.massDry   == 1500,   'expected massDry   = 1500 kg');
    assert(bst.massProp  == 30000,  'expected massProp  = 30000 kg');
    assert(bst.thrustVac == 950000, 'expected thrustVac = 950000 N');
    assert(bst.Isp       == 260,    'expected Isp       = 260 s (vacuum)');
    assert(bst.Aexit     == 1.25,   'expected Aexit     = 1.25 m^2');
    assert(veh.mass      == 900,    'expected payload mass = 900 kg');
    assert(c.g0          == 9.80665,'expected g0 = 9.80665 m/s^2');

%% ---- Vacuum thrust ------------------------------------------------------
%% At P = 0 the back-pressure debit vanishes identically, so the delivered
%% thrust must be the vacuum thrust to the last bit, not merely to tolerance:
    [Tvac,mdotVac] = coorbital.prop.constThrust(0,0,bst);
    assert(Tvac == 950000,'vacuum thrust must equal thrustVac exactly');

%% ---- Back-pressure debit ------------------------------------------------
%% Two different "sea levels" are in play and they are NOT the same number, so
%% both are tested deliberately rather than one being assumed to cover the
%% other:
%%
%%   (a) The US Standard Atmosphere 1976 sea-level pressure, 101325 Pa. This is
%%       the condition against which a motor's sea-level thrust is quoted in
%%       the literature, so it is the number a reader will check the booster
%%       sizing against.
%%
%%   (b) coorbital.atmos.expAtmos(0) = 87909.92 Pa. The library atmosphere is
%%       isothermal at 250 K, so P(0) = rho0*Rair*T0 = 1.225*287.0528*250,
%%       some 13% below (a). This is the pressure the equations of motion will
%%       actually hand this model at lift-off, so it is the one that sets the
%%       trajectory.
%%
%% (a) hand arithmetic: 1.25 m^2 * 101325 Pa = 126656.25 N debit;
%%                      950000 - 126656.25   = 823343.75 N
      [T76,mdot76] = coorbital.prop.constThrust(0,101325,bst);
    assert(abs(T76 - 823343.75) < 1e-6, ...
        'thrust at the US76 sea-level pressure should be 823343.75 N');

%% (b) hand arithmetic: 1.225*287.0528*250   = 87909.92 Pa;
%%                      1.25 m^2 * 87909.92  = 109887.40 N debit;
%%                      950000 - 109887.40   = 840112.60 N
          [~,Pexp] = coorbital.atmos.expAtmos(0);
    assert(abs(Pexp - 87909.92) < 1e-9, ...
        'expAtmos sea-level pressure should be 87909.92 Pa');
    [Texp,mdotExp] = coorbital.prop.constThrust(0,Pexp,bst);
    assert(abs(Texp - 840112.60) < 1e-6, ...
        'thrust at the expAtmos sea-level pressure should be 840112.60 N');

%% Thrust falls monotonically as ambient pressure rises, and the debit is the
%% only thing that changed:
    assert(Texp < Tvac && T76 < Texp, ...
        'thrust must decrease as ambient pressure increases');

%% ---- Mass flow ----------------------------------------------------------
%% Hand arithmetic: exhaust speed Isp*g0 = 260*9.80665 = 2549.729 m/s;
%%                  mdot = thrustVac/(Isp*g0)
%%                       = 950000/2549.729 = 372.5886162803969 kg/s
      mdotExpected = 372.5886162803969;
    assert(abs(mdotVac - mdotExpected) < 1e-9, ...
        'vacuum mdot should be 372.5886162803969 kg/s');
    assert(mdotVac > 0, ...
        'mdot is positive by convention; the equations of motion apply the minus');

%% The choked-nozzle assumption is the whole reason this model is
%% self-consistent: back pressure changes the delivered thrust, NOT the flow.
%% Bit-identical, not merely close -- any pressure dependence at all is wrong:
    assert(mdot76 == mdotVac && mdotExp == mdotVac, ...
        'mdot must not depend on ambient pressure');

%% ---- Burn time ----------------------------------------------------------
%% Hand arithmetic: massProp/mdot = 30000/372.5886162803969
%%                                = 80.51775789473683 s
          burnTime = bst.massProp./mdotVac;
    assert(abs(burnTime - 80.51775789473683) < 1e-9, ...
        'burn time should be 80.51775789473683 s');
    assert(burnTime > 0,'burn time must be positive');

%% ---- Booster sizing sanity ----------------------------------------------
%% These check the performance figures quoted in the boosterDefaults header,
%% so the header cannot drift away from the numbers underneath it.
%%
%% Hand arithmetic: liftoff mass   = 900 + 1500 + 30000 = 32400 kg
%%                  burnout mass   = 900 + 1500         =  2400 kg
%%                  liftoff weight = 32400*9.80665      = 317735.46 N
%%                  T/W            = 950000/317735.46   = 2.989908649163679
          mLiftoff = veh.mass + bst.massDry + bst.massProp;
          mBurnout = veh.mass + bst.massDry;
    assert(mLiftoff == 32400,'liftoff mass should be 32400 kg');
    assert(mBurnout == 2400, 'burnout mass should be 2400 kg');
       liftoffTtoW = bst.thrustVac./(mLiftoff.*c.g0);
    assert(abs(liftoffTtoW - 2.989908649163679) < 1e-9, ...
        'liftoff T/W should be 2.989908649163679');
    assert(liftoffTtoW > 2 && liftoffTtoW < 3, ...
        'liftoff T/W must sit in the design band of roughly 2 to 3');

%% Ideal (loss-free) Tsiolkovsky delta-V.
%% Hand arithmetic: mass ratio = 32400/2400     = 13.5;
%%                  dV = 2549.729*log(13.5)     = 6636.153368978423 m/s
           dVIdeal = bst.Isp.*c.g0.*log(mLiftoff./mBurnout);
    assert(abs(dVIdeal - 6636.153368978423) < 1e-9, ...
        'ideal delta-V should be 6636.153368978423 m/s');
    assert(dVIdeal > 6000, ...
        'ideal delta-V must leave room for losses and still reach entry speed');

%% ---- Time independence --------------------------------------------------
%% Thrust is constant through the burn, so the time argument must not move
%% either output:
  [Tlate,mdotLate] = coorbital.prop.constThrust(1e4,Pexp,bst);
    assert(Tlate == Texp && mdotLate == mdotExp, ...
        'constThrust must not depend on t');

%% ---- Clamped branch -----------------------------------------------------
%% Drive the back-pressure debit far past the vacuum thrust with an absurd
%% exit area: 1e6 m^2 * 101325 Pa = 1.01325e11 N of debit against 9.5e5 N of
%% vacuum thrust. The model has left its domain of validity, so the motor must
%% report as not running -- BOTH outputs zero, together. A zero thrust with a
%% live propellant flow would silently drain the mass state for nothing:
         bstAbsurd = bst;
   bstAbsurd.Aexit = 1e6;
[Tclamp,mdotClamp] = coorbital.prop.constThrust(0,101325,bstAbsurd);
    assert(Tclamp == 0,'clamped thrust must be exactly zero, never negative');
    assert(mdotClamp == 0, ...
        'mdot must be exactly zero whenever the thrust clamp fires');

%% ---- Vectorised evaluation ----------------------------------------------
%% A column of pressures must return the same answers the scalar calls gave,
%% element for element, so the equations of motion can be run over a whole
%% time grid:
              Pcol = [0; Pexp; 101325];
    [Tcol,mdotCol] = coorbital.prop.constThrust(0,Pcol,bst);
    assert(isequal(size(Tcol),[3 1]) && isequal(size(mdotCol),[3 1]), ...
        'vectorised call must preserve the shape of P');
    assert(isequal(Tcol,[Tvac; Texp; T76]), ...
        'vectorised thrust must match the scalar evaluations');
    assert(isequal(mdotCol,[mdotVac; mdotVac; mdotVac]), ...
        'vectorised mdot must match the scalar evaluations');
end
