function eomFn = massConstant(baseEom)
%% Purpose:
%
%  Adapt a six-state atmospheric-flight EOM to the seven-state vector used by
%  a chain that includes a powered phase, by appending dm/dt = 0. Unpowered
%  phases carry mass so that every phase in a chain shares one state vector,
%  which removes the need for a state mapping at the junctions.
%
%% WARNING -- x(7) is CARRIED, not READ, and that is a trap:
%
%  The wrapped equations of motion are six-state. They take the vehicle mass
%  from veh.mass and they never look at x(7). The seventh component is
%  bookkeeping: it rides along so the chain shares one state vector, and it
%  drives nothing.
%
%  So there are two sources of truth for one physical quantity, and nothing
%  about the arithmetic makes them agree. A chain that jettisons a booster
%  through a phase link -- dropping x(7) from the stack mass to the payload
%  mass -- while continuing to hand the same veh struct to the next phase will
%  keep flying at the OLD weight, silently, with a mass history that says
%  otherwise. Nothing diverges, nothing warns, and the trajectory looks
%  entirely plausible.
%
%  This wrapper therefore refuses to run unless the two agree, to
%  1e-9 relative (never tighter than 1e-9 absolute), which is loose enough to
%  clear the residual an ODE burnout-event solve leaves on the mass state --
%  measured at 3.4e-10 kg on a 900 kg payload -- and tight enough that a real
%  staging jump of hundreds of kilograms cannot slip through. THE CALLER MUST
%  REBUILD veh AFTER EVERY STAGING EVENT.
%
%  Rebuilding veh.mass alone is NOT enough and must not be automated here by
%  injecting x(7) into a copy of veh: that would silently leave Sref, CL and
%  LD describing the jettisoned stack, which is worse than the bug it fixes
%  because it looks repaired. Separation changes the whole airframe. The
%  structural answer is a per-phase vehicle carried on the phase struct; until
%  that exists, this guard is what stands in for it.
%
%% Inputs:
%
%  baseEom          Function handle             Six-state EOM with the
%                                               signature xdot = f(t,x,u,veh,env)
%
%% Outputs:
%
%  eomFn            Function handle             Seven-state EOM with the same
%                                               signature; component 7 is mass
%                                               (kg) and its derivative is
%                                               zero. Raises
%                                               coorbital:massConstant:massMismatch
%                                               if x(7) disagrees with
%                                               veh.mass, see the WARNING above
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Guard x(7) against veh.mass                    08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
                fD = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
                xD = [c.rE + 40e3; 0; 0; 4000; deg2rad(-2); deg2rad(90); veh.mass];
             xdotD = fD(0,xD,[0;0],veh,env);
    fprintf('dm/dt = %.1f kg/s over a %d-state vector\n',xdotD(7),numel(xdotD));
             eomFn = [];
    return;
end

             eomFn = @(t,x,u,veh,env) massChecked(baseEom,t,x,u,veh,env);
end

function xdot = massChecked(baseEom,t,x,u,veh,env)
%% Purpose:
%
%  Append dm/dt = 0 to a six-state EOM, having first confirmed that the mass
%  the caller is CARRYING in x(7) is the mass the wrapped equations will
%  actually FLY, veh.mass. See the WARNING in the parent header: those two are
%  independent, the equations read only the second, and a staging jump applied
%  to the first alone is invisible without this check.
%
%  The check is inside the derivative call, not at wrap time, because the
%  divergence is created mid-chain by a phase link and there is nothing to
%  compare when the handle is built.
%
%% Inputs:
%
%  baseEom          Function handle             Six-state EOM,
%                                               xdot = f(t,x,u,veh,env)
%
%  t                scalar                      Time since phase start (s)
%
%  x                [7 x 1]                     State; x(7) is mass (kg)
%
%  u                [nu x 1]                    Control, passed through
%
%  veh              Struct                      Vehicle parameters; the mass
%                                               (kg) field is read here and is
%                                               what the base EOM will divide by
%
%  env              Struct                      Environment, passed through
%
%% Outputs:
%
%  xdot             [7 x 1]                     State derivative; component 7
%                                               is exactly zero
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The contract, checked before anything is believed:
if numel(x) < 7
    error('coorbital:massConstant:stateWidth', ...
        ['A %d-state vector was handed to the seven-state wrapper; the mass ' ...
         'must be component 7.'],numel(x));
end
if ~isfield(veh,'mass')
    error('coorbital:massConstant:noMass', ...
        ['The vehicle struct has no mass field. The wrapped equations of ' ...
         'motion divide by veh.mass, so it must be present and must equal ' ...
         'the carried mass x(7) = %.9f kg.'],x(7));
end
if ~isscalar(veh.mass)
    error('coorbital:massConstant:massNotScalar', ...
        ['veh.mass has %d elements; a point-mass vehicle carries one mass, ' ...
         'and a non-scalar would turn the comparison below into an ' ...
         'array test that quietly passes on a partial match.'],numel(veh.mass));
end

%% Relative tolerance, floored so a near-zero veh.mass cannot make it
%% vanishingly tight. Sized to clear the residual an ODE event solve leaves on
%% the mass state, not to police the integrator:
              mVeh = veh.mass(1);
               tol = 1e-9;
if mVeh > 1
               tol = tol.*mVeh;
end
if abs(x(7) - mVeh) > tol
    error('coorbital:massConstant:massMismatch', ...
        ['The carried mass x(7) = %.9f kg disagrees with veh.mass = ' ...
         '%.9f kg by %.6g kg, over a %.6g kg budget. The wrapped equations ' ...
         'of motion divide by veh.mass and NEVER read x(7), so this flight ' ...
         'would run silently at the wrong weight. Rebuild the vehicle struct ' ...
         '-- mass AND aerodynamics -- after every staging event.'], ...
        x(7),mVeh,abs(x(7) - mVeh),tol);
end

              xdot = [baseEom(t,x(1:6),u,veh,env); 0];
end
