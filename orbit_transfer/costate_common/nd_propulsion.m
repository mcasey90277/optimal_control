function p = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar)
%% Purpose:
%
%   THE nondimensional propulsion conversion for the CR3BP campaigns: a
%   thruster in engineering units (N, s, kg) becomes the two numbers every
%   solver here actually takes -- the ND exhaust speed c and the ND thrust
%   acceleration at unit mass fraction T.
%
%   About twenty files wrote these three lines out by hand, identically:
%
%       g0  = 9.80665*tStar^2/(1000*lStar);
%       cnd = (ispS/tStar)*g0;
%       Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);
%
%   One of them already carried the comment "same ND thrust/exhaust
%   conversion thrust_ladder_library uses", which is how a copy announces
%   itself. The expressions here are those, character for character, so
%   every migrated site keeps its value to the last bit.
%
%  ASSUMPTIONS / NOTES:
%
% • `Tnd` is the thrust acceleration at MASS FRACTION ONE, matching
%   pumpkyn's tfMin convention: the flown acceleration is Tnd/m(t).
% • `ndT` is the same conversion as a closure over thrust, for the rung
%   ladders that sweep thrust at fixed Isp and mass.
% • The defaults are the Earth-Moon canonical units every campaign here
%   uses; pass them explicitly for any other system.
%
%% Inputs:
%
%  thrustN                  double                  Thrust (N), positive
%                                                   scalar; [] for a RUNG
%                                                   LADDER, which has an Isp
%                                                   and a mass but no single
%                                                   thrust -- c and ndT are
%                                                   still defined, Tnd is NaN
%  ispS                     double                  Specific impulse (s);
%                                                   [] when the caller
%                                                   already HAS c (e.g.
%                                                   arclength_thrust takes
%                                                   it as an argument) --
%                                                   cnd is then NaN and ndC
%                                                   still exact
%  m0kg                     double                  Initial (wet) mass (kg)
%  lStar                    double (optional)       Characteristic length
%                                                   (km) [389703.264829278]
%  tStar                    double (optional)       Characteristic time (s)
%                                                   [382981.289129055]
%
%% Outputs:
%
%  p                        struct                  .g0 ND sea-level gravity
%                                                   .cnd ND exhaust speed
%                                                   .Tnd ND thrust accel at
%                                                   m = 1 .ndT closure
%                                                   (thrust N -> Tnd) .ndC
%                                                   closure (Isp s -> c_nd),
%                                                   for a continuation over
%                                                   Isp
%                                                   .lStar .tStar .ispS
%                                                   .m0kg .thrustN (echoed,
%                                                   so a consumer can carry
%                                                   the identity)
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
             pD = nd_propulsion(0.070, 900, 150);
    fprintf('70 mN, Isp 900 s, 150 kg -> T_nd = %.6e, c_nd = %.4f (c = %.4f km/s)\n', ...
            pD.Tnd, pD.cnd, pD.cnd*pD.lStar/pD.tStar);
              p = pD;
    return;
end

if nargin < 4 || isempty(lStar), lStar = 389703.264829278; end
if nargin < 5 || isempty(tStar), tStar = 382981.289129055; end

%% Refuse malformed input by name -- a vector thrust would silently build a
%% vector Tnd and propagate a shape error into a solver:
       noThrust = isempty(thrustN);          % a rung ladder: see Inputs
          noIsp = isempty(ispS);             % c supplied by the caller
              v = {m0kg, lStar, tStar};
              n = {'m0kg', 'lStar', 'tStar'};
if ~noThrust, v = [{thrustN} v];  n = [{'thrustN'} n]; end
if ~noIsp,    v = [{ispS}    v];  n = [{'ispS'}    n]; end
for k = 1:numel(v)
    if ~isscalar(v{k})
        error('nd_propulsion:scalar', '%s must be a scalar (use the ndT closure for a thrust ladder)', n{k});
    end
    if ~(isnumeric(v{k}) && isreal(v{k}) && isfinite(v{k}) && v{k} > 0)
        error('nd_propulsion:positive', '%s must be a positive finite real, got %g', n{k}, v{k});
    end
end

%% The three expressions, exactly as every call site wrote them:
           p.g0 = 9.80665*tStar^2/(1000*lStar);
if noIsp
          p.cnd = NaN;                         % undefined without an Isp
else
          p.cnd = (ispS/tStar)*p.g0;
end
if noThrust
          p.Tnd = NaN;                         % undefined without a thrust
else
          p.Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);
end
          p.ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);
          p.ndC = @(ispX) (ispX/tStar)*p.g0;

        p.lStar = lStar;     p.tStar = tStar;
         p.m0kg = m0kg;       p.ispS = NaN;        p.thrustN = NaN;
if ~noIsp,    p.ispS    = ispS;    end
if ~noThrust, p.thrustN = thrustN; end
end
