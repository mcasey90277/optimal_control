function [c,ach] = terminalConstraint(rVec,vVec,des)
%% Purpose:
%
%  Measure a burnout state against the five terminal constraints that PEG and
%  VOA both enforce: radius, speed, flight-path angle, and the two components
%  that fix the orbital plane. Returns a DIMENSIONLESS residual vector that is
%  zero on the target manifold, so a shooting or least-squares solver can
%  drive it to zero without rescaling and neither algorithm has to re-derive
%  the five.
%
%  The mathematics is docs/hgv_dynamics_note.tex, Section "The five terminal
%  constraints", Eqs. (qbasis) and (termcon); it is not re-derived here. This
%  file evaluates it and guards it -- and the guards are the substance, not
%  the packaging, because two of them catch states the five constraints as
%  specified would score as a perfect hit.
%
%% THE NULL SET OF THE SPECIFIED FIVE HAS TWO SHEETS, AND THIS FILE SERVES ONE:
%
%  Residuals four and five ask that the achieved angular momentum have no
%  component on q1 and none on q2, which says h_f is PARALLEL TO THE LINE
%  through wHat -- not that it points the same way along it. A retrograde
%  burnout at the right radius, the right speed, the right flight-path angle
%  and in the right plane, flying round the wrong way, satisfies all five
%  EXACTLY:
%
%      c = [0  2.2e-16  0  0  0]',  and the plane angle is 180 degrees.
%
%  That is a defect in the specification, inherited from the source deck, and
%  this is the right place to catch it: both algorithms call this function in
%  their inner loop, and a dispersion campaign scoring runs on the norm of c
%  would record a wrong-way insertion as a perfect one. So the achieved
%  angular momentum is required to lie in the SAME half-space as wHat, and
%  h_f . wHat <= 0 raises. The branch is therefore selected HERE, not by the
%  caller, and a caller wanting the retrograde sheet must negate des.wHat.
%
%  A throw is preferred to a large residual, but BY ANALOGY rather than by
%  established contract, and the distinction matters because the solver that
%  will call this does not exist yet. The analogy is coorbital.util.aimSolve,
%  which treats a residual that throws as an infeasible point and shrinks the
%  step -- which is the behaviour wanted at the boundary between the sheets.
%  But aimSolve does NOT do that uniformly: it shrinks only on a TRIAL point,
%  hard-refuses with jacobianFailed on a finite-difference PROBE, and lets a
%  throw at the initial guess propagate uncaught. So there is no library-wide
%  rule here to appeal to. Whether a PEG least-squares or a VOA shooting
%  iteration should steer away from this throw or die on it is an OPEN DESIGN
%  QUESTION, and whichever builds first must decide it deliberately rather
%  than inherit it from this comment.
%
%% Inputs:
%
%  rVec             [3 x 1]                     Achieved burnout position (m).
%                                               A row vector is accepted and
%                                               reshaped
%
%  vVec             [3 x 1]                     Achieved burnout velocity
%                                               (m/s), same frame as rVec. A
%                                               row vector is accepted and
%                                               reshaped
%
%  des              Struct                      Desired burnout condition:
%                                               Rd     scalar (m), radius
%                                               Vd     scalar (m/s), speed
%                                               gammaD scalar (rad), flight
%                                                      path angle, positive up
%                                               wHat   [3 x 1], normal of the
%                                                      desired orbit plane,
%                                                      h/|h|; need not arrive
%                                                      normalised. ITS SIGN IS
%                                                      MEANINGFUL -- it picks
%                                                      the direction of
%                                                      travel, see the note
%                                                      above
%                                               bRef   [3 x 1], OPTIONAL roll
%                                                      reference fixing which
%                                                      way the in-plane basis
%                                                      points. It cancels from
%                                                      the constraint SET --
%                                                      it moves c(4) and c(5)
%                                                      individually and not
%                                                      the surface c = 0 they
%                                                      define, nor
%                                                      hypot(c(4),c(5)). Its
%                                                      LENGTH is irrelevant,
%                                                      only its direction.
%                                                      Absent or [] picks the
%                                                      coordinate axis least
%                                                      aligned with wHat
%
%% Outputs:
%
%  c                [5 x 1]                     Dimensionless residual, zero
%                                               on the target manifold:
%                                               c(1) radius error, r/Rd - 1
%                                               c(2) speed error, V/Vd - 1
%                                               c(3) flight-path-angle error
%                                               c(4) plane error on q1
%                                               c(5) plane error on q2
%
%  ach              Struct                      What the state actually is, at
%                                               full precision, so a caller
%                                               reporting a miss does not have
%                                               to recompute it:
%                                               rMag   (m)
%                                               vMag   (m/s)
%                                               gamma  (rad)
%                                               hMag   (m^2/s)
%                                               planeAngle (rad), angle
%                                                      between the achieved
%                                                      angular momentum and
%                                                      des.wHat; strictly less
%                                                      than pi/2 on any state
%                                                      this function accepts
%
%                                               Raises, rather than returning
%                                               a plausible number:
%                                               coorbital:terminalConstraint:badVector
%                                                      rVec or vVec is not
%                                                      three finite real
%                                                      components
%                                               coorbital:terminalConstraint:badTarget
%                                                      des is missing a field,
%                                                      or Rd, Vd or gammaD is
%                                                      not a usable scalar
%                                               coorbital:terminalConstraint:radialBurnout
%                                                      the DESIRED flight-path
%                                                      angle is radial, so the
%                                                      desired angular
%                                                      momentum is zero and
%                                                      the plane residuals
%                                                      have no scale
%                                               coorbital:terminalConstraint:badPlane
%                                                      des.wHat is malformed or
%                                                      cannot be normalised
%                                               coorbital:terminalConstraint:badRef
%                                                      des.bRef is malformed,
%                                                      cannot be normalised, or
%                                                      is parallel to des.wHat
%                                               coorbital:terminalConstraint:degenerateState
%                                                      the ACHIEVED position,
%                                                      speed or angular
%                                                      momentum is zero, which
%                                                      would leave three of
%                                                      the five residuals
%                                                      reading as satisfied
%                                               coorbital:terminalConstraint:retrograde
%                                                      the achieved angular
%                                                      momentum is in the
%                                                      opposite half-space to
%                                                      des.wHat, see the note
%                                                      above
%
%% Notes:
%
%  WHAT THIS IS NOT. These are orbital-INSERTION constraints. A ground target
%  is not among them, and mapping "hit this latitude and longitude" into
%  (Rd,Vd,gammaD,wHat) is a modelling decision the caller makes -- for a
%  ballistic boost, by the Keplerian free-flight range of
%  docs/hgv_dynamics_note.tex. This function takes the answer; it does not
%  make it.
%
%  THE FRAME IS THE CALLER'S, AND THAT IS DELIBERATE. Inputs are CARTESIAN
%  position and velocity, in whatever frame the caller's guidance problem is
%  posed in -- for PEG and VOA that is the inertial frame of the two-body
%  equations they integrate. This function never sees a library state and
%  never asks what the frame is; it only requires that rVec, vVec and des.wHat
%  are expressed in the SAME one. Converting the library's seven-state
%  spherical planet-relative vector into an inertial Cartesian pair is a
%  separate utility, because the velocity conversion v_inertial = v_rel +
%  omega x r belongs with whoever owns the chain's rotation rate, not with the
%  one function both guidance laws share.
%
%  ONE CONTINUOUS FREEDOM SURVIVES BY CONSTRUCTION. Rotating the whole
%  terminal state rigidly about des.wHat leaves all five residuals' content
%  unchanged -- c(1:3) exactly, and hypot(c(4),c(5)) exactly -- because the
%  five say nothing about WHERE in the target orbit the vehicle arrives. That
%  is not a defect; it is the freedom whose transversality condition VOA's
%  seventh boundary condition expresses. It is distinct from the discrete
%  two-sheet freedom above, which IS a defect and is refused.
%  tests/test_terminalConstraint.m asserts both.
%
%  THE JACOBIAN IS THE CALLER'S. None is returned and none is implied. PEG's
%  five-by-five solve will finite-difference this five times per iteration and
%  VOA's shooting seven times, which is affordable because the residual itself
%  costs no propagation. On the prograde sheet, away from the guards, c is
%  smooth in both vectors, so a central difference is well conditioned. At the
%  guards it throws rather than returning a large value, and a caller
%  differencing across a guard will see the throw, not a wrong derivative.
%
%% References:
%   [1] docs/hgv_dynamics_note.tex, Sections sec:termcon and sec:voa.
%   [2] Jaggers, R.F., "An Explicit Solution to the Exoatmospheric Powered
%       Flight Guidance and Trajectory Optimization Problem for Rocket
%       Propelled Vehicles," AIAA Guidance and Control Conference, 1977.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
           desD.Rd = c.rE + 180e3;
           desD.Vd = 7100;
       desD.gammaD = deg2rad(3);
         desD.wHat = [0; -sin(deg2rad(50)); cos(deg2rad(50))];
%% A state built to sit exactly on the manifold, then the same state with one
%% per cent more radius, then the same state flown the wrong way round:
              uHat = [1; 0; 0];
              sHat = cross(desD.wHat,uHat);
              sHat = sHat./sqrt(sHat.'*sHat);
               rOn = desD.Rd.*uHat;
               vOn = desD.Vd.*(sin(desD.gammaD).*uHat ...
                               + cos(desD.gammaD).*sHat);
           [cOn,~] = coorbital.guide.terminalConstraint(rOn,vOn,desD);
          [cOff,~] = coorbital.guide.terminalConstraint(1.01.*rOn,vOn,desD);
              vRet = desD.Vd.*(sin(desD.gammaD).*uHat ...
                               - cos(desD.gammaD).*sHat);
    fprintf('on  the manifold: c = [%9.2e %9.2e %9.2e %9.2e %9.2e]''\n',cOn);
    fprintf('one per cent high: c = [%9.2e %9.2e %9.2e %9.2e %9.2e]''\n',cOff);
    try
        coorbital.guide.terminalConstraint(rOn,vRet,desD);
    catch errD
        fprintf('retrograde: refused, %s\n',errD.identifier);
    end
    [c,ach] = deal([]);
    return;
end

%% Nothing below can be believed unless the inputs are numbers. A NaN slipping
%% into a residual makes a shooting solver wander rather than fail, which is
%% the expensive failure mode, so the vectors are checked first:
if ~isnumeric(rVec) || numel(rVec) ~= 3 || ~isreal(rVec) || ~all(isfinite(rVec))
    error('coorbital:terminalConstraint:badVector', ...
        'The position must be three finite real components; got %s.', ...
        mat2str(rVec(:).'));
end
if ~isnumeric(vVec) || numel(vVec) ~= 3 || ~isreal(vVec) || ~all(isfinite(vVec))
    error('coorbital:terminalConstraint:badVector', ...
        'The velocity must be three finite real components; got %s.', ...
        mat2str(vVec(:).'));
end
              rVec = rVec(:);
              vVec = vVec(:);

%% The desired condition. Rd and Vd divide, so they must be positive; gammaD
%% at plus or minus ninety degrees is a radial burnout, for which the angular
%% momentum vanishes and an orbit plane does not exist -- the normalisation
%% Rd*Vd*abs(cos(gammaD)) is exactly the desired |h|, and it is zero there:
         reqFields = {'Rd','Vd','gammaD','wHat'};
    for kf = 1:numel(reqFields)
        if ~isfield(des,reqFields{kf})
            error('coorbital:terminalConstraint:badTarget', ...
                'des is missing the field %s.',reqFields{kf});
        end
    end
                Rd = des.Rd;
                Vd = des.Vd;
            gammaD = des.gammaD;
if ~isscalar(Rd) || ~isreal(Rd) || ~isfinite(Rd) || Rd <= 0 || ...
   ~isscalar(Vd) || ~isreal(Vd) || ~isfinite(Vd) || Vd <= 0
    error('coorbital:terminalConstraint:badTarget', ...
        ['The desired radius and speed must be finite positive scalars; ' ...
         'got Rd = %s m and Vd = %s m/s.'],mat2str(Rd),mat2str(Vd));
end
if ~isscalar(gammaD) || ~isreal(gammaD) || ~isfinite(gammaD)
    error('coorbital:terminalConstraint:badTarget', ...
        'The desired flight-path angle must be a finite real scalar; got %s.', ...
        mat2str(gammaD));
end
if abs(cos(gammaD)) < 1e-8
    error('coorbital:terminalConstraint:radialBurnout', ...
        ['A desired flight-path angle of %.6f deg is radial: the desired ' ...
         'angular momentum Rd*Vd*cos(gammaD) is zero, so the two plane ' ...
         'residuals have no scale and the orbit plane is undefined.'], ...
        rad2deg(gammaD));
end

%% The desired plane normal, normalised here so a caller may hand over an
%% unnormalised angular momentum vector:
              wHat = des.wHat;
if ~isnumeric(wHat) || numel(wHat) ~= 3 || ~isreal(wHat) || ~all(isfinite(wHat))
    error('coorbital:terminalConstraint:badPlane', ...
        'des.wHat must be three finite real components; got %s.', ...
        mat2str(wHat(:).'));
end
              wHat = wHat(:);
              wMag = sqrt(wHat.'*wHat);
if wMag < 1e-12
    error('coorbital:terminalConstraint:badPlane', ...
        ['des.wHat has magnitude %.6g and cannot be normalised into a plane ' ...
         'normal.'],wMag);
end
              wHat = wHat./wMag;

%% The roll reference. It fixes which way q1 points inside the target plane
%% and cancels from the constraint set; absent, the coordinate axis least
%% aligned with wHat is the choice that is furthest from the degenerate one:
if ~isfield(des,'bRef') || isempty(des.bRef)
          [~,kMin] = min(abs(wHat));
              bRef = zeros(3,1);
        bRef(kMin) = 1;
else
              bRef = des.bRef;
    if ~isnumeric(bRef) || numel(bRef) ~= 3 || ~isreal(bRef) || ...
       ~all(isfinite(bRef))
        error('coorbital:terminalConstraint:badRef', ...
            'des.bRef must be three finite real components; got %s.', ...
            mat2str(bRef(:).'));
    end
              bRef = bRef(:);
end

%% NORMALISE THE REFERENCE BEFORE TESTING IT. Only its direction matters, so
%% the degeneracy test has to be on the SINE of the angle it makes with wHat.
%% An absolute threshold on the raw cross product would reject a short but
%% perfectly transverse bRef with a diagnosis -- "parallel to wHat" -- that is
%% the opposite of the truth:
              bMag = sqrt(bRef.'*bRef);
if bMag < 1e-300
    error('coorbital:terminalConstraint:badRef', ...
        'des.bRef has magnitude %.6g and has no direction.',bMag);
end
              bHat = bRef./bMag;
             q1Raw = cross(bHat,wHat);
            sinSep = sqrt(q1Raw.'*q1Raw);
if sinSep < 1e-8
    error('coorbital:terminalConstraint:badRef', ...
        ['des.bRef is parallel to des.wHat to within %.3g deg, so the ' ...
         'in-plane basis is undefined. Any bRef off the plane normal will ' ...
         'do; its length does not matter.'],rad2deg(asin(sinSep)));
end
                q1 = q1Raw./sinSep;
                q2 = cross(wHat,q1);

%% Achieved magnitudes. THE MAGNITUDES ARE GUARDED, not merely the
%% finiteness: a zero position vector is finite, real and three components
%% long, and it returns c = [-1  0  -sin(gammaD)  0  0] -- three of the five
%% residuals reading as satisfied, and no NaN anywhere in c to give it away.
%% A solver sees only c:
              rMag = sqrt(rVec.'*rVec);
              vMag = sqrt(vVec.'*vVec);
if rMag < eps(Rd) || vMag < eps(Vd)
    error('coorbital:terminalConstraint:degenerateState', ...
        ['The achieved state is degenerate: |r| = %.6g m and |v| = %.6g ' ...
         'm/s, against a desired %.6g m and %.6g m/s. The flight-path ' ...
         'angle is undefined and three of the five residuals would read as ' ...
         'satisfied.'],rMag,vMag,Rd,Vd);
end

%% Angular momentum, and the desired magnitude that makes the two plane
%% residuals dimensionless direction-cosine errors:
              hVec = cross(rVec,vVec);
              hMag = sqrt(hVec.'*hVec);
              hDes = Rd.*Vd.*abs(cos(gammaD));
if hMag <= eps(hDes)
    error('coorbital:terminalConstraint:degenerateState', ...
        ['The achieved velocity is radial: |r x v| = %.6g against a desired ' ...
         '%.6g, so the achieved plane is undefined and both plane residuals ' ...
         'would be zero for the wrong reason.'],hMag,hDes);
end

%% THE SECOND SHEET. See the note in the header: residuals four and five
%% cannot tell parallel from antiparallel, so the half-space is fixed here:
             hDotW = hVec.'*wHat;
if hDotW <= 0
    error('coorbital:terminalConstraint:retrograde', ...
        ['The achieved angular momentum is %.4f deg from des.wHat, in the ' ...
         'opposite half-space. The five constraints cannot see the ' ...
         'difference -- a retrograde burnout at the right radius, speed, ' ...
         'flight-path angle and plane satisfies all five exactly -- so it ' ...
         'is refused here. Negate des.wHat to target the retrograde ' ...
         'sheet.'],rad2deg(atan2(hMag,hDotW)));
end

%% The five, in the order of the source: radius, speed, flight-path angle,
%% then the two plane components:
             rDotV = rVec.'*vVec;
                 c = [ rMag./Rd - 1; ...
                       vMag./Vd - 1; ...
                       rDotV./(Rd.*Vd) - sin(gammaD); ...
                       (q1.'*hVec)./hDes; ...
                       (q2.'*hVec)./hDes ];

%% Diagnostics. Every quantity below is well defined, the guards above having
%% removed the cases in which it would not be:
          ach.rMag = rMag;
          ach.vMag = vMag;
         ach.gamma = asin(max(-1,min(1,rDotV./(rMag.*vMag))));
          ach.hMag = hMag;
    ach.planeAngle = atan2(sqrt(sum(cross(hVec,wHat).^2)),hDotW);
end
