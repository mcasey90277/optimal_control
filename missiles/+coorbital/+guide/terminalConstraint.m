function [c,ach] = terminalConstraint(rVec,vVec,des)
%% Purpose:
%
%  Measure a burnout state against the five terminal constraints that PEG and
%  VOA both enforce: radius, speed, flight-path angle, and the two components
%  that fix the orbital plane. Returns a DIMENSIONLESS residual vector that is
%  zero exactly on the target manifold, so a shooting or least-squares solver
%  can drive it to zero without rescaling and neither algorithm has to
%  re-derive the five.
%
%  The mathematics is docs/hgv_dynamics_note.tex, Section "The five terminal
%  constraints", Eqs. (qbasis) and (termcon); it is not re-derived here. This
%  file evaluates it and guards it.
%
%% THE FRAME IS THE CALLER'S, AND THAT IS DELIBERATE:
%
%  Inputs are CARTESIAN position and velocity, in whatever frame the caller's
%  guidance problem is posed in -- for PEG and VOA that is the inertial frame
%  of the two-body equations they integrate. This function never sees a
%  library state and never asks what the frame is; it only requires that
%  rVec, vVec and des.wHat are expressed in the SAME one.
%
%  It takes Cartesian vectors rather than the library's seven-state spherical
%  vector on purpose. Converting the library state to an inertial Cartesian
%  pair needs an inertial frame, and this library has not defined one: with
%  env.omegaE nonzero the state's longitude is Earth-fixed and its speed is
%  planet-relative, and there is no epoch anywhere in the library that says
%  where the inertial x-axis points. Baking a guess about that into the one
%  interface both guidance laws share would put the guess somewhere it could
%  not be seen. The conversion belongs with the caller that owns the epoch.
%
%% Inputs:
%
%  rVec             [3 x 1]                     Achieved burnout position (m)
%
%  vVec             [3 x 1]                     Achieved burnout velocity
%                                               (m/s), same frame as rVec
%
%  des              Struct                      Desired burnout condition:
%                                               Rd     scalar (m), radius
%                                               Vd     scalar (m/s), speed
%                                               gammaD scalar (rad), flight
%                                                      path angle, positive up
%                                               wHat   [3 x 1], unit normal of
%                                                      the desired orbit
%                                                      plane, h/|h|; need not
%                                                      arrive normalised
%                                               bRef   [3 x 1], OPTIONAL roll
%                                                      reference fixing which
%                                                      way the in-plane basis
%                                                      points. It cancels from
%                                                      the constraint SET --
%                                                      it moves c(4) and c(5)
%                                                      individually and not
%                                                      the surface c = 0 they
%                                                      define, nor
%                                                      hypot(c(4),c(5)).
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
%                                                      des.wHat
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
%  ONE ROTATIONAL FREEDOM SURVIVES BY CONSTRUCTION. Rotating the whole
%  terminal state rigidly about des.wHat leaves all five residuals' content
%  unchanged -- c(1:3) exactly, and hypot(c(4),c(5)) exactly -- because the
%  five say nothing about WHERE in the target orbit the vehicle arrives. That
%  is not a defect; it is the freedom whose transversality condition VOA's
%  seventh boundary condition expresses. tests/test_terminalConstraint.m
%  asserts it.
%
%% References:
%   [1] docs/hgv_dynamics_note.tex, Sections sec:termcon and sec:voa.
%   [2] Jaggers, R.F., "An Explicit Solution to the Exoatmospheric Powered
%       Flight Guidance and Trajectory Optimization Problem for Rocket
%       Propelled Vehicles," AIAA Guidance and Control Conference, 1977.
%
%% Revision History:
%  Michael Casey                                                08/09/2026
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
%% per cent more radius:
              uHat = [1; 0; 0];
              sHat = cross(desD.wHat,uHat);
              sHat = sHat./sqrt(sHat.'*sHat);
              rOn  = desD.Rd.*uHat;
              vOn  = desD.Vd.*(sin(desD.gammaD).*uHat + cos(desD.gammaD).*sHat);
           [cOn,~] = coorbital.guide.terminalConstraint(rOn,vOn,desD);
          [cOff,~] = coorbital.guide.terminalConstraint(1.01.*rOn,vOn,desD);
    fprintf('on  the manifold: c = [%9.2e %9.2e %9.2e %9.2e %9.2e]''\n',cOn);
    fprintf('one per cent high: c = [%9.2e %9.2e %9.2e %9.2e %9.2e]''\n',cOff);
    [c,ach] = deal([]);
    return;
end

%% Nothing below can be believed unless the inputs are numbers. A NaN slipping
%% into a residual makes a shooting solver wander rather than fail, which is
%% the expensive failure mode, so the vectors are checked first:
if ~isnumeric(rVec) || numel(rVec) ~= 3 || ~isreal(rVec) || ~all(isfinite(rVec))
    error('coorbital:terminalConstraint:badVector', ...
        ['The position must be three finite real components; got %s.'], ...
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

%% The in-plane orthonormal pair. A bRef parallel to wHat leaves nothing to
%% cross with and would silently produce a q1 of arbitrary direction:
             q1Raw = cross(bRef,wHat);
             q1Mag = sqrt(q1Raw.'*q1Raw);
if q1Mag < 1e-8.*max(1,sqrt(bRef.'*bRef))
    error('coorbital:terminalConstraint:badRef', ...
        ['des.bRef is parallel to des.wHat to within %.3g, so the in-plane ' ...
         'basis is undefined. Any bRef off the plane normal will do.'],q1Mag);
end
                q1 = q1Raw./q1Mag;
                q2 = cross(wHat,q1);

%% Achieved quantities:
              rMag = sqrt(rVec.'*rVec);
              vMag = sqrt(vVec.'*vVec);
              hVec = cross(rVec,vVec);
              hMag = sqrt(hVec.'*hVec);
              rDotV = rVec.'*vVec;

%% The desired angular momentum magnitude, which is what makes the two plane
%% residuals dimensionless direction-cosine errors:
              hDes = Rd.*Vd.*abs(cos(gammaD));

%% The five, in the order of the source: radius, speed, flight-path angle,
%% then the two plane components:
                 c = [ rMag./Rd - 1; ...
                       vMag./Vd - 1; ...
                       rDotV./(Rd.*Vd) - sin(gammaD); ...
                       (q1.'*hVec)./hDes; ...
                       (q2.'*hVec)./hDes ];

%% Diagnostics. gamma and the plane angle are undefined at a zero vector, so
%% they are reported as NaN rather than as a division:
    if rMag > 0 && vMag > 0
             gammaF = asin(max(-1,min(1,rDotV./(rMag.*vMag))));
    else
             gammaF = NaN;
    end
    if hMag > 0
              hCross = cross(hVec,wHat);
          planeAngle = atan2(sqrt(hCross.'*hCross),hVec.'*wHat);
    else
          planeAngle = NaN;
    end
          ach.rMag = rMag;
          ach.vMag = vMag;
         ach.gamma = gammaF;
          ach.hMag = hMag;
    ach.planeAngle = planeAngle;
end
