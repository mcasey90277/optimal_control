function u = pitchProgram(t,x,sched)
%% Purpose:
%
%  Evaluate a prescribed pitch-attitude boost guidance law at the current
%  time and state. Linear interpolation of the commanded pitch ATTITUDE
%  theta between grid points, held constant outside the grid -- identical
%  semantics to coorbital.guide.prescribed. Unlike prescribed, this law
%  reads the state: the commanded angle of attack is the difference
%  between the commanded attitude and the current flight-path angle -- the
%  classic prescribed-pitch boost law, where the vehicle is flown to an
%  attitude schedule and the angle of attack falls out of whatever gap opens
%  between that attitude and the velocity vector. If sched.alphaMax is
%  present, |alpha| is clamped to it (both bounds); if absent, NO clamping is
%  applied. Supply alphaMax whenever theta can run far ahead of gamma --
%  e.g. a 90 deg commanded attitude at liftoff, when gamma is still near
%  zero, would otherwise hand the equations of motion a physically absurd
%  angle of attack.
%
%% THE ATTITUDE RELATION IS NOT alpha = theta - gamma WHEN BANKED:
%
%  alpha is measured from the velocity vector in the plane set by the bank
%  angle sigma, the same convention lift and thrust already use in
%  coorbital.eom.boost3DOF. The thrust/body axis is therefore
%
%      bHat = cos(alpha)*vHat + sin(alpha)*(cos(sigma)*eGamma + sin(sigma)*ePsi)
%
%  and projecting it onto the local radial -- vHat.rHat = sin(gamma),
%  eGamma.rHat = cos(gamma), ePsi.rHat = 0 -- gives the pitch attitude as
%
%      sin(theta) = sin(gamma)*cos(alpha) + cos(gamma)*cos(sigma)*sin(alpha)
%                 = R*sin(alpha + phi)
%
%      R   = hypot(sin(gamma),cos(gamma)*cos(sigma))
%      phi = atan2(sin(gamma),cos(gamma)*cos(sigma))
%
%  so the angle of attack that ACHIEVES the commanded attitude is
%
%      alpha = asin(sin(theta)/R) - phi.
%
%  At sigma = 0 the thrust axis lies in the vertical plane, R is one, phi is
%  gamma, and this collapses exactly to the familiar alpha = theta - gamma;
%  that branch is written out separately below so the reduction is bit-exact
%  and every unbanked trajectory in this library is unchanged to its last
%  digit. Banked, only cos(sigma) of the incidence lifts the nose, so the
%  commanded attitude costs MORE angle of attack than the plain difference --
%  and beyond |sin(theta)| > R it cannot be bought at all, which raises
%  coorbital:pitchProgram:unreachableAttitude rather than commanding an
%  attitude the vehicle does not reach.
%
%  The principal arcsine branch is the one continuous with the unbanked law
%  over the documented |theta| <= 90 deg attitude range, which is where every
%  schedule in this library lives.
%
%% Inputs:
%
%  t                scalar                      Phase clock (s): the phase's
%                                               OWN tspan value, not rebased
%                                               to zero. A phase given
%                                               tspan = [10 50] is evaluated
%                                               at 10 through 50. See TWO
%                                               CLOCKS in
%                                               coorbital.prop.phaseRun
%
%  x                [7 x 1]                     State. Only x(5) = gamma
%                                               (rad) flight-path angle is
%                                               read.
%
%  sched            Struct                      Schedule:
%                                               tGrid [1 x K] (s)
%                                               theta [1 x K] (rad)
%                                                     commanded pitch
%                                                     attitude
%                                               sigma [1 x K] (rad)
%                                               alphaMax scalar (rad),
%                                                     OPTIONAL: clamps
%                                                     |alpha| to this value
%                                                     if the field is
%                                                     present; no clamping
%                                                     if absent
%
%% Outputs:
%
%  u                [2 x 1]                     [alpha; sigma] (rad). Raises
%                                               coorbital:pitchProgram:unreachableAttitude
%                                               when a banked geometry cannot
%                                               reach the commanded pitch
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Correct attitude relation at nonzero bank      08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
            schedD = struct('tGrid',[0 30 60],'theta',deg2rad([80 45 10]), ...
                            'sigma',deg2rad([0 -45 45]),'alphaMax',deg2rad(20));
                tD = linspace(-10,70,321)';
              xD = [0;0;0;0;deg2rad(10);0;0];
                uD = zeros(numel(tD),2);
    for kd = 1:numel(tD)
             uD(kd,:) = coorbital.guide.pitchProgram(tD(kd),xD,schedD).';
    end
    figure('color',[1 1 1]);
    plot(tD,rad2deg(uD),'linewidth',1.5); grid on;
    xlabel('time (s)'); ylabel('control (deg)');
    legend('angle of attack','bank angle','location','best');
    title(['Pitch program at gamma = 10 deg: alpha is theta(t) - gamma ' ...
           'unbanked, and more than that when banked']);
                 u = [];
    return;
end

%% Interpolate the commanded attitude and bank angle, clamping outside the
%  grid, exactly as coorbital.guide.prescribed does. Single-point schedules
%  are constant; interp1 needs two points:
if numel(sched.tGrid) == 1
             theta = sched.theta(1);
             sigma = sched.sigma(1);
else
             theta = interp1(sched.tGrid,sched.theta,t,'linear','extrap');
             sigma = interp1(sched.tGrid,sched.sigma,t,'linear','extrap');
    if t <= sched.tGrid(1)
             theta = sched.theta(1);
             sigma = sched.sigma(1);
    elseif t >= sched.tGrid(end)
             theta = sched.theta(end);
             sigma = sched.sigma(end);
    end
end

%% State-dependent angle of attack. See THE ATTITUDE RELATION in the header:
%  what the schedule commands is a pitch ATTITUDE, and what the equations of
%  motion accept is an incidence measured from the velocity vector in the
%  banked plane. The two coincide only at zero bank:
             gamma = x(5);
              cSig = cos(sigma);
if cSig == 1
%% Zero bank. R is exactly one and phi is exactly gamma, so the general
%% expression below collapses to the classic difference. It is written out
%% rather than routed through the general form because the two agree only to
%% about 1e-15 rad in floating point -- hypot(sin(gamma),cos(gamma)) is not
%% bit-exactly one, and atan2(sin(gamma),cos(gamma)) is not bit-exactly gamma
%% -- while every trajectory this library ships flies unbanked and is pinned
%% to its last printed digit. The two branches are the same expression
%% analytically; only the rounding differs:
             alpha = theta - gamma;
else
%% Banked: sin(theta) = R*sin(alpha + phi), inverted on the principal branch:
              rAmp = hypot(sin(gamma),cos(gamma).*cSig);
            phiRef = atan2(sin(gamma),cos(gamma).*cSig);
    if abs(sin(theta)) > rAmp
        error('coorbital:pitchProgram:unreachableAttitude', ...
            ['A commanded pitch attitude of %.6f deg cannot be reached at ' ...
             'gamma = %.6f deg with %.6f deg of bank. Banked, the thrust ' ...
             'axis is out of the vertical plane, so the radial component it ' ...
             'can produce tops out at %.9f while the commanded attitude ' ...
             'demands %.9f. Reduce the bank, reduce the commanded attitude, ' ...
             'or fly this segment unbanked.'], ...
            rad2deg(theta),rad2deg(gamma),rad2deg(sigma),rAmp,abs(sin(theta)));
    end
             alpha = asin(sin(theta)./rAmp) - phiRef;
end

%% Clamp |alpha| to alphaMax if the field is present; pass through unclamped
%  if it is absent:
if isfield(sched,'alphaMax')
             alpha = max(-sched.alphaMax,min(sched.alphaMax,alpha));
end

                 u = [alpha; sigma];
end
