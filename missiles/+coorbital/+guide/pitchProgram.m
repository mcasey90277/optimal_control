function u = pitchProgram(t,x,sched)
%% Purpose:
%
%  Evaluate a prescribed pitch-attitude boost guidance law at the current
%  time and state. Linear interpolation of the commanded pitch ATTITUDE
%  theta between grid points, held constant outside the grid -- identical
%  semantics to coorbital.guide.prescribed. Unlike prescribed, this law
%  reads the state: the commanded angle of attack is the difference
%  between the commanded attitude and the current flight-path angle,
%  alpha = theta(t) - gamma, gamma = x(5) -- the classic prescribed-pitch
%  boost law, where the vehicle is flown to an attitude schedule and the
%  angle of attack falls out of whatever gap opens between that attitude
%  and the velocity vector. If sched.alphaMax is present, |alpha| is
%  clamped to it (both bounds); if absent, NO clamping is applied. Supply
%  alphaMax whenever theta can run far ahead of gamma -- e.g. a 90 deg
%  commanded attitude at liftoff, when gamma is still near zero, would
%  otherwise hand the equations of motion a physically absurd angle of
%  attack.
%
%% Inputs:
%
%  t                scalar                      Time since phase start (s)
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
%  u                [2 x 1]                     [alpha; sigma] (rad)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
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
    title('Pitch program: alpha = theta(t) - gamma, gamma = 10 deg fixed');
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

%% State-dependent angle of attack -- the classic prescribed-pitch law:
             gamma = x(5);
             alpha = theta - gamma;

%% Clamp |alpha| to alphaMax if the field is present; pass through unclamped
%  if it is absent:
if isfield(sched,'alphaMax')
             alpha = max(-sched.alphaMax,min(sched.alphaMax,alpha));
end

                 u = [alpha; sigma];
end
