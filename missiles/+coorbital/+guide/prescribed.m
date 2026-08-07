function u = prescribed(t,x,sched)
%% Purpose:
%
%  Evaluate a prescribed control schedule at the current time. Linear
%  interpolation between grid points, held constant outside the grid.
%
%% Inputs:
%
%  t                scalar                      Time since phase start (s)
%
%  x                [6 x 1]                     State. Unused; present so the
%                                               signature matches closed-loop
%                                               guidance laws.
%
%  sched            Struct                      Schedule:
%                                               tGrid [1 x K] (s)
%                                               alpha [1 x K] (rad)
%                                               sigma [1 x K] (rad)
%
%% Outputs:
%
%  u                [2 x 1]                     [alpha; sigma] (rad)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
            schedD = struct('tGrid',[0 30 60],'alpha',deg2rad([5 12 8]), ...
                            'sigma',deg2rad([0 -45 45]));
                tD = linspace(-10,70,321)';
                uD = zeros(numel(tD),2);
    for kd = 1:numel(tD)
             uD(kd,:) = coorbital.guide.prescribed(tD(kd),[],schedD).';
    end
    figure('color',[1 1 1]);
    plot(tD,rad2deg(uD),'linewidth',1.5); grid on;
    xlabel('time (s)'); ylabel('control (deg)');
    legend('angle of attack','bank angle','location','best');
    title('Prescribed schedule: linear inside the grid, clamped outside');
                 u = [];
    return;
end

%% Single-point schedules are constant; interp1 needs two points:
if numel(sched.tGrid) == 1
                 u = [sched.alpha(1); sched.sigma(1)];
    return;
end

%% Interpolate, clamping outside the grid:
             alpha = interp1(sched.tGrid,sched.alpha,t,'linear','extrap');
             sigma = interp1(sched.tGrid,sched.sigma,t,'linear','extrap');
    if t <= sched.tGrid(1)
             alpha = sched.alpha(1);
             sigma = sched.sigma(1);
    elseif t >= sched.tGrid(end)
             alpha = sched.alpha(end);
             sigma = sched.sigma(end);
    end
                 u = [alpha; sigma];
end
