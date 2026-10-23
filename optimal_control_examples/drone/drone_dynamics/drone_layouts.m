function [x0, xT, geom] = drone_layouts(which, M, r)
%% Purpose:
%
%   Initial and final positions for the multi-agent problems, as 3 x M
%   arrays (2D layouts carry z = 0), with the geometry struct the clearance
%   and barrier functions read.
%
%   'paper'   Zhang, Liu, Darbon & Karniadakis (arXiv:2610.30388) Table 1
%             by M: M <= 16 one circle of radius 0.5; M = 32 two rings of
%             16 at radii [1/3, 0.5]; M = 64 two rings of 32 at
%             [0.33, 0.5]; M = 128 four rings of 32 at [0.35 0.40 0.45
%             0.50]. Targets antipodal. r = 0.02 for M <= 32, else 0.01.
%             The inner radii are INFERRED from the paper's own J_LB
%             column (1.156 at M = 32 needs 1/3; 2.297 at M = 64 needs
%             0.33); its text says "0.33" for both.
%   'circle'  one circle of radius 0.5, any M
%   'rings'   two rings of M/2 at radii [1/3, 0.5], z = 0 (M must be even)
%   'sphere'  OUR 3D layout: M points on a Fibonacci sphere of radius 0.4,
%             antipodal targets, ball obstacle at the origin, R = 0.1
%   'headon'  the two-agent closed-form oracle: (-0.4,0,0) <-> (0.4,0,0)
%
%% Inputs:
%
%  which                    char                    Layout name
%
%  M                        double                  Agents (ignored for
%                                                   'headon')
%
%  r                        double (optional)       Agent radius; default
%                                                   per layout
%
%% Outputs:
%
%  x0                       [3 x M]                 Initial positions
%
%  xT                       [3 x M]                 Final positions
%
%  geom                     struct                  .M .r .d (2|3) .which
%                                                   .eps [1e-5] .l [1e-4]
%                                                   .ball ([] or .c .R)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 8 circle:
     [a, b] = drone_layouts('paper', 8);
     figure('color', [1 1 1]);  plot(a(1,:), a(2,:), 'o', b(1,:), b(2,:), 'x');  axis equal;  grid on
     return
end
if nargin < 3, r = []; end
if nargin < 2, M = 2; end

d = 2;  ball = [];
switch lower(which)
    case 'paper'
        if M <= 16,      radii = 0.5;                 per = M;
        elseif M == 32,  radii = [1/3, 0.5];          per = 16;
        elseif M == 64,  radii = [0.33, 0.5];         per = 32;
        elseif M == 128, radii = [0.35 0.40 0.45 0.50]; per = 32;
        else, error('drone_layouts:M', 'the paper has no layout for M = %d', M);
        end
        if isempty(r), if M <= 32, r = 0.02; else, r = 0.01; end, end
        x0 = rings(radii, per);
    case 'circle'
        if isempty(r), r = 0.02; end
        x0 = rings(0.5, M);
    case 'rings'
        nRing = 2;
        if mod(M, nRing) ~= 0
            error('drone_layouts:M', 'M = %d is not divisible by %d rings', M, nRing);
        end
        if isempty(r), r = 0.02; end
        x0 = rings([1/3, 0.5], M/nRing);
    case 'sphere'
        if isempty(r), r = 0.02; end
        d = 3;
        ball = struct('c', [0; 0; 0], 'R', 0.1);
        x0 = fibonacci_sphere(M, 0.4);
    case 'headon'
        M = 2;  if isempty(r), r = 0.02; end
        x0 = [-0.4 0.4; 0 0; 0 0];
    otherwise
        error('drone_layouts:which', 'unknown layout ''%s''', which);
end
xT = -x0;
geom = struct('M', M, 'r', r, 'd', d, 'which', lower(which), 'eps', 1e-5, 'l', 1e-4, 'ball', ball);
end

% ------------------------------------------------------------------------
function X = rings(radii, per)
%% Purpose:
%
%   per points on each ring, equally spaced from angle 0, z = 0.
%
X = zeros(3, per*numel(radii));
ang = 2*pi*(0:per-1)/per;
for k = 1:numel(radii)
    cols = (k-1)*per + (1:per);
    X(:, cols) = [radii(k)*cos(ang); radii(k)*sin(ang); zeros(1, per)];
end
end

function X = fibonacci_sphere(M, radius)
%% Purpose:
%
%   M approximately equally spaced points on a sphere (golden-angle
%   spiral).
%
X = zeros(3, M);
gold = pi*(3 - sqrt(5));
for k = 1:M
    z = 1 - 2*(k - 0.5)/M;
    rho = sqrt(1 - z^2);
    X(:, k) = radius*[rho*cos(gold*(k-1)); rho*sin(gold*(k-1)); z];
end
end
