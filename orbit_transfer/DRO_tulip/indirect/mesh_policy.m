function [N, K, useSundman] = mesh_policy(revs, meshMode)
%% Purpose:
%
%   THE single place where the discretization is chosen. Given the expected
%   size of a transfer -- measured in DRO revolutions, the quantity that
%   actually controls difficulty -- returns the collocation mesh size N, the
%   multiple-shooting segment count K, and whether to use the Sundman-type
%   mesh.
%
%   The rules, and where they come from:
%
%     N  ~ 250 nodes per revolution (floor 400, cap 1200). The v2 thrust
%          ladder solved sub-revolution transfers cleanly at N = 400; the v1
%          low-thrust library used N = 800 over 3.6-7.5 revolutions. Both
%          are consistent with ~100-250 nodes per revolution, and the flown-
%          control verification gate is the real accuracy check either way.
%
%     K  ~ 10 multiple-shooting segments per revolution (floor 12, cap 80).
%          Short segments are what defeat the ~1000x error amplification of
%          single shooting; more revolutions need more segments.
%
%     Sundman from 1.5 revolutions up. The v2 ladder needed no Sundman mesh
%          through 1.15 revolutions; the v1 library used it at 3.6+. The
%          crossover is untested in between -- which is exactly what the
%          low-thrust pilot measures by running both modes.
%
%% Inputs:
%
%  revs                     double                  Expected transfer length
%                                                   in DRO revolutions
%                                                   (= t_f in ND time, since
%                                                   the tau = 1 DRO has unit
%                                                   period)
%
%  meshMode                 char                    'auto'    use the rules
%                                                   'uniform' force Sundman
%                                                             OFF
%                                                   'sundman' force Sundman
%                                                             ON
%
%% Outputs:
%
%  N                        double                  Collocation intervals
%
%  K                        double                  Multiple-shooting
%                                                   segments
%
%  useSundman               logical                 Sundman-type mesh on?
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, meshMode = 'auto'; end

       N = min(max(400, round(250*revs)), 1200);
       K = min(max(12,  round(10*revs)),  80);

switch lower(meshMode)
    case 'uniform',  useSundman = false;
    case 'sundman',  useSundman = true;
    otherwise,       useSundman = revs >= 1.5;
end

end
