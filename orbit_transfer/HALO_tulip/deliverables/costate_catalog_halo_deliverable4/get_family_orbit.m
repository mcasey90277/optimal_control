function [tau, rv, info] = get_family_orbit(family, p)
%% Purpose:
%
%   THE one place a costate campaign turns a family name + parameters into a
%   propagated periodic orbit. Every engine (thrust ladders, densifiers,
%   surveys, pickers' examples) builds its endpoints through this helper, so
%   adding a family here makes it available to the whole pipeline -- this is
%   what generalized the DRO->tulip machinery to Halo and beyond.
%
%   All orbits come from pumpkyn / pumpkynPie catalog getters, refined by
%   cont_np and propagated one period by pumpkyn.cr3bp.prop -- the same
%   construction the pumpkynPie demos use. Nothing is written into pumpkyn.
%
%   FAMILIES (p fields in brackets):
%     'dro'    [.tau]                 tau IS the period (ND)
%     'tulip'  [.Np, .pm]             period locked to 2*pi*(Np-2)/(Np-1);
%                                     pm = -1 branch (+1 is its z-mirror)
%     'halo'   [.tau, .Lpt, .pm]      period tau (ND), Lagrange point 1-3,
%                                     pm = +1 north / -1 south
%     'dpo'    [.tau]                 distant prograde orbit, period tau
%     'lyapunov' [.tau, .Lpt]         planar Lyapunov about L1-L3
%
%% Inputs:
%
%  family                   char                    Family name (above)
%
%  p                        struct                  Family parameters (above)
%                                                   plus optional .muStar
%                                                   [0.012150585609624] and
%                                                   .contTol [1e-12]
%
%% Outputs:
%
%  tau                      [M x 1]                 Time along one period
%                                                   (ND)
%
%  rv                       [M x 6]                 Rotating-frame states
%                                                   over one period
%
%  info                     struct                  .periodND, .family, and
%                                                   the resolved parameters
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: one of each family, overplotted:
     figure('Color','w'); hold on
     specs = { 'dro',   struct('tau',1.0), 'k';
               'tulip', struct('Np',7,'pm',-1), 'r';
               'halo',  struct('tau',(2*2*pi)/9,'Lpt',2,'pm',-1), 'b';
               'dpo',   struct('tau',1.0), 'm' };
     for k = 1:size(specs,1)
         try
             [~, x] = get_family_orbit(specs{k,1}, specs{k,2});
             plot3(x(:,1), x(:,2), x(:,3), ['-' specs{k,3}]);
         catch e
             fprintf('%s: %s\n', specs{k,1}, e.message);
         end
     end
     axis equal; grid on; legend(specs(:,1)); view(-35,25);
     return;
end

muStar  = fieldd(p, 'muStar', 0.012150585609624);
contTol = fieldd(p, 'contTol', 1e-12);

switch lower(family)
    case 'dro'
        tau0 = p.tau;
        [~, rv0] = pumpkynPie.cr3bp.getDRO(tau0);
    case 'tulip'
        tau0 = 2*pi*(p.Np-2)/(p.Np-1);
        [~, rv0] = pumpkyn.cr3bp.getTulip(tau0, p.Np, p.pm);
    case 'halo'
        [tau0, rv0] = pumpkynPie.cr3bp.getHalo(p.tau, p.Lpt, p.pm);
    case 'dpo'
        [tau0, rv0] = pumpkynPie.cr3bp.getDPO(p.tau);
    case 'lyapunov'
        [tau0, rv0] = pumpkynPie.cr3bp.getLyapunov(p.tau, p.Lpt);
    otherwise
        error('get_family_orbit:family', 'unknown family ''%s''', family);
end

rv0 = pumpkyn.cr3bp.cont_np(rv0, tau0, muStar, contTol);
[tau, rv] = pumpkyn.cr3bp.prop(tau0, rv0, muStar);
info = struct('family', lower(family), 'periodND', tau(end), 'params', p);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
