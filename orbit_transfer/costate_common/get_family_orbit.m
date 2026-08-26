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
%     'gto'    [.orientDeg, .sma_km, .ecc]   PSEUDO-family: an algebraic
%                                     (unpropagated) GTO departure locus at
%                                     fixed orientation, sampled uniform in
%                                     eccentric anomaly; tau is a LOCUS
%                                     PARAMETER (time since perigee over one
%                                     Kepler period), never an epoch offset
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
    case 'gto'
        % GTO pseudo-family (Stage B spec 2026-08-26): the locus of departure
        % states at FIXED orientation, parametrized by TIME fraction (mean
        % anomaly) over one Kepler period. Algebraic construction -- no
        % propagation; orientation is a sheet key, not a flow consequence.
        % orientDeg: Earth->Moon line to perigee direction, rotation sense.
        % tau is a LOCUS PARAMETER (time-since-perigee / period), NEVER an
        % epoch offset -- the pseudo-family carries no flight direction.
        smaKm  = fieldd(p, 'sma_km', (6378+350 + 6378+35786)/2);
        eccGto = fieldd(p, 'ecc', (35786-350)/(2*((6378+350+6378+35786)/2)));
        orient = fieldd(p, 'orientDeg', -25);
        lStar  = 389703.264829278;
        tStar  = 382981.289129055;
        muE = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
        % REVIEW FIX (GPT+Gemini, critical): sample uniform in ECCENTRIC
        % anomaly (dense) so perigee is resolved on an e~0.72 ellipse; tau
        % comes from Kepler's equation M(E). orb2eci's element 6 is TRUE
        % anomaly (per its own header, unconditionally -- the trailing flag
        % `dim3` only controls array-shape flattening, not the anomaly
        % convention), so E is converted to true anomaly analytically
        % (exact, no Newton solve needed in this forward direction) before
        % each call.
        M  = 2001;                                  % odd: sample lands on apogee
        E_ = linspace(0, 2*pi, M);
        mAnom  = E_ - eccGto*sin(E_);                % Kepler: M(E), monotone
        nuAnom = 2*atan2(sqrt(1+eccGto)*sin(E_/2), sqrt(1-eccGto)*cos(E_/2));
        Tkep   = 2*pi*sqrt((smaKm)^3/muE);           % s
        rvLoc  = zeros(M, 6);
        for km = 1:M
            [rr, vv] = pumpkyn.cr3bp.orb2eci(muE, ...
                         [smaKm, eccGto, 0, orient*pi/180, 0, nuAnom(km)], 2);
            rvLoc(km,:) = pumpkyn.cr3bp.fromPCI(0, [rr, vv], muStar, ...
                         tStar, lStar, 1);
        end
        tau  = (mAnom.'/(2*pi)) * (Tkep/tStar);
        rv   = rvLoc;
        info = struct('family', 'gto', 'periodND', tau(end), 'params', p, ...
                      'sma_km', smaKm, 'ecc', eccGto, 'orientDeg', orient);
        return;                                     % algebraic -- no cont_np/prop
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
